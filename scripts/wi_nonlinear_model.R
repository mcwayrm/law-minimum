# ============================================================
# Wisconsin County-Level Nonlinear Soil Phosphorus Model
#
# Model:
#   P = (a0 + a1*X1) * [1 - (b0 + b1*X1) * exp(-exp(c0+c1*X1) * X2)]
#
#   y  = P_mean       (mean soil phosphorus, ppm)
#   X1 = OM_mean      (mean organic matter, %)
#   X2 = N_scaled     (farm N fertilizer, millions of kg)
#
# Note: c is reparametrized as exp(c0 + c1*X1) to guarantee the
#       decay rate is strictly positive. c0 and c1 are estimated
#       on the log scale.
#
# Data:
#   Soil: WI-DATCP-soil-summary-combined-1995-2019.xlsx
#   N:    N-P_from_fertilizer_1950-2017-july23-2020.xlsx
#         (census years 1987-2017, linearly interpolated to annual)
#
# Author: Generated via Claude
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(minpack.lm)    # Levenberg-Marquardt NLS

# ── 1. Load Soil Data ──────────────────────────────────────────────────────────
soil_raw <- read_excel("WI-DATCP-soil-summary-combined-1995-2019.xlsx")

om_data <- soil_raw |>
  filter(measure == "OM") |>
  select(county_fips, county, year, OM_mean = mean)

p_data <- soil_raw |>
  filter(measure == "P") |>
  select(county_fips, county, year, P_mean = mean)

soil_merged <- inner_join(om_data, p_data, by = c("county_fips", "county", "year"))

# ── 2. Load & Interpolate N Fertilizer Data ────────────────────────────────────
farm_raw <- read_excel("N-P_from_fertilizer_1950-2017-july23-2020.xlsx",
                       sheet = "farm")

wi_farm <- farm_raw |>
  filter(State == "WI") |>
  select(STCOFIPS, CountyName, matches("farmfertN-kg-")) |>
  pivot_longer(
    cols        = matches("farmfertN-kg-"),
    names_to    = "year_col",
    values_to   = "N_kg"
  ) |>
  mutate(year = as.integer(sub(".*-(\\d{4})$", "\\1", year_col))) |>
  select(county_fips = STCOFIPS, year, N_kg) |>
  arrange(county_fips, year)

# Linear interpolation from census years to annual 1995-2017
n_annual <- wi_farm |>
  group_by(county_fips) |>
  complete(year = 1987:2017) |>
  mutate(
    N_kg = approx(
      x    = year[!is.na(N_kg)],
      y    = N_kg[!is.na(N_kg)],
      xout = year
    )$y
  ) |>
  filter(year >= 1995, year <= 2017) |>
  ungroup() |>
  select(county_fips, year, N_kg)

# ── 3. Merge & Clean ──────────────────────────────────────────────────────────
panel <- soil_merged |>
  inner_join(n_annual, by = c("county_fips", "year")) |>
  filter(!is.na(OM_mean), !is.na(N_kg), !is.na(P_mean),
         OM_mean > 0, N_kg > 0, P_mean > 0) |>
  mutate(N_scaled = N_kg / 1e6)   # millions of kg for numerical stability

cat(sprintf("Panel: %d obs | %d counties | %d–%d\n",
            nrow(panel), n_distinct(panel$county),
            min(panel$year), max(panel$year)))

# ── 4. Nonlinear Model Function ───────────────────────────────────────────────
# c is modelled as exp(c0 + c1*X1) to keep decay rate strictly > 0
model_fn <- function(X1, X2, a0, a1, b0, b1, c0, c1) {
  A <- a0 + a1 * X1
  B <- b0 + b1 * X1
  C <- exp(c0 + c1 * X1)          # reparametrized decay
  A * (1 - B * exp(-C * X2))
}

# ── 5. Estimation via Levenberg-Marquardt ─────────────────────────────────────
# Starting values derived from global search (differential evolution in Python)
start_vals <- list(
  a0 =  32.1,
  a1 =   3.85,
  b0 =  -2.02,
  b1 =   0.60,
  c0 =  -8.31,
  c1 =   2.20
)

fit <- nlsLM(
  P_mean ~ model_fn(OM_mean, N_scaled, a0, a1, b0, b1, c0, c1),
  data    = panel,
  start   = start_vals,
  control = nls.lm.control(maxiter = 2000, ftol = 1e-12, ptol = 1e-12)
)

# ── 6. Parameter Table ────────────────────────────────────────────────────────
s          <- summary(fit)
params     <- coef(fit)
se_vals    <- s$coefficients[, "Std. Error"]
t_vals     <- s$coefficients[, "t value"]
p_vals     <- s$coefficients[, "Pr(>|t|)"]
ci_lo      <- params - 1.96 * se_vals
ci_hi      <- params + 1.96 * se_vals

param_table <- data.frame(
  Parameter   = c("a0", "a1", "b0", "b1", "c0 (log-scale)", "c1 (log-scale)"),
  Estimate    = round(params,   5),
  Std_Error   = round(se_vals,  5),
  t_value     = round(t_vals,   3),
  p_value     = round(p_vals,   4),
  CI_95_Low   = round(ci_lo,    5),
  CI_95_High  = round(ci_hi,    5),
  row.names   = NULL
)

cat("\n=== Parameter Estimates ===\n")
print(param_table)

y_hat    <- fitted(fit)
resid_v  <- residuals(fit)
ss_res   <- sum(resid_v^2)
ss_tot   <- sum((panel$P_mean - mean(panel$P_mean))^2)
r2       <- 1 - ss_res / ss_tot
rmse     <- sqrt(mean(resid_v^2))

cat(sprintf("\nPseudo R²   = %.4f\nRMSE        = %.4f ppm\nResidual SE = %.4f ppm\nn obs       = %d\n",
            r2, rmse, s$sigma, nrow(panel)))

# Effective decay at representative OM values
cat("\nEffective decay c = exp(c0 + c1*OM) at key OM values:\n")
for (om_val in c(2, 3, 4, 5)) {
  c_eff <- exp(params["c0"] + params["c1"] * om_val)
  cat(sprintf("  OM = %.1f  →  c_eff = %.4f\n", om_val, c_eff))
}

panel <- panel |>
  mutate(P_fitted = y_hat, resid = resid_v)

# ── 7. Plots ──────────────────────────────────────────────────────────────────
BLUE  <- "#2C7BB6"
RED   <- "#D73027"
LBLUE <- "#91BFDB"

## A — Fitted vs Observed
p1 <- ggplot(panel, aes(x = P_mean, y = P_fitted)) +
  geom_point(alpha = 0.25, size = 1.3, color = BLUE) +
  geom_abline(slope = 1, intercept = 0, color = RED, linetype = "dashed", linewidth = 1.1) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("Pseudo R² = %.3f\nRMSE = %.1f ppm\nn = %d", r2, rmse, nrow(panel)),
           hjust = -0.1, vjust = 1.3, size = 3.7, color = "#333333") +
  labs(title    = "A — Fitted vs Observed Soil P",
       x = "Observed Soil P (ppm)", y = "Fitted Soil P (ppm)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.background = element_rect(fill = "#f9f9f9"))

## B — Response over OM at three N levels
N_quants <- quantile(panel$N_scaled, c(0.10, 0.50, 0.90))
OM_seq   <- seq(quantile(panel$OM_mean, 0.02), quantile(panel$OM_mean, 0.98), length.out = 300)

surface_b <- expand.grid(OM_mean = OM_seq,
                          N_label = c("Low N (10th pct)","Median N (50th pct)","High N (90th pct)")) |>
  mutate(N_scaled = dplyr::case_when(
    N_label == "Low N (10th pct)"  ~ N_quants[1],
    N_label == "Median N (50th pct)" ~ N_quants[2],
    N_label == "High N (90th pct)" ~ N_quants[3]
  ),
  P_fitted = model_fn(OM_mean, N_scaled, params["a0"], params["a1"],
                      params["b0"], params["b1"], params["c0"], params["c1"]),
  N_label = factor(N_label, levels = c("Low N (10th pct)","Median N (50th pct)","High N (90th pct)")))

p2 <- ggplot(surface_b, aes(x = OM_mean, y = P_fitted, color = N_label)) +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = c(LBLUE, BLUE, RED)) +
  labs(title = "B — Fitted P vs OM | Three N Levels",
       x = "Organic Matter — X₁ (%)", y = "Fitted Soil P (ppm)",
       color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "#f9f9f9"))

## C — Response over N at three OM levels
OM_quants <- quantile(panel$OM_mean, c(0.10, 0.50, 0.90))
N_seq     <- seq(quantile(panel$N_scaled, 0.02), quantile(panel$N_scaled, 0.98), length.out = 300)

surface_c <- expand.grid(N_scaled = N_seq,
                          OM_label = c("Low OM (10th pct)","Median OM (50th pct)","High OM (90th pct)")) |>
  mutate(OM_mean = dplyr::case_when(
    OM_label == "Low OM (10th pct)"    ~ OM_quants[1],
    OM_label == "Median OM (50th pct)" ~ OM_quants[2],
    OM_label == "High OM (90th pct)"   ~ OM_quants[3]
  ),
  P_fitted = model_fn(OM_mean, N_scaled, params["a0"], params["a1"],
                      params["b0"], params["b1"], params["c0"], params["c1"]),
  OM_label = factor(OM_label, levels = c("Low OM (10th pct)","Median OM (50th pct)","High OM (90th pct)")))

p3 <- ggplot(surface_c, aes(x = N_scaled, y = P_fitted, color = OM_label)) +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = c(LBLUE, BLUE, RED)) +
  labs(title = "C — Fitted P vs N Fertilizer | Three OM Levels",
       x = "N Fertilizer — X₂ (million kg)", y = "Fitted Soil P (ppm)",
       color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "#f9f9f9"))

## D — Mean residual by year
resid_yr <- panel |>
  group_by(year) |>
  summarise(mean_r = mean(resid), se_r = sd(resid) / sqrt(n()), .groups = "drop")

p4 <- ggplot(resid_yr, aes(x = year, y = mean_r)) +
  geom_ribbon(aes(ymin = mean_r - 1.96*se_r, ymax = mean_r + 1.96*se_r),
              fill = BLUE, alpha = 0.22) +
  geom_line(color = BLUE, linewidth = 1.4) +
  geom_point(color = BLUE, size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = RED, linewidth = 1.1) +
  labs(title    = "D — Mean Residual by Year",
       subtitle = "Shaded band = 95% CI around mean residual",
       x = "Year", y = "Mean Residual (ppm)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "#f9f9f9"))

combined <- grid.arrange(p1, p2, p3, p4, ncol = 2,
  top = grid::textGrob(
    "Wisconsin County Panel — Nonlinear Soil Phosphorus Model (1995–2017)\nP = (a\u2080+a\u2081·OM)[1−(b\u2080+b\u2081·OM)·exp(−exp(c\u2080+c\u2081·OM)·N)]",
    gp = grid::gpar(fontsize = 13, fontface = "bold")))

ggsave("wi_nonlinear_model_results.png", combined,
       width = 14, height = 10, dpi = 150, bg = "white")

write.csv(param_table, "wi_model_parameter_table.csv", row.names = FALSE)
cat("\nFiles saved: wi_nonlinear_model_results.png, wi_model_parameter_table.csv\n")
