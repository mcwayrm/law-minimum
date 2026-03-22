# ============================================================
# Wisconsin County-Level Corn Production CES Model
#
# Model (Constant Elasticity of Substitution):
#   y_{ct} = A_{ct} * [delta * K_{ct}^(-rho) + (1-delta) * N_{ct}^(-rho)]^(-1/rho)
#
#   A_{ct} = a0 + a1 * OM   (Hicks-neutral TFP shifter)
#   K_{ct} = b0 + b1 * OM   (effective capital/OM composite)
#
# Reparameterizations for identification:
#   delta  in (0,1) via  delta = 1 / (1 + exp(-d))    [logistic]
#   rho > -1 enforced by bounding during estimation
#
#   y  = corn_grain_yield  (county corn grain yield, bu/acre)
#   X1 = OM_mean           (mean soil organic matter, %)
#   X2 = N_rate_scaled     (N fertilizer per planted corn acre, / 100)
#
# Data:
#   Soil: WI-DATCP-soil-summary-combined-1995-2019.xlsx
#   N:    N-P_from_fertilizer_1950-2017-july23-2020.xlsx
#         (census years 1987-2017, linearly interpolated to annual)
#   Yield: USDA NASS Quick Stats
#
# Author: Ryan McWay
# ============================================================

library(readxl)
library(readr)
library(rnassqs)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(minpack.lm)    # Levenberg-Marquardt NLS

set.seed(123)

# ── 1. Load Soil Data ──────────────────────────────────────────────────────────
setwd("E:/law-minimum/")
soil_raw <- read_excel("./data/WI-DATCP-soil-summary-combined-1995-2019.xlsx")

om_data <- soil_raw |>
  filter(measure == "OM") |>
  select(county_fips, county, year, OM_mean = `mean`) |>
  mutate(year        = as.character(year),
         county_fips = as.character(county_fips))

# ── 2. Load & Interpolate N Fertilizer Data ────────────────────────────────────
farm_raw <- read_excel("./data/N-P_from_fertilizer_1950-2017-july23-2020.xlsx",
                       sheet = "farm")

wi_farm <- farm_raw |>
  filter(State == "WI") |>
  select(STCOFIPS, CountyName, matches("farmfertN-kg-")) |>
  pivot_longer(
    cols      = matches("farmfertN-kg-"),
    names_to  = "year_col",
    values_to = "N_kg"
  ) |>
  mutate(year = as.integer(sub(".*-(\\d{4})$", "\\1", year_col))) |>
  select(county_fips = STCOFIPS, year, N_kg) |>
  mutate(county_fips = as.character(county_fips)) |>
  arrange(county_fips, year)

n_annual <- wi_farm |>
  group_by(county_fips) |>
  complete(year = 1987:2017) |>
  mutate(
    N_kg = approx(
      x    = as.numeric(year[!is.na(N_kg)]),
      y    = N_kg[!is.na(N_kg)],
      xout = as.numeric(year)
    )$y
  ) |>
  filter(year >= 1995, year <= 2017) |>
  ungroup() |>
  select(county_fips, year, N_kg) |>
  mutate(year = as.character(year))

# ── 3. Load Corn Yield & Acres Data ───────────────────────────────────────────
NASSQS_TOKEN <- trimws(readr::read_file("./scripts/nassqs-api-key.txt"))
nassqs_auth(key = NASSQS_TOKEN)

params_yield <- list(group_desc            = "FIELD CROPS",
                     short_desc            = "CORN, GRAIN - YIELD, MEASURED IN BU / ACRE",
                     agg_level_desc        = "COUNTY",
                     state_alpha           = "WI",
                     reference_period_desc = "YEAR")
df_wi_corn_yield <- setDT(nassqs(params_yield))
df_wi_corn_yield <- df_wi_corn_yield[, .(state_fips_code, county_code, county_name, year, Value)]
setnames(df_wi_corn_yield, "Value", "corn_grain_yield")
df_wi_corn_yield[, county_fips      := paste0(state_fips_code, county_code)]
df_wi_corn_yield[, year             := as.character(year)]
df_wi_corn_yield[, corn_grain_yield := as.numeric(gsub(",", "", corn_grain_yield))]

params_acres <- list(group_desc            = "FIELD CROPS",
                     short_desc            = "CORN - ACRES PLANTED",
                     agg_level_desc        = "COUNTY",
                     state_alpha           = "WI",
                     reference_period_desc = "YEAR")
df_wi_corn_acres <- setDT(nassqs(params_acres))
df_wi_corn_acres <- df_wi_corn_acres[, .(state_fips_code, county_code, county_name, year, Value)]
setnames(df_wi_corn_acres, "Value", "corn_acres")
df_wi_corn_acres[, county_fips := paste0(state_fips_code, county_code)]
df_wi_corn_acres[, year        := as.character(year)]
df_wi_corn_acres[, corn_acres  := as.numeric(gsub(",", "", corn_acres))]

# ── 4. Merge & Clean ──────────────────────────────────────────────────────────
panel <- om_data |>
  inner_join(n_annual,          by = c("county_fips", "year")) |>
  inner_join(df_wi_corn_yield,  by = c("county_fips", "year")) |>
  inner_join(df_wi_corn_acres,  by = c("county_fips", "year")) |>
  filter(!is.na(OM_mean), !is.na(N_kg), !is.na(corn_grain_yield), !is.na(corn_acres),
         OM_mean > 0, N_kg > 0, corn_grain_yield > 0, corn_acres > 0) |>
  mutate(
    N_rate_kg_acre = N_kg / corn_acres,
    N_rate_scaled  = N_rate_kg_acre / 100
  )

cat(sprintf("Panel: %d obs | %d counties | %s–%s\n",
            nrow(panel), n_distinct(panel$county),
            min(panel$year), max(panel$year)))

# ── 5. CES Model Function ──────────────────────────────────────────────────────
# y = A * [delta * K^(-rho) + (1-delta) * N^(-rho)]^(-1/rho)
# A    = a0 + a1*OM
# K    = b0 + b1*OM
# delta in (0,1) via logistic:  delta = 1 / (1 + exp(-d))
# rho is estimated directly; NLS bounds will keep rho > -1 and K,A > 0
ces_fn <- function(X1, X2, a0, a1, b0, b1, d, rho) {
  A     <- a0 + a1 * X1
  K     <- b0 + b1 * X1
  delta <- 1 / (1 + exp(-d))          # logistic reparameterization of share
  inside <- delta * K^(-rho) + (1 - delta) * X2^(-rho)
  A * inside^(-1 / rho)
}

# ── 6. Estimation via Levenberg-Marquardt ─────────────────────────────────────
# Starting values:
#   a0 ~ 100 (bu/acre baseline),  a1 small positive
#   b0 ~ 1 (OM composite ~1 at avg OM), b1 small
#   d  = 0  → delta = 0.5  (equal shares)
#   rho = 0.5  → sigma = 2  (moderate substitutability)
start_vals <- list(
  a0  = 100,
  a1  =  5,
  b0  =  1.0,
  b1  =  0.1,
  d   =  0.0,
  rho =  0.5
)

# Lower/upper bounds: ensure A > 0, K > 0, rho > -1 (rho = -1 is Leontief limit)
lower_bounds <- c(a0 = 0,    a1 = -Inf, b0 = 1e-4, b1 = -Inf, d = -10, rho = -0.99)
upper_bounds <- c(a0 = Inf,  a1 =  Inf, b0 =  Inf, b1 =  Inf, d =  10, rho =  10)

fit <- nlsLM(
  corn_grain_yield ~ ces_fn(OM_mean, N_rate_scaled, a0, a1, b0, b1, d, rho),
  data    = panel,
  start   = start_vals,
  lower   = lower_bounds,
  upper   = upper_bounds,
  control = nls.lm.control(maxiter = 2000, ftol = 1e-12, ptol = 1e-12)
)

# ── 7. Parameter Table ────────────────────────────────────────────────────────
s       <- summary(fit)
coefs   <- coef(fit)
se_vals <- s$coefficients[, "Std. Error"]
t_vals  <- s$coefficients[, "t value"]
p_vals  <- s$coefficients[, "Pr(>|t|)"]
ci_lo   <- coefs - 1.96 * se_vals
ci_hi   <- coefs + 1.96 * se_vals

# Derived quantities
delta_est <- 1 / (1 + exp(-coefs["d"]))
sigma_est <- 1 / (1 + coefs["rho"])

param_table <- data.frame(
  Parameter  = c("a0", "a1", "b0", "b1", "d (logit-delta)", "rho"),
  Estimate   = round(coefs,   5),
  Std_Error  = round(se_vals, 5),
  t_value    = round(t_vals,  3),
  p_value    = round(p_vals,  4),
  CI_95_Low  = round(ci_lo,   5),
  CI_95_High = round(ci_hi,   5),
  row.names  = NULL
)

cat("\n=== CES Parameter Estimates ===\n")
print(param_table)
cat(sprintf("\nDerived  delta = %.4f  (share of K in CES nest)\n", delta_est))
cat(sprintf("Derived  sigma = %.4f  (elasticity of substitution between K and N)\n", sigma_est))

y_hat   <- fitted(fit)
resid_v <- residuals(fit)
ss_res  <- sum(resid_v^2)
ss_tot  <- sum((panel$corn_grain_yield - mean(panel$corn_grain_yield))^2)
r2      <- 1 - ss_res / ss_tot
rmse    <- sqrt(mean(resid_v^2))

cat(sprintf("\nPseudo R²   = %.4f\nRMSE        = %.4f bu/acre\nResidual SE = %.4f bu/acre\nn obs       = %d\n",
            r2, rmse, s$sigma, nrow(panel)))

panel <- panel |>
  mutate(yield_fitted = y_hat, resid = resid_v)

# ── 8. Plots ──────────────────────────────────────────────────────────────────
BLUE  <- "#2C7BB6"
RED   <- "#D73027"
LBLUE <- "#91BFDB"

## A — Fitted vs Observed
p1 <- ggplot(panel, aes(x = corn_grain_yield, y = yield_fitted)) +
  geom_point(alpha = 0.25, size = 1.3, color = BLUE) +
  geom_abline(slope = 1, intercept = 0, color = RED, linetype = "dashed", linewidth = 1.1) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("Pseudo R\u00b2 = %.3f\nRMSE = %.1f bu/acre\nn = %d", r2, rmse, nrow(panel)),
           hjust = -0.1, vjust = 1.3, size = 3.7, color = "#333333") +
  labs(title = "A \u2014 Fitted vs Observed Corn Grain Yield (CES)",
       x = "Observed Corn Yield (bu/acre)", y = "Fitted Corn Yield (bu/acre)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.background = element_rect(fill = "#f9f9f9"))

## B — Response over OM at three N levels
N_quants <- quantile(panel$N_rate_scaled, c(0.10, 0.50, 0.90))
OM_seq   <- seq(quantile(panel$OM_mean, 0.02), quantile(panel$OM_mean, 0.98), length.out = 300)

surface_b <- expand.grid(
  OM_mean = OM_seq,
  N_label = c("Low N (10th pct)", "Median N (50th pct)", "High N (90th pct)")
) |>
  mutate(
    N_rate_scaled = dplyr::case_when(
      N_label == "Low N (10th pct)"    ~ N_quants[1],
      N_label == "Median N (50th pct)" ~ N_quants[2],
      N_label == "High N (90th pct)"   ~ N_quants[3]
    ),
    N_rate_kg_acre = 100 * N_rate_scaled,
    yield_fitted   = ces_fn(OM_mean, N_rate_scaled,
                            coefs["a0"], coefs["a1"], coefs["b0"], coefs["b1"],
                            coefs["d"],  coefs["rho"]),
    N_label = factor(N_label, levels = c("Low N (10th pct)", "Median N (50th pct)", "High N (90th pct)"))
  )

p2 <- ggplot(surface_b, aes(x = OM_mean, y = yield_fitted, color = N_label)) +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = c(LBLUE, BLUE, RED)) +
  labs(title = "B \u2014 Fitted Yield vs OM | Three N Levels (CES)",
       x = "Organic Matter (%)", y = "Fitted Corn Yield (bu/acre)", color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "#f9f9f9"))

## C — Response over N at three OM levels
OM_quants <- quantile(panel$OM_mean, c(0.10, 0.50, 0.90))
N_seq     <- seq(quantile(panel$N_rate_kg_acre, 0.02), quantile(panel$N_rate_kg_acre, 0.98), length.out = 300)

surface_c <- expand.grid(
  N_rate_kg_acre = N_seq,
  OM_label = c("Low OM (10th pct)", "Median OM (50th pct)", "High OM (90th pct)")
) |>
  mutate(
    OM_mean = dplyr::case_when(
      OM_label == "Low OM (10th pct)"    ~ OM_quants[1],
      OM_label == "Median OM (50th pct)" ~ OM_quants[2],
      OM_label == "High OM (90th pct)"   ~ OM_quants[3]
    ),
    N_rate_scaled = N_rate_kg_acre / 100,
    yield_fitted  = ces_fn(OM_mean, N_rate_scaled,
                           coefs["a0"], coefs["a1"], coefs["b0"], coefs["b1"],
                           coefs["d"],  coefs["rho"]),
    OM_label = factor(OM_label, levels = c("Low OM (10th pct)", "Median OM (50th pct)", "High OM (90th pct)"))
  )

p3 <- ggplot(surface_c, aes(x = N_rate_kg_acre, y = yield_fitted, color = OM_label)) +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = c(LBLUE, BLUE, RED)) +
  labs(title = "C \u2014 Fitted Yield vs N Fertilizer | Three OM Levels (CES)",
       x = "N Fertilizer (kg per planted corn acre)", y = "Fitted Corn Yield (bu/acre)", color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "#f9f9f9"))

## D — Mean residual by year
resid_yr <- panel |>
  mutate(year_num = as.integer(year)) |>
  group_by(year_num) |>
  summarise(mean_r = mean(resid), se_r = sd(resid) / sqrt(n()), .groups = "drop")

p4 <- ggplot(resid_yr, aes(x = year_num, y = mean_r)) +
  geom_ribbon(aes(x = year_num, ymin = mean_r - 1.96 * se_r, ymax = mean_r + 1.96 * se_r),
              fill = BLUE, alpha = 0.22, inherit.aes = FALSE) +
  geom_line(color = BLUE, linewidth = 1.4) +
  geom_point(color = BLUE, size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = RED, linewidth = 1.1) +
  scale_x_continuous(breaks = seq(1995, 2020, by = 5)) +
  labs(title    = "D \u2014 Mean Residual by Year (CES)",
       subtitle = "Shaded band = 95% CI around mean residual",
       x = "Year", y = "Mean Residual (bu/acre)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.background = element_rect(fill = "#f9f9f9"))

combined <- grid.arrange(p1, p2, p3, p4, ncol = 2,
  top = grid::textGrob(
    sprintf("Wisconsin CES Model (1995\u20132017)  |  delta=%.3f  rho=%.3f  sigma=%.3f",
            delta_est, coefs["rho"], sigma_est),
    gp = grid::gpar(fontsize = 12, fontface = "bold")))

ggsave("./output/figures/figure_wi_ces_plot_A.png",     p1, width = 7,  height = 5,  dpi = 150, bg = "white")
ggsave("./output/figures/figure_wi_ces_plot_B.png",     p2, width = 7,  height = 5,  dpi = 150, bg = "white")
ggsave("./output/figures/figure_wi_ces_plot_C.png",     p3, width = 7,  height = 5,  dpi = 150, bg = "white")
ggsave("./output/figures/figure_wi_ces_plot_D.png",     p4, width = 7,  height = 5,  dpi = 150, bg = "white")
ggsave("./output/figures/figure_wi_ces_visuals.png", combined, width = 14, height = 10, dpi = 150, bg = "white")

## E — Contour map
OM_contour_seq <- seq(min(panel$OM_mean, na.rm = TRUE), max(panel$OM_mean, na.rm = TRUE), length.out = 140)
N_contour_seq  <- seq(min(panel$N_rate_kg_acre, na.rm = TRUE), max(panel$N_rate_kg_acre, na.rm = TRUE), length.out = 140)

surface_contour <- expand.grid(OM_mean = OM_contour_seq, N_rate_kg_acre = N_contour_seq) |>
  mutate(
    N_rate_scaled = N_rate_kg_acre / 100,
    yield_fitted  = ces_fn(OM_mean, N_rate_scaled,
                           coefs["a0"], coefs["a1"], coefs["b0"], coefs["b1"],
                           coefs["d"],  coefs["rho"])
  )

p5 <- ggplot(surface_contour, aes(x = OM_mean, y = N_rate_kg_acre, fill = yield_fitted)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(data = surface_contour,
               aes(x = OM_mean, y = N_rate_kg_acre, z = yield_fitted),
               color = "white", alpha = 0.6, linewidth = 0.3, inherit.aes = FALSE) +
  geom_point(data = panel, aes(x = OM_mean, y = N_rate_kg_acre),
             inherit.aes = FALSE, color = "black", alpha = 0.20, size = 0.8) +
  scale_fill_gradient(low = BLUE, high = RED) +
  labs(title    = "E \u2014 CES Predicted Yield Contours by OM and N",
       subtitle = "Points show county-year observations",
       x = "Organic Matter (%)", y = "N Fertilizer (kg per planted corn acre)",
       fill = "Fitted\nYield") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.background = element_rect(fill = "#f9f9f9"))

ggsave("./output/figures/figure_wi_ces_contour_om_n.png", p5, width = 10, height = 7, dpi = 180, bg = "white")

# ── 9. Parameter Table Export ──────────────────────────────────────────────────
param_table_labeled <- param_table |>
  mutate(
    Parameter_Display = c(
      "$\\alpha_0$",
      "$\\alpha_1$",
      "$\\beta_0$",
      "$\\beta_1$",
      "$d$ (logit $\\delta$)",
      "$\\rho$"
    ),
    Parameter_Label = c(
      "Baseline Yield Shifter (bu/acre)",
      "Yield Shifter Sensitivity to OM",
      "Capital Composite Intercept",
      "Capital Composite Sensitivity to OM",
      "Logit Share Parameter (logit $\\delta$)",
      "Substitution Parameter"
    )
  ) |>
  select(Parameter, Parameter_Display, Parameter_Label,
         Estimate, Std_Error, t_value, p_value, CI_95_Low, CI_95_High)

write.csv(param_table_labeled, "./output/tables/table_wi_ces_parameters_labeled.csv", row.names = FALSE)
write.csv(param_table,         "./output/tables/table_wi_ces_parameters.csv",         row.names = FALSE)

latex_lines <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{CES Model Parameter Estimates for Wisconsin Corn Yield}",
  "\\label{tab:wi_ces_params}",
  "\\small",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{llrrrrrr}",
  "\\toprule",
  "Parameter & Interpretation & Estimate & Std. Error & t-value & p-value & 95\\% CI Low & 95\\% CI High \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(param_table_labeled))) {
  r <- param_table_labeled[i, ]
  latex_lines <- c(latex_lines,
    sprintf("%s & %s & %.5f & %.5f & %.3f & %.4f & %.5f & %.5f \\\\",
            r$Parameter_Display, r$Parameter_Label,
            r$Estimate, r$Std_Error, r$t_value, r$p_value, r$CI_95_Low, r$CI_95_High))
}

latex_lines <- c(
  latex_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "}",
  "\\begin{tablenotes}[flushleft]",
  "\\footnotesize",
  sprintf("\\item Notes: Derived $\\delta = %.4f$ (input share), $\\sigma = %.4f$ (elasticity of substitution). $\\delta$ recovered via $\\delta = 1/(1+e^{-d})$.", delta_est, sigma_est),
  "\\end{tablenotes}",
  "\\end{table}"
)

writeLines(latex_lines, "./output/tables/table_wi_ces_parameters.tex")
cat("\nOutputs saved to output/figures and output/tables.\n")

