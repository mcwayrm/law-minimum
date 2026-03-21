# ============================================================
# Wisconsin County-Level Corn Production Nonlinear Model
# Stochastic / Resampling Version
#
# Model:
#   Yield = (a0 + a1*X1) * [1 - (b0 + b1*X1) * exp(-exp(c0+c1*X1) * X2)]
#
#   y  = corn_grain_yield  (county corn grain yield, bu/acre)
#   X1 = OM_mean           (mean soil organic matter, %)
#   X2 = N_rate_scaled     (farm N fertilizer, 100 kg per planted corn acre)
#
# Added stochastic components:
#   1. Random train/test split
#   2. Randomized multi-start estimation for nlsLM
#   3. Bootstrap resampling for parameter confidence intervals
#
# Reproducibility is controlled via set.seed(SEED).
# ============================================================

library(readxl)
library(readr)
library(rnassqs)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(minpack.lm)

SEED <- 123
TRAIN_PROP <- 0.80
N_RANDOM_STARTS <- 12
N_BOOT <- 20

set.seed(SEED)

# ── 1. Load Soil Data ──────────────────────────────────────────────────────────
setwd("E:/law-minimum/")
soil_raw <- read_excel("./data/WI-DATCP-soil-summary-combined-1995-2019.xlsx")

om_data <- soil_raw |>
  filter(measure == "OM") |>
  select(county_fips, county, year, OM_mean = `mean`) |>
  mutate(
    year = as.character(year),
    county_fips = as.character(county_fips)
  )

# ── 2. Load & Interpolate N Fertilizer Data ───────────────────────────────────
farm_raw <- read_excel("./data/N-P_from_fertilizer_1950-2017-july23-2020.xlsx",
                       sheet = "farm")

wi_farm <- farm_raw |>
  filter(State == "WI") |>
  select(STCOFIPS, CountyName, matches("farmfertN-kg-")) |>
  pivot_longer(
    cols = matches("farmfertN-kg-"),
    names_to = "year_col",
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
      x = as.numeric(year[!is.na(N_kg)]),
      y = N_kg[!is.na(N_kg)],
      xout = as.numeric(year)
    )$y
  ) |>
  filter(year >= 1995, year <= 2017) |>
  ungroup() |>
  select(county_fips, year, N_kg) |>
  mutate(year = as.character(year))

# ── 3. Load Corn Yield Data ───────────────────────────────────────────────────
NASSQS_TOKEN <- trimws(readr::read_file("./scripts/nassqs-api-key.txt"))
nassqs_auth(key = NASSQS_TOKEN)

params_query <- list(
  group_desc = "FIELD CROPS",
  short_desc = c("CORN, GRAIN - YIELD, MEASURED IN BU / ACRE"),
  agg_level_desc = c("COUNTY"),
  state_alpha = "WI",
  reference_period_desc = c("YEAR")
)

df_wi_corn_yield <- setDT(nassqs(params_query))
df_wi_corn_yield <- df_wi_corn_yield[, c(
  "state_alpha", "state_fips_code", "county_code", "county_name", "year", "Value"
)]
setnames(df_wi_corn_yield, old = "Value", new = "corn_grain_yield")
df_wi_corn_yield[, county_fips := paste0(state_fips_code, county_code)]
df_wi_corn_yield[, year := as.character(year)]
df_wi_corn_yield[, corn_grain_yield := as.numeric(gsub(",", "", corn_grain_yield))]

params <- list(group_desc = "FIELD CROPS", 
                short_desc = c("CORN - ACRES PLANTED"), 
                agg_level_desc = c("COUNTY"),
                state_alpha = "WI",
                reference_period_desc = c("YEAR"))
df_wi_corn_acres <- setDT(nassqs(params))
df_wi_corn_acres <- df_wi_corn_acres[, c("state_alpha", "state_fips_code", 
    "county_code", "county_name", "year", "Value")]
setnames(df_wi_corn_acres, old = "Value", new = "corn_acres")
df_wi_corn_acres[, county_fips := paste0(state_fips_code, county_code)]
df_wi_corn_acres[, year := as.character(year)]
df_wi_corn_acres[, corn_acres := as.numeric(gsub(",", "", corn_acres))]

# ── 4. Merge & Clean ──────────────────────────────────────────────────────────
panel <- om_data |>
  inner_join(n_annual, by = c("county_fips", "year")) |>
  inner_join(df_wi_corn_yield, by = c("county_fips", "year")) |>
  inner_join(df_wi_corn_acres, by = c("county_fips", "year")) |>
  filter(
    !is.na(OM_mean),
    !is.na(N_kg),
    !is.na(corn_grain_yield),
    !is.na(corn_acres),
    OM_mean > 0,
    N_kg > 0,
    corn_grain_yield > 0,
    corn_acres > 0
  ) |>
  mutate(
    N_rate_kg_acre = N_kg / corn_acres,
    N_rate_scaled = N_rate_kg_acre / 100
  )

cat(sprintf(
  "Panel: %d obs | %d counties | %s-%s\n",
  nrow(panel), n_distinct(panel$county), min(panel$year), max(panel$year)
))

# ── 5. Random Train/Test Split ────────────────────────────────────────────────
train_n <- floor(TRAIN_PROP * nrow(panel))
train_idx <- sample.int(nrow(panel), size = train_n, replace = FALSE)

panel <- panel |>
  mutate(sample_split = if_else(row_number() %in% train_idx, "Train", "Test"))

train_panel <- panel[train_idx, , drop = FALSE]
test_panel <- panel[-train_idx, , drop = FALSE]

cat(sprintf(
  "Train/Test split: %d train | %d test | seed = %d\n",
  nrow(train_panel), nrow(test_panel), SEED
))

# ── 6. Nonlinear Model Function ───────────────────────────────────────────────
model_fn <- function(X1, X2, a0, a1, b0, b1, c0, c1) {
  A <- a0 + a1 * X1
  B <- b0 + b1 * X1
  C <- exp(c0 + c1 * X1)
  A * (1 - B * exp(-C * X2))
}

predict_from_params <- function(data, coeffs) {
  model_fn(
    X1 = data$OM_mean,
    X2 = data$N_rate_scaled,
    a0 = unname(coeffs["a0"]),
    a1 = unname(coeffs["a1"]),
    b0 = unname(coeffs["b0"]),
    b1 = unname(coeffs["b1"]),
    c0 = unname(coeffs["c0"]),
    c1 = unname(coeffs["c1"])
  )
}

calc_metrics <- function(actual, predicted) {
  residuals <- actual - predicted
  rss <- sum(residuals^2)
  tss <- sum((actual - mean(actual))^2)

  data.frame(
    RMSE = sqrt(mean(residuals^2)),
    Pseudo_R2 = if (tss > 0) 1 - rss / tss else NA_real_
  )
}

fit_nls_safe <- function(data, start_vals) {
  tryCatch(
    nlsLM(
      corn_grain_yield ~ model_fn(OM_mean, N_rate_scaled, a0, a1, b0, b1, c0, c1),
      data = data,
      start = start_vals,
      control = nls.lm.control(maxiter = 1000, ftol = 1e-12, ptol = 1e-12)
    ),
    error = function(e) NULL
  )
}

randomize_start <- function(base_start) {
  list(
    a0 = max(20, rnorm(1, mean = base_start$a0, sd = 25)),
    a1 = rnorm(1, mean = base_start$a1, sd = 4),
    b0 = rnorm(1, mean = base_start$b0, sd = 0.5),
    b1 = rnorm(1, mean = base_start$b1, sd = 0.2),
    c0 = rnorm(1, mean = base_start$c0, sd = 0.8),
    c1 = rnorm(1, mean = base_start$c1, sd = 0.4)
  )
}

# ── 7. Randomized Multi-Start Estimation ─────────────────────────────────────
base_start <- list(
  a0 = 120,
  a1 = 8,
  b0 = -1.5,
  b1 = 0.4,
  c0 = -4.0,
  c1 = 0.8
)

candidate_starts <- c(
  list(base_start),
  replicate(N_RANDOM_STARTS, randomize_start(base_start), simplify = FALSE)
)

fit_attempts <- lapply(seq_along(candidate_starts), function(i) {
  current_fit <- fit_nls_safe(train_panel, candidate_starts[[i]])
  if (is.null(current_fit)) {
    return(NULL)
  }

  list(
    attempt = i,
    start_vals = candidate_starts[[i]],
    fit = current_fit,
    rss = sum(residuals(current_fit)^2)
  )
})

fit_attempts <- Filter(Negate(is.null), fit_attempts)

if (length(fit_attempts) == 0) {
  stop("No nlsLM fits converged across randomized starting values.")
}

best_attempt_id <- which.min(vapply(fit_attempts, function(x) x$rss, numeric(1)))
best_attempt <- fit_attempts[[best_attempt_id]]
fit <- best_attempt$fit

cat(sprintf(
  "Multi-start convergence: %d of %d attempts succeeded | best attempt = %d\n",
  length(fit_attempts), length(candidate_starts), best_attempt$attempt
))

# ── 8. Parameter Table and Performance ───────────────────────────────────────
s <- summary(fit)
fit_params <- coef(fit)
se_vals <- s$coefficients[, "Std. Error"]
t_vals <- s$coefficients[, "t value"]
p_vals <- s$coefficients[, "Pr(>|t|)"]
ci_lo <- fit_params - 1.96 * se_vals
ci_hi <- fit_params + 1.96 * se_vals

train_pred <- predict_from_params(train_panel, fit_params)
test_pred <- predict_from_params(test_panel, fit_params)
full_pred <- predict_from_params(panel, fit_params)

train_metrics <- calc_metrics(train_panel$corn_grain_yield, train_pred)
test_metrics <- calc_metrics(test_panel$corn_grain_yield, test_pred)

metrics_table <- bind_rows(
  cbind(Sample = "Train", train_metrics),
  cbind(Sample = "Test", test_metrics)
)

cat("\n=== Sample Performance ===\n")
print(metrics_table)

# ── 9. Bootstrap Resampling ──────────────────────────────────────────────────
boot_mat <- matrix(
  NA_real_,
  nrow = N_BOOT,
  ncol = length(fit_params),
  dimnames = list(NULL, names(fit_params))
)

boot_start <- as.list(fit_params)

for (boot_id in seq_len(N_BOOT)) {
  boot_idx <- sample.int(nrow(train_panel), size = nrow(train_panel), replace = TRUE)
  boot_data <- train_panel[boot_idx, , drop = FALSE]
  boot_fit <- fit_nls_safe(boot_data, boot_start)

  if (!is.null(boot_fit)) {
    boot_mat[boot_id, ] <- coef(boot_fit)
  }
}

boot_success <- complete.cases(boot_mat)
boot_success_n <- sum(boot_success)

cat(sprintf("Bootstrap convergence: %d of %d resamples succeeded\n", boot_success_n, N_BOOT))

if (boot_success_n > 1) {
  boot_ci <- t(apply(boot_mat[boot_success, , drop = FALSE], 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  colnames(boot_ci) <- c("Bootstrap_CI_95_Low", "Bootstrap_CI_95_High")
  boot_ci <- as.data.frame(boot_ci)
  boot_ci$Parameter <- rownames(boot_ci)
  rownames(boot_ci) <- NULL
} else {
  boot_ci <- data.frame(
    Parameter = names(fit_params),
    Bootstrap_CI_95_Low = NA_real_,
    Bootstrap_CI_95_High = NA_real_
  )
}

param_table <- data.frame(
  Parameter = c("a0", "a1", "b0", "b1", "c0 (log-scale)", "c1 (log-scale)"),
  Estimate = round(fit_params, 5),
  Std_Error = round(se_vals, 5),
  t_value = round(t_vals, 3),
  p_value = round(p_vals, 4),
  CI_95_Low = round(ci_lo, 5),
  CI_95_High = round(ci_hi, 5),
  row.names = NULL
)

param_table$Parameter_Key <- c("a0", "a1", "b0", "b1", "c0", "c1")
param_table <- param_table |>
  left_join(boot_ci, by = c("Parameter_Key" = "Parameter")) |>
  mutate(
    Bootstrap_CI_95_Low = round(Bootstrap_CI_95_Low, 5),
    Bootstrap_CI_95_High = round(Bootstrap_CI_95_High, 5)
  ) |>
  select(-Parameter_Key)

cat("\n=== Parameter Estimates ===\n")
print(param_table)

full_metrics <- calc_metrics(panel$corn_grain_yield, full_pred)

cat(sprintf(
  "\nFull-sample pseudo R2 = %.4f\nFull-sample RMSE = %.4f bu/acre\n",
  full_metrics$Pseudo_R2,
  full_metrics$RMSE
))

cat("\nEffective decay c = exp(c0 + c1*OM) at key OM values:\n")
for (om_val in c(2, 3, 4, 5)) {
  c_eff <- exp(fit_params["c0"] + fit_params["c1"] * om_val)
  cat(sprintf("  OM = %.1f  ->  c_eff = %.4f\n", om_val, c_eff))
}

panel <- panel |>
  mutate(
    sample_split = factor(sample_split, levels = c("Train", "Test")),
    yield_fitted = full_pred,
    resid = corn_grain_yield - full_pred
  )

# ── 10. Plots ─────────────────────────────────────────────────────────────────
BLUE <- "#2C7BB6"
RED <- "#D73027"
LBLUE <- "#91BFDB"
GRAY <- "#6E6E6E"

p1 <- ggplot(panel, aes(x = corn_grain_yield, y = yield_fitted, color = sample_split)) +
  geom_point(alpha = 0.35, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, color = RED, linetype = "dashed", linewidth = 1.1) +
  scale_color_manual(values = c(BLUE, GRAY)) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = sprintf(
      "Train RMSE = %.1f\nTest RMSE = %.1f\nBootstrap n = %d",
      train_metrics$RMSE,
      test_metrics$RMSE,
      boot_success_n
    ),
    hjust = -0.1,
    vjust = 1.3,
    size = 3.7,
    color = "#333333"
  ) +
  labs(
    title = "A — Fitted vs Observed Corn Yield",
    x = "Observed Corn Yield (bu/acre)",
    y = "Fitted Corn Yield (bu/acre)",
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.background = element_rect(fill = "#f9f9f9")
  )

N_quants <- quantile(panel$N_rate_scaled, c(0.10, 0.50, 0.90))
OM_seq <- seq(quantile(panel$OM_mean, 0.02), quantile(panel$OM_mean, 0.98), length.out = 300)

surface_b <- expand.grid(
  OM_mean = OM_seq,
  N_label = c("Low N (10th pct)", "Median N (50th pct)", "High N (90th pct)")
) |>
  mutate(
    N_rate_scaled = dplyr::case_when(
      N_label == "Low N (10th pct)" ~ N_quants[1],
      N_label == "Median N (50th pct)" ~ N_quants[2],
      N_label == "High N (90th pct)" ~ N_quants[3]
    ),
    N_rate_kg_acre = 100 * N_rate_scaled,
    yield_fitted = model_fn(OM_mean, N_rate_scaled, fit_params["a0"], fit_params["a1"], fit_params["b0"], fit_params["b1"], fit_params["c0"], fit_params["c1"]),
    N_label = factor(N_label, levels = c("Low N (10th pct)", "Median N (50th pct)", "High N (90th pct)"))
  )

p2 <- ggplot(surface_b, aes(x = OM_mean, y = yield_fitted, color = N_label)) +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = c(LBLUE, BLUE, RED)) +
  labs(
    title = "B — Fitted Corn Yield vs OM | Three N Levels",
    x = "Organic Matter (%)",
    y = "Fitted Corn Yield (bu/acre)",
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "#f9f9f9")
  )

OM_quants <- quantile(panel$OM_mean, c(0.10, 0.50, 0.90))
N_seq <- seq(quantile(panel$N_rate_kg_acre, 0.02), quantile(panel$N_rate_kg_acre, 0.98), length.out = 300)

surface_c <- expand.grid(
  N_rate_kg_acre = N_seq,
  OM_label = c("Low OM (10th pct)", "Median OM (50th pct)", "High OM (90th pct)")
) |>
  mutate(
    OM_mean = dplyr::case_when(
      OM_label == "Low OM (10th pct)" ~ OM_quants[1],
      OM_label == "Median OM (50th pct)" ~ OM_quants[2],
      OM_label == "High OM (90th pct)" ~ OM_quants[3]
    ),
    N_rate_scaled = N_rate_kg_acre / 100,
    yield_fitted = model_fn(OM_mean, N_rate_scaled, fit_params["a0"], fit_params["a1"], fit_params["b0"], fit_params["b1"], fit_params["c0"], fit_params["c1"]),
    OM_label = factor(OM_label, levels = c("Low OM (10th pct)", "Median OM (50th pct)", "High OM (90th pct)"))
  )

p3 <- ggplot(surface_c, aes(x = N_rate_kg_acre, y = yield_fitted, color = OM_label)) +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = c(LBLUE, BLUE, RED)) +
  labs(
    title = "C — Fitted Yield vs N Fertilizer | Three OM Levels",
    x = "N Fertilizer (kg per planted corn acre)",
    y = "Fitted Corn Yield (bu/acre)",
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "#f9f9f9")
  )

resid_yr <- panel |>
  mutate(year_num = as.integer(year)) |>
  group_by(year_num) |>
  summarise(mean_r = mean(resid), se_r = sd(resid) / sqrt(n()), .groups = "drop")

p4 <- ggplot(resid_yr, aes(x = year_num, y = mean_r)) +
  geom_ribbon(
    aes(x = year_num, ymin = mean_r - 1.96 * se_r, ymax = mean_r + 1.96 * se_r),
    fill = BLUE,
    alpha = 0.22,
    inherit.aes = FALSE
  ) +
  geom_line(color = BLUE, linewidth = 1.4) +
  geom_point(color = BLUE, size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = RED, linewidth = 1.1) +
  scale_x_continuous(breaks = seq(1995, 2020, by = 5)) +
  labs(
    title = "D — Mean Residual by Year",
    subtitle = "Shaded band = 95% CI around mean residual",
    x = "Year",
    y = "Mean Residual (bu/acre)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "#f9f9f9")
  )

combined <- grid.arrange(
  p1, p2, p3, p4,
  ncol = 2,
  top = grid::textGrob(
    "Wisconsin County Panel — Stochastic Corn Yield Model (1995–2017)\nRandom split + randomized multi-start + bootstrap",
    gp = grid::gpar(fontsize = 13, fontface = "bold")
  )
)

ggsave(
  "./output/figures/figure_wi_model_visuals_stochastic.png",
  combined,
  width = 14,
  height = 10,
  dpi = 150,
  bg = "white"
)

write.csv(param_table, "./output/tables/table_wi_model_parameters_stochastic.csv", row.names = FALSE)
write.csv(metrics_table, "./output/tables/table_wi_model_metrics_stochastic.csv", row.names = FALSE)
