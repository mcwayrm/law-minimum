# Title: Matching WI soil data to Corn data
# Blame: Ryan McWay
# Purpose: Playing with matching corn yield to WI historical soil data

###############################################################################
# SECTION: Preamble
###############################################################################
# Clear the environment
rm(list = ls())
# Set up the directory
setwd("E:/law-minimum")
# Load Packages
library("rio")
library("rnassqs")
library("data.table")
library("tidyverse")
library("ggplot2")
library("tibble")
library("fixest")
library("modelsummary")
library("broom")
library("vtable")
library("gmm")
library("cowplot")

# NASS provides API keys at https://quickstats.nass.usda.gov/api
NASSQS_TOKEN <- readr::read_file("./scripts/nassqs-api-key.txt")
# Set api key before requesting data
nassqs_auth(key = NASSQS_TOKEN)

###############################################################################
# SECTION: Clean, Merge Data
###############################################################################

# WI SOC ground truth data
df_soil_wi <- setDT(import("./data/soc-ground-truth/WI-soil-summary/WI-DATCP-soil-summary-combined-1995-2019.xlsx"))
    # Rename variables
    setnames(df_soil_wi, c("county_fips"), 
                    c("fips_code"))
    df_soil_wi <- df_soil_wi[, fips_code := as.character(fips_code)] 
    df_soil_wi <- df_soil_wi[, year := as.character(year)] 
# Keep only mean
df_soil_wi <- df_soil_wi[, c("measure", "fips_code", "year", "mean")]
    # NOTE: Maybe add in samples_submitted in future


# Seperate Data for each soil measure 
df_soc_wi <- df_soil_wi[measure == "OM",] # Organic Matter %
    df_soc_wi <- df_soc_wi[, -c("measure")]
    setnames(df_soc_wi, old = "mean", new = "OM")
df_b_wi <- df_soil_wi[measure == "B",] # Boron (ppm)
    df_b_wi <- df_b_wi[, -c("measure")]
    setnames(df_b_wi, old = "mean", new = "B")
df_bph_wi <- df_soil_wi[measure == "BPH",] # Buffer pH
    df_bph_wi <- df_bph_wi[, -c("measure")]
    setnames(df_bph_wi, old = "mean", new = "BPH")
df_ca_wi <- df_soil_wi[measure == "Ca",] # Calcium (ppm)
    df_ca_wi <- df_ca_wi[, -c("measure")]
    setnames(df_ca_wi, old = "mean", new = "Ca")
df_k_wi <- df_soil_wi[measure == "K",] # Potassium (ppm)
    df_k_wi <- df_k_wi[, -c("measure")]
    setnames(df_k_wi, old = "mean", new = "K")
df_mn_wi <- df_soil_wi[measure == "M",] # Manganese (ppm)
    df_mn_wi <- df_mn_wi[, -c("measure")]
    setnames(df_mn_wi, old = "mean", new = "M")
df_mg_wi <- df_soil_wi[measure == "Mg",] # Magnesium (ppm)
    df_mg_wi <- df_mg_wi[, -c("measure")]
    setnames(df_mg_wi, old = "mean", new = "Mg") 
df_p_wi <- df_soil_wi[measure == "P",] # Phosphorus (ppm)
    df_p_wi <- df_p_wi[, -c("measure")]
    setnames(df_p_wi, old = "mean", new = "P")
df_ph_wi <- df_soil_wi[measure == "PH",] # Soil pH
    df_ph_wi <- df_ph_wi[, -c("measure")]
    setnames(df_ph_wi, old = "mean", new = "PH")
df_s_wi <- df_soil_wi[measure == "S",] # Sulfur (ppm)
    df_s_wi <- df_s_wi[, -c("measure")]
    setnames(df_s_wi, old = "mean", new = "S")
df_zn_wi <- df_soil_wi[measure == "Zn",] # Zinc (ppm)
    df_zn_wi <- df_zn_wi[, -c("measure")]
    setnames(df_zn_wi, old = "mean", new = "Zn")

# TODO: Combine the soil measures via merge on year and fips_code
# df_soil_wi <- dcast(df_soil_wi, fips_code + year ~ measure, value.var = "mean")

# Bring in annual corn yield data
params <- list(group_desc = "FIELD CROPS", 
                short_desc = c("CORN, GRAIN - YIELD, MEASURED IN BU / ACRE"), 
                agg_level_desc = c("COUNTY"),
                state_alpha = "WI",
                reference_period_desc = c("YEAR"))
df_wi_corn_yield <- setDT(nassqs(params))
df_wi_corn_yield <- df_wi_corn_yield[, c("state_alpha", "state_fips_code", 
    "county_code", "county_name", "year", "Value")]
setnames(df_wi_corn_yield, old = "Value", new = "corn_grain_yield")
df_wi_corn_yield[, fips_code := paste0(state_fips_code, county_code)]
df_wi_corn_yield[, year := as.character(year)]

params <- list(group_desc = "FIELD CROPS", 
                short_desc = c("CORN, GRAIN - YIELD, MEASURED IN BU / NET PLANTED ACRE"), 
                agg_level_desc = c("COUNTY"),
                state_alpha = "WI",
                reference_period_desc = c("YEAR"))
df_wi_corn_yield_net <- setDT(nassqs(params))
df_wi_corn_yield_net <- df_wi_corn_yield_net[, c("state_alpha", "state_fips_code", 
    "county_code", "county_name", "year", "Value")]
setnames(df_wi_corn_yield_net, old = "Value", new = "corn_grain_net_yield")
df_wi_corn_yield_net[, fips_code := paste0(state_fips_code, county_code)]
df_wi_corn_yield_net[, year := as.character(year)]

params <- list(group_desc = "FIELD CROPS", 
                short_desc = c("CORN, SILAGE - YIELD, MEASURED IN TONS / ACRE"), 
                agg_level_desc = c("COUNTY"),
                state_alpha = "WI",
                reference_period_desc = c("YEAR"))
df_wi_silage_yield <- setDT(nassqs(params))
df_wi_silage_yield <- df_wi_silage_yield[, c( "state_alpha", "state_fips_code", 
    "county_code", "county_name", "year", "Value")]
setnames(df_wi_silage_yield, old = "Value", new = "corn_silage_yield")
df_wi_silage_yield[, fips_code := paste0(state_fips_code, county_code)]
df_wi_silage_yield[, year := as.character(year)]

# df_wi_corn <- left_join(df_wi_corn_yield, df_wi_corn_yield_net, 
#     by = c("fips_code", "year"))
# df_wi_corn <- left_join(df_wi_corn_yield, df_wi_silage_yield, 
#     by = c("fips_code", "year"))

# Merge corn yields with soil data
df_wi <- merge(df_wi_corn_yield, df_soc_wi, by = c("fips_code", "year"))
length(unique(df_wi$fips_code))
length(unique(df_wi$year))
df_wi <- merge(df_wi, df_b_wi, by = c("fips_code", "year"))
df_wi <- merge(df_wi, df_bph_wi, by = c("fips_code", "year"))
df_wi <- merge(df_wi, df_ca_wi, by = c("fips_code", "year"))
df_wi <- merge(df_wi, df_k_wi, by = c("fips_code", "year"))
df_wi <- merge(df_wi, df_mn_wi, by = c("fips_code", "year"))
df_wi <- merge(df_wi, df_mg_wi, by = c("fips_code", "year"))
df_wi <- merge(df_wi, df_p_wi, by = c("fips_code", "year"))
df_wi <- merge(df_wi, df_ph_wi, by = c("fips_code", "year"))
df_wi <- merge(df_wi, df_s_wi, by = c("fips_code", "year"))
df_wi <- merge(df_wi, df_zn_wi, by = c("fips_code", "year"))
length(unique(df_wi$fips_code))
length(unique(df_wi$year))
    # NOTE: Losing some years. But no counties. 

###############################################################################
# SECTION: Simple Regressions
###############################################################################

# Understand the data
colnames(df_wi)
sumtable(df_wi)
table(df_wi$year)

model_basic1 <- feols(
    corn_grain_yield ~ OM * B| fips_code + year,
    data = df_wi
)
model_basic2 <- feols(
    corn_grain_yield ~ OM * BPH| fips_code + year,
    data = df_wi
)
model_basic3 <- feols(
    corn_grain_yield ~ OM * Ca| fips_code + year,
    data = df_wi
)
model_basic4 <- feols(
    corn_grain_yield ~ OM * K| fips_code + year,
    data = df_wi
)
model_basic5 <- feols(
    corn_grain_yield ~ OM * M| fips_code + year,
    data = df_wi
)
model_basic6 <- feols(
    corn_grain_yield ~ OM * Mg| fips_code + year,
    data = df_wi
)
model_basic7 <- feols(
    corn_grain_yield ~ OM * P| fips_code + year,
    data = df_wi
)
model_basic8 <- feols(
    corn_grain_yield ~ OM * PH| fips_code + year,
    data = df_wi
)
model_basic9 <- feols(
    corn_grain_yield ~ OM * S| fips_code + year,
    data = df_wi
)
model_basic10 <- feols(
    corn_grain_yield ~ OM * Zn| fips_code + year,
    data = df_wi
)

modelsummary(list(model_basic1, model_basic2, model_basic3, model_basic4, model_basic5, model_basic6, model_basic7, 
    model_basic8, model_basic9, model_basic10), stars = TRUE)
# NOTE: 
    # OM signifiant for Ca, K, P, PH, Zn
    # Interaction signifcant for BPH, Ca, K, Mg, P, PH, Zn
modelsummary(list( model_basic3, model_basic4,  model_basic7, 
    model_basic8, model_basic10), stars = TRUE)

###############################################################################
# SECTION: Estimating Law of the Minimum
###############################################################################


# Standardize covariates to reduce collinearity and scale issues
df_wi$OM_s <- as.numeric(scale(df_wi$OM))
df_wi$Ca_s <- as.numeric(scale(df_wi$Ca))
df_wi$K_s <- as.numeric(scale(df_wi$K))
df_wi$P_s <- as.numeric(scale(df_wi$P))
df_wi$PH_s <- as.numeric(scale(df_wi$PH))
df_wi$Zn_s <- as.numeric(scale(df_wi$Zn))

# Structural model: y = a(1 - b exp((c OM))) *(1- d exp(e N)))
# Reparameterize a,b,d as exp(.) to enforce positivity
gmm_moments <- function(theta, data) {
    a <- exp(theta[1])
    c1 <- exp(theta[2])
    b1 <- theta[3]
    c2 <- exp(theta[4])
    b2 <- theta[5]
    c3 <- exp(theta[6])
    b3 <- theta[7]
    c4 <- exp(theta[8])
    b4 <- theta[9]
    c5 <- exp(theta[10])
    b5 <- theta[11]
    c6 <- exp(theta[12])
    b6 <- theta[13]

    y <- data$corn_grain_yield
    OM <- data$OM_s
    Ca <- data$Ca_s
    K <- data$K_s
    P <- data$P_s
    PH <- data$PH_s
    Zn <- data$Zn_s

    # Soil components
    comp_OM <- (1 - c1 * exp(b1 * OM))
    comp_Ca <- (1 - c2 * exp(b2 * Ca))
    comp_K <- (1 - c3 * exp(b3 * K))
    comp_P <- (1 - c4 * exp(b4 * P))
    comp_PH <- (1 - c5 * exp(b5 * PH))
    comp_Zn <- (1 - c6 * exp(b6 * Zn))
    y_hat <- a * comp_OM * comp_Ca * comp_K * comp_P * comp_PH * comp_Zn 
    resid <- y - y_hat

    Z <- cbind(1, OM, Ca, K, P, PH, Zn, 
        OM^2, Ca^2, K^2, P^2, PH^2, Zn^2)
    resid * Z
}

# a0 <- mean(df_wi$corn_grain_yield, na.rm = TRUE)
a0 <- max(df_wi$corn_grain_yield, na.rm = TRUE)
start_vals <- c(a = log(a0),
                c1 = log(0.1),
                b1 = 0.0,
                c2 = log(0.1),
                b2 = 0.0,
                c3 = log(0.1),
                b3 = 0.0,
                c4 = log(0.1),
                b4 = 0.0,
                c5 = log(0.1),
                b5 = 0.0,
                c6 = log(0.1),
                b6 = 0.0)

gmm_fit <- gmm(
    g = gmm_moments,
    x = df_wi,
    t0 = start_vals,
    type = "iterative"
)

summary(gmm_fit)

# Visualize fitted model against OM and B
coef_theta <- coef(gmm_fit)
a_hat <- exp(coef_theta[1])
c1_hat <- exp(coef_theta[2])
b1_hat <- coef_theta[3]
c2_hat <- exp(coef_theta[4])
b2_hat <- coef_theta[5]
c3_hat <- exp(coef_theta[6])
b3_hat <- coef_theta[7]
c4_hat <- exp(coef_theta[8])
b4_hat <- coef_theta[9]
c5_hat <- exp(coef_theta[10])
b5_hat <- coef_theta[11]
c6_hat <- exp(coef_theta[12])
b6_hat <- coef_theta[13]

df_wi$corn_grain_fitted <- a_hat * (1 - c1_hat * exp( b1_hat * df_wi$OM_s)) *
    (1 - c2_hat * exp( b2_hat * df_wi$Ca_s)) *
    (1 - c3_hat * exp( b3_hat * df_wi$K_s)) * 
    (1 - c4_hat * exp( b4_hat * df_wi$P_s)) * 
    (1 - c5_hat * exp( b5_hat * df_wi$PH_s)) * 
    (1 - c6_hat * exp( b6_hat * df_wi$Zn_s))

summary(df_wi$corn_grain_fitted)
summary(df_wi$corn_grain_yield)
ggplot(df_wi, aes(x = OM, y = Ca, color = corn_grain_fitted)) +
    geom_point(alpha = 0.7, size = 1.8) +
    scale_color_gradient(low = "#2C7BB6", high = "#D7191C") +
    labs(title = "Fitted Corn Yield by OM and Ca",
        x = "Organic Matter (OM)",
        y = "Calcium (Ca)",
        color = "Fitted Yield") +
    theme_minimal()



# Prepare Contour map grid (standardized)
om_mean <- mean(df_wi$OM, na.rm = TRUE)
om_sd <- sd(df_wi$OM, na.rm = TRUE)
om_seq <- seq(min(df_wi$OM, na.rm = TRUE), max(df_wi$OM, na.rm = TRUE), length.out = 10)
ca_mean <- mean(df_wi$Ca, na.rm = TRUE)
ca_sd <- sd(df_wi$Ca, na.rm = TRUE)
ca_seq <- seq(min(df_wi$Ca, na.rm = TRUE), max(df_wi$Ca, na.rm = TRUE), length.out = 10)
k_mean <- mean(df_wi$K, na.rm = TRUE)
k_sd <- sd(df_wi$K, na.rm = TRUE)
k_seq <- seq(min(df_wi$K, na.rm = TRUE), max(df_wi$K, na.rm = TRUE), length.out = 10)
p_mean <- mean(df_wi$P, na.rm = TRUE)
p_sd <- sd(df_wi$P, na.rm = TRUE)
p_seq <- seq(min(df_wi$P, na.rm = TRUE), max(df_wi$P, na.rm = TRUE), length.out = 10)
ph_mean <- mean(df_wi$PH, na.rm = TRUE)
ph_sd <- sd(df_wi$PH, na.rm = TRUE)
ph_seq <- seq(min(df_wi$PH, na.rm = TRUE), max(df_wi$PH, na.rm = TRUE), length.out = 10)
zn_mean <- mean(df_wi$Zn, na.rm = TRUE)
zn_sd <- sd(df_wi$Zn, na.rm = TRUE)
zn_seq <- seq(min(df_wi$Zn, na.rm = TRUE), max(df_wi$Zn, na.rm = TRUE), length.out = 10)


grid <- expand.grid(OM = om_seq, Ca = ca_seq, K = k_seq, P = p_seq, PH = ph_seq, Zn = zn_seq)
grid$OM_s <- (grid$OM - om_mean) / om_sd
grid$Ca_s <- (grid$Ca - ca_mean) / ca_sd
grid$K_s <- (grid$K - k_mean) / k_sd
grid$P_s <- (grid$P - p_mean) / p_sd
grid$PH_s <- (grid$PH - ph_mean) / ph_sd
grid$Zn_s <- (grid$Zn - zn_mean) / zn_sd

grid$corn_grain_fitted <- a_hat * (1 - c1_hat * exp( b1_hat * grid$OM_s)) *
    (1 - c2_hat * exp( b2_hat * grid$Ca_s)) *
    (1 - c3_hat * exp( b3_hat * grid$K_s)) * 
    (1 - c4_hat * exp( b4_hat * grid$P_s)) * 
    (1 - c5_hat * exp( b5_hat * grid$PH_s)) * 
    (1 - c6_hat * exp( b6_hat * grid$Zn_s))


# Contour Map: Calcium
p1 <- ggplot(grid, aes(x = OM, y = Ca, fill = corn_grain_fitted)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = corn_grain_fitted), color = "white", alpha = 0.5) +
    scale_fill_gradient(low = "#2C7BB6", high = "#D7191C") +
    geom_point(data = df_wi, aes(x = OM, y = Ca), inherit.aes = FALSE,
                color = "black", alpha = 0.2, size = 0.9) +
    labs(title = "Fitted Corn Yield by OM and Ca",
        x = "Organic Matter (OM)",
        y = "Calcium (Ca)",
        fill = "Fitted Yield") +
    theme_minimal()

# Contour Map: Potasium
p2 <- ggplot(grid, aes(x = OM, y = K, fill = corn_grain_fitted)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = corn_grain_fitted), color = "white", alpha = 0.5) +
    scale_fill_gradient(low = "#2C7BB6", high = "#D7191C") +
    geom_point(data = df_wi, aes(x = OM, y = K), inherit.aes = FALSE,
                color = "black", alpha = 0.2, size = 0.9) +
    labs(title = "Fitted Corn Yield by OM and K",
        x = "Organic Matter (OM)",
        y = "Potasium (K)",
        fill = "Fitted Yield") +
    theme_minimal()

# Contour Map: P
p3 <- ggplot(grid, aes(x = OM, y = P, fill = corn_grain_fitted)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = corn_grain_fitted), color = "white", alpha = 0.5) +
    scale_fill_gradient(low = "#2C7BB6", high = "#D7191C") +
    geom_point(data = df_wi, aes(x = OM, y = P), inherit.aes = FALSE,
                color = "black", alpha = 0.2, size = 0.9) +
    labs(title = "Fitted Corn Yield by OM and P",
        x = "Organic Matter (OM)",
        y = "Phosphrous (P)",
        fill = "Fitted Yield") +
    theme_minimal()

# Contour Map: Acidity
p4 <- ggplot(grid, aes(x = OM, y = PH, fill = corn_grain_fitted)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = corn_grain_fitted), color = "white", alpha = 0.5) +
    scale_fill_gradient(low = "#2C7BB6", high = "#D7191C") +
    geom_point(data = df_wi, aes(x = OM, y = PH), inherit.aes = FALSE,
                color = "black", alpha = 0.2, size = 0.9) +
    labs(title = "Fitted Corn Yield by OM and PH",
        x = "Organic Matter (OM)",
        y = "Acidity (PH)",
        fill = "Fitted Yield") +
    theme_minimal()

# Contour Map: Zinc
p5 <- ggplot(grid, aes(x = OM, y = Zn, fill = corn_grain_fitted)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = corn_grain_fitted), color = "white", alpha = 0.5) +
    scale_fill_gradient(low = "#2C7BB6", high = "#D7191C") +
    geom_point(data = df_wi, aes(x = OM, y = Zn), inherit.aes = FALSE,
                color = "black", alpha = 0.2, size = 0.9) +
    labs(title = "Fitted Corn Yield by OM and Zn",
        x = "Organic Matter (OM)",
        y = "Zinc (Zn)",
        fill = "Fitted Yield") +
    theme_minimal()


plot_grid(p1, p2, p3, p4, p5, ncol = 2)

summary(df_wi$corn_grain_fitted)


#######
# Not full interaction 
###### 

# Contour Map: Calcium
grid <- expand.grid(OM = om_seq, Ca = ca_seq)
grid$OM_s <- (grid$OM - om_mean) / om_sd
grid$Ca_s <- (grid$Ca - ca_mean) / ca_sd
grid$corn_grain_fitted <- a_hat * (1 - exp(-c1_hat * (b1_hat + grid$OM_s))) *
    (1 - exp(-c2_hat * (b2_hat + grid$Ca_s))) 
p1 <- ggplot(grid, aes(x = OM, y = Ca, fill = corn_grain_fitted)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = corn_grain_fitted), color = "white", alpha = 0.5) +
    scale_fill_gradient(low = "#2C7BB6", high = "#D7191C") +
    geom_point(data = df_wi, aes(x = OM, y = Ca), inherit.aes = FALSE,
                color = "black", alpha = 0.2, size = 0.9) +
    labs(title = "Fitted Corn Yield by OM and Ca",
        x = "Organic Matter (OM)",
        y = "Calcium (Ca)",
        fill = "Fitted Yield") +
    theme_minimal()

# Contour Map: Potasium
grid <- expand.grid(OM = om_seq, K = k_seq)
grid$OM_s <- (grid$OM - om_mean) / om_sd
grid$K_s <- (grid$K - k_mean) / k_sd
grid$corn_grain_fitted <- a_hat * (1 - exp(-c1_hat * (b1_hat + grid$OM_s))) *
    (1 - exp(-c3_hat * (b3_hat + grid$K_s))) 
p2 <- ggplot(grid, aes(x = OM, y = K, fill = corn_grain_fitted)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = corn_grain_fitted), color = "white", alpha = 0.5) +
    scale_fill_gradient(low = "#2C7BB6", high = "#D7191C") +
    geom_point(data = df_wi, aes(x = OM, y = K), inherit.aes = FALSE,
                color = "black", alpha = 0.2, size = 0.9) +
    labs(title = "Fitted Corn Yield by OM and K",
        x = "Organic Matter (OM)",
        y = "Potasium (K)",
        fill = "Fitted Yield") +
    theme_minimal()

# Contour Map: P
grid <- expand.grid(OM = om_seq, P = p_seq)
grid$OM_s <- (grid$OM - om_mean) / om_sd
grid$P_s <- (grid$P - p_mean) / p_sd
grid$corn_grain_fitted <- a_hat * (1 - exp(-c1_hat * (b1_hat + grid$OM_s))) *
    (1 - exp(-c4_hat * (b4_hat + grid$P_s))) 
p3 <- ggplot(grid, aes(x = OM, y = P, fill = corn_grain_fitted)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = corn_grain_fitted), color = "white", alpha = 0.5) +
    scale_fill_gradient(low = "#2C7BB6", high = "#D7191C") +
    geom_point(data = df_wi, aes(x = OM, y = P), inherit.aes = FALSE,
                color = "black", alpha = 0.2, size = 0.9) +
    labs(title = "Fitted Corn Yield by OM and P",
        x = "Organic Matter (OM)",
        y = "Phosphrous (P)",
        fill = "Fitted Yield") +
    theme_minimal()

# Contour Map: Acidity
grid <- expand.grid(OM = om_seq, PH = ph_seq)
grid$OM_s <- (grid$OM - om_mean) / om_sd
grid$PH_s <- (grid$PH - ph_mean) / ph_sd
grid$corn_grain_fitted <- a_hat * (1 - exp(-c1_hat * (b1_hat + grid$OM_s))) *
    (1 - exp(-c5_hat * (b5_hat + grid$PH_s))) 
p4 <- ggplot(grid, aes(x = OM, y = PH, fill = corn_grain_fitted)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = corn_grain_fitted), color = "white", alpha = 0.5) +
    scale_fill_gradient(low = "#2C7BB6", high = "#D7191C") +
    geom_point(data = df_wi, aes(x = OM, y = PH), inherit.aes = FALSE,
                color = "black", alpha = 0.2, size = 0.9) +
    labs(title = "Fitted Corn Yield by OM and PH",
        x = "Organic Matter (OM)",
        y = "Acidity (PH)",
        fill = "Fitted Yield") +
    theme_minimal()

# Contour Map: Zinc
grid <- expand.grid(OM = om_seq, Zn = zn_seq)
grid$OM_s <- (grid$OM - om_mean) / om_sd
grid$Zn_s <- (grid$Zn - zn_mean) / zn_sd
grid$corn_grain_fitted <- a_hat * (1 - exp(-c1_hat * (b1_hat + grid$OM_s))) *
    (1 - exp(-c6_hat * (b6_hat + grid$Zn_s))) 
p5 <- ggplot(grid, aes(x = OM, y = Zn, fill = corn_grain_fitted)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = corn_grain_fitted), color = "white", alpha = 0.5) +
    scale_fill_gradient(low = "#2C7BB6", high = "#D7191C") +
    geom_point(data = df_wi, aes(x = OM, y = Zn), inherit.aes = FALSE,
                color = "black", alpha = 0.2, size = 0.9) +
    labs(title = "Fitted Corn Yield by OM and Zn",
        x = "Organic Matter (OM)",
        y = "Zinc (Zn)",
        fill = "Fitted Yield") +
    theme_minimal()

plot_grid(p1, p2, p3, p4, p5, ncol = 2)
