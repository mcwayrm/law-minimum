# Production Capacity Estimation
# 2/18/26
# Nicholas Gallagher and Ryan McWay

# install.packages("rnassqs")
library("rnassqs")
library(readxl)

## NASS County data
# NASS provides API keys at https://quickstats.nass.usda.gov/api
NASSQS_TOKEN = 

# Set api key before requesting data
nassqs_auth(key = NASSQS_TOKEN)

############ NASS Queries ##########
till_vars <- c("PRACTICES, LAND USE, CROPLAND, CONSERVATION TILLAGE, (EXCL NO-TILL) - ACRES",
               "PRACTICES, LAND USE, CROPLAND, CONVENTIONAL TILLAGE - ACRES",
               "PRACTICES, LAND USE, CROPLAND, CONSERVATION TILLAGE, NO-TILL - ACRES")

till_params <- list(short_desc = till_vars,
                    agg_level_desc = c("COUNTY","STATE"),
                    domain_desc = "TOTAL")

till_df <- nassqs(till_params)

# add fips column to NASS data
till_df$fips <- rep(0,nrow(till_df))

# populate fips column one row at a time. This is an clean but slow approach. Tidyverse has other options
for(i in 1:nrow(till_df)){
  till_df[i,"fips"] <-  as.numeric(paste(till_df[i,"state_fips_code"],till_df[i,"county_ansi"],sep = ""))
}

# N Data
N_use <- data.frame(read_excel("N-P_from_fertilizer_1950-2017-july23-2020.xlsx", sheet = "farm", range = "B1:K3067"))[]

# Target years
years <- c(1987,1992,1997,2002,2007,2012,2017,2022)

# years of nitrogen data
Nyears <- c(1987,1992,1997,2002,2007,2012,2017)

# rename N use dataframe
names(N_use) <- c("fips","County","State",Nyears)


## Corn
# Target variables for corn
corn_vars <- c(
  "CORN - ACRES PLANTED",
  "CORN, GRAIN - YIELD, MEASURED IN BU / ACRE"
)

corn_params <- list(commodity_desc = "CORN", 
                    year = years,
                    short_desc = corn_vars,
                    agg_level_desc = "COUNTY",
                    reference_period_desc = "YEAR",
                    domain_desc = "TOTAL")

# Query Data
corn_df <- nassqs(corn_params)

# remove others
corn_df <- corn_df[corn_df$county_code != 998,]

# add fips column to NASS data
corn_df$fips <- rep(0,nrow(corn_df))

# populate fips column one row at a time. This is an clean but slow approach. Tidyverse has other options
for(i in 1:nrow(corn_df)){
  corn_df[i,"fips"] <-  as.numeric(paste(corn_df[i,"state_fips_code"],corn_df[i,"county_ansi"],sep = ""))
}

## CROPLAND
crop_vars <- "AG LAND, CROPLAND - ACRES"
  
crop_params <- list(year = years,
                    short_desc = crop_vars,
                    agg_level_desc = "COUNTY",
                    reference_period_desc = "YEAR",
                    domain_desc = "TOTAL")

cropland_df <- nassqs(crop_params)

# add fips column to NASS data
cropland_df$fips <- rep(0,nrow(cropland_df))

# populate fips column one row at a time. This is an clean but slow approach. Tidyverse has other options
for(i in 1:nrow(cropland_df)){
  cropland_df[i,"fips"] <-  as.numeric(paste(cropland_df[i,"state_fips_code"],cropland_df[i,"county_ansi"],sep = ""))
}
#### Regression Dataframe

reg_data <- expand.grid(years,unique(corn_df$fips))
names(reg_data) <- c("year","fips")

reg_data$yield <- rep(0, nrow(reg_data))
reg_data$corn_acres <- rep(0, nrow(reg_data))
reg_data$crop_acres <- rep(0, nrow(reg_data))
reg_data$n_use <- rep(0, nrow(reg_data))
reg_data$cnv_till <- rep(0, nrow(reg_data))
reg_data$red_till <- rep(0, nrow(reg_data))
reg_data$no_till <- rep(0, nrow(reg_data))

for (f in unique(reg_data$fips)){
  for (y in years){
    yld <- corn_df[(corn_df$fips == f & corn_df$year == y & corn_df$short_desc == "CORN, GRAIN - YIELD, MEASURED IN BU / ACRE"),"Value"]
    if (length(yld)==0){
      reg_data[reg_data$fips == f & reg_data$year == y,"yield"] <- NA
    } else {
      reg_data[reg_data$fips == f & reg_data$year == y,"yield"] <- yld
    }
    corn_ac <- corn_df[(corn_df$fips == f & corn_df$year == y & corn_df$short_desc == "CORN - ACRES PLANTED"), "Value"]
    if (length(corn_ac)==0){
      reg_data[reg_data$fips == f & reg_data$year == y,"corn_acres"] <- NA
    } else {
      reg_data[reg_data$fips == f & reg_data$year == y,"corn_acres"] <- corn_ac
    }
    crop_ac <- cropland_df[(cropland_df$fips == f & cropland_df$year == y & cropland_df$short_desc == "AG LAND, CROPLAND - ACRES"), "Value"]
    if (length(crop_ac)==0){
      reg_data[reg_data$fips == f & reg_data$year == y,"crop_acres"] <- NA
    } else {
      reg_data[reg_data$fips == f & reg_data$year == y,"crop_acres"] <- crop_ac
    }
    n <- N_use[N_use$fips == f,as.character(y)]
    if (length(n)==0){
      reg_data[reg_data$fips == f & reg_data$year == y,"n_use"] <- NA
    } else {
      reg_data[reg_data$fips == f & reg_data$year == y,"n_use"] <- n*2.20462
    }
    cnv_till <- till_df[(till_df$fips == f & till_df$year == y & till_df$short_desc == "PRACTICES, LAND USE, CROPLAND, CONVENTIONAL TILLAGE - ACRES"),"Value"]
    if (length(cnv_till)==0){
      reg_data[reg_data$fips == f & reg_data$year == y,"cnv_till"] <- NA
    } else {
      reg_data[reg_data$fips == f & reg_data$year == y,"cnv_till"] <- cnv_till
    }
    red_till <- till_df[(till_df$fips == f & till_df$year == y & till_df$short_desc == "PRACTICES, LAND USE, CROPLAND, CONSERVATION TILLAGE, (EXCL NO-TILL) - ACRES"),"Value"]
    if (length(red_till)==0){
      reg_data[reg_data$fips == f & reg_data$year == y,"red_till"] <- NA
    } else {
      reg_data[reg_data$fips == f & reg_data$year == y,"red_till"] <- red_till
    }
    no_till <- till_df[(till_df$fips == f & till_df$year == y & till_df$short_desc == "PRACTICES, LAND USE, CROPLAND, CONSERVATION TILLAGE, NO-TILL - ACRES"),"Value"]
    if (length(no_till)==0){
      reg_data[reg_data$fips == f & reg_data$year == y,"no_till"] <- NA
    } else {
      reg_data[reg_data$fips == f & reg_data$year == y,"no_till"] <- no_till
    }
  }
}

reg_data$n_use_ac_corn <- reg_data$n_use/reg_data$corn_acres
reg_data$n_use_ac_all <- reg_data$n_use/reg_data$crop_acres

reg_data$perc_no_till <- reg_data$no_till/(reg_data$no_till+reg_data$cnv_till+reg_data$red_till)
reg_data$perc_red_till <- (reg_data$no_till + reg_data$red_till)/(reg_data$no_till+reg_data$cnv_till+reg_data$red_till)

reg_data$year <- as.factor(reg_data$year)
reg_data$fips <- as.factor(reg_data$fips)

summary(reg_data$corn_acres)
summary(reg_data$crop_acres)

reg_data_large_ac <- reg_data[reg_data$corn_acres/reg_data$crop_acres >= 0.40 & reg_data$corn_acres >= 100000,]

reg_data_large_ac2 <- reg_data_large_ac[rowSums(is.na(reg_data_large_ac)) != ncol(reg_data_large_ac) & !(is.na(reg_data_large_ac$n_use_ac_corn)),]

reg_data_large_ac3 <- reg_data[reg_data$corn_acres/reg_data$crop_acres >= 0.20,]


summary(reg_data_large_ac2$n_use_ac_corn)

######### Regressions

n_response_reg1 <- lm(yield ~ n_use_ac_corn + I(n_use_ac_corn^2), data = reg_data_large_ac2)

summary(n_response_reg1)

x_grid <- seq(min(reg_data_large_ac2$n_use_ac_corn), max(reg_data_large_ac2$n_use_ac_corn), length.out = 200)
y_hat  <- predict(n_response_reg1, newdata = data.frame(n_use_ac_corn = x_grid))

# Plot data
plot(reg_data_large_ac2$n_use_ac_corn, reg_data_large_ac2$yield, pch = 19)

# Add fitted quadratic
lines(x_grid, y_hat, lwd = 2)

###

n_response_reg2 <- lm(yield ~ n_use_ac_corn*perc_no_till, data = reg_data_large_ac3)

summary(n_response_reg2)

summary(reg_data_large_ac2$perc_no_till)


#### Comparison reg

reg_test <- lm(yield ~ perc_no_till + year + fips - 1, data = reg_data_large_ac3)

summary(reg_test)

######### Soil data for WI
WI_soil <- data.frame(read_excel("WI-DATCP-soil-summary-combined-1995-2019.xlsx", sheet = "Sheet1"))

# get WI data
reg_data$fips_num <- as.numeric(levels(reg_data$fips))[reg_data$fips]
WI_reg_data <- reg_data[reg_data$fips_num>=55000 & reg_data$fips_num<56000,]

unique(WI_soil$measure)

WI_reg_data$OM <- rep(NA,nrow(WI_reg_data))
for (f in unique(WI_reg_data$fips)){
  for (y in years){
    om <- WI_soil[WI_soil$county_fips == f & WI_soil$year == y & WI_soil$measure == "OM","mean"]
    if(length(om) != 0){
      WI_reg_data[WI_reg_data$fips == f & WI_reg_data$year == y,"OM"] <- om
    }
  }
}

WI_reg_data <- WI_reg_data[WI_reg_data$OM <= 3.5,]

WI_reg_om <- lm(yield ~ n_use_ac_all*OM, data = WI_reg_data)

summary(WI_reg_om)

WI_reg_om_fe <- lm(yield ~ n_use_ac_all*OM + year + fips - 1, data = WI_reg_data)

summary(WI_reg_om_fe)

WI_reg_om_sq <- lm(yield ~ n_use_ac_all*OM + I(n_use_ac_all^2), data = WI_reg_data)

summary(WI_reg_om_sq)

plot(WI_reg_data$yield,WI_reg_data$OM)

########## GMM Approach
library("gmm")
library("ggplot2")


# Standardize covariates to reduce collinearity and scale issues
WI_reg_data$OM_s <- as.numeric(scale(WI_reg_data$OM))
WI_reg_data$N_s <- as.numeric(scale(WI_reg_data$n_use_ac_all))

vars <- c("year","fips","yield","n_use_ac_all","OM","OM_s","N_s")
WI_reg_data_small <- na.omit(WI_reg_data[,vars])

# Structural model: y = a(1 - b*exp(-c * OM))*(1 - d*exp(-e * N)))
# Reparameterize a,b,d as exp(.) to enforce positivity
gmm_moments <- function(theta, data) {
  a <- theta[1]
  b <- theta[2]
  c <- theta[3]
  d <- theta[4]
  e <- theta[5]
  
  y <- data$yield
  OM <- data$OM_s
  N <- data$N_s
  
  y_hat <- a * (1 - b*exp(-c * OM)) * (1 - d*exp(-e * N))
  resid <- y - y_hat
  
  Z <- cbind(1, OM, N, OM^2, N^2)
  resid * Z
}

a0 <- max(WI_reg_data_small$yield, na.rm = TRUE)
start_vals <- c(a = a0,
                b = 0.1,
                c = 0.1,
                d = 0.1,
                e = 0.1)

gmm_fit <- gmm(
  g = gmm_moments,
  x = WI_reg_data_small,
  t0 = start_vals,
  type = "iterative",
  method = "BFGS"
)

summary(gmm_fit)

coef_theta <- coef(gmm_fit)
a_hat <- coef_theta[1]
b_hat <- coef_theta[2]
c_hat <- coef_theta[3]
d_hat <- coef_theta[4]
e_hat <- coef_theta[5]

WI_reg_data_small$yield_fitted <- a_hat * (1 - b_hat*exp(-c_hat*WI_reg_data_small$OM_s)) *
  (1-d_hat*exp(-e_hat*WI_reg_data_small$N_s))

ggplot(WI_reg_data_small, aes(x = OM, y = n_use_ac_all, color = yield_fitted)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_gradient(low = "#2C7BB6", high = "#D7191C") +
  labs(title = "Fitted Corn Yield by OM and B",
       x = "Organic Matter (OM)",
       y = "Boron (B)",
       color = "Fitted Yield") +
  theme_minimal()



om_mean <- mean(WI_reg_data_small$OM, na.rm = TRUE)
om_sd <- sd(WI_reg_data_small$OM, na.rm = TRUE)
n_mean <- mean(WI_reg_data_small$n_use_ac_all, na.rm = TRUE)
n_sd <- sd(WI_reg_data_small$n_use_ac_all, na.rm = TRUE)

om_seq <- seq(min(WI_reg_data_small$OM, na.rm = TRUE), max(WI_reg_data_small$OM, na.rm = TRUE), length.out = 60)
n_seq <- seq(min(WI_reg_data_small$n_use_ac_all, na.rm = TRUE), max(WI_reg_data_small$n_use_ac_all, na.rm = TRUE), length.out = 60)
grid <- expand.grid(OM = om_seq, N = n_seq)
grid$OM_s <- (grid$OM - om_mean) / om_sd
grid$N_s <- (grid$N - n_mean) / n_sd
grid$yield_fitted <- a_hat * (1 - b_hat*exp(-c_hat * grid$OM_s)) * (1 - d_hat*exp(-e_hat * grid$N_s))

ggplot(grid, aes(x = OM, y = N, fill = yield_fitted)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = yield_fitted), color = "white", alpha = 0.5) +
  scale_fill_gradient(low = "#2C7BB6", high = "#D7191C") +
  geom_point(data = WI_reg_data_small, aes(x = OM, y = n_use_ac_all), inherit.aes = FALSE,
             color = "black", alpha = 0.2, size = 0.9) +
  labs(title = "Fitted Yield Surface with Observations",
       x = "Organic Matter (OM)",
       y = "Nitrogen (N)",
       fill = "Fitted Yield") +
  theme_minimal()
