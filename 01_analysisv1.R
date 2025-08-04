
################################################################################
### Excess Heat Mortality Estimation: US Weekly CDC Data (2018–2024)
################################################################################

# Clear environment
rm(list = ls()); cat("\014");

# Load necessary libraries
suppressMessages(library(tidyverse))
suppressMessages(library(lubridate))
suppressMessages(library(dlnm))
suppressMessages(library(splines))
suppressMessages(library(mixmeta))
suppressMessages(library(tsModel))
suppressMessages(library(MASS))

################################################################################
### Load and Prepare Data
################################################################################

# Replace with your data path or import method
final_250k <- readRDS("weekly_250k.RDS")

# Rename and convert relevant variables
DATATABLE <- final_250k %>%
  rename(
    mort = X_deaths,
    temp = mean_temp_C,
    date = X_week_end_date,
    popu = population_estimate,
    location = StCoFIPS
  ) %>%
  mutate(date = as.Date(date))

# Region vector (county FIPS codes)
vREG <- unique(DATATABLE$location)
nREG <- length(vREG)

################################################################################
### Define Periods and Constants
################################################################################

baseline_years <- c(2018, 2019, 2021, 2022)
DATE1_CALI = as.Date("2018-01-04")
DATE2_CALI = as.Date("2022-12-29")
DATE1_PRED = as.Date("2023-01-05")
DATE2_PRED = as.Date("2024-12-26")
DATE1_SU23 = as.Date("2023-05-04")
DATE2_SU23 = as.Date("2023-09-28")
DATE1_SU24 = as.Date("2024-05-02")
DATE2_SU24 = as.Date("2024-09-26")

# Spline & DLNM setup
VAR_FUN = "ns"
VAR_PRC = c(10,50,90)/100
MIN_LAG = 0
MAX_LAG = 3
DF_SEAS = 8
PRED_PRC = seq(0,1,by=0.01)
# Temperature Percentile Range for the Minimum Mortality Temperature
MIN_PMMT =   5 / 100; if( any( 0 > MIN_PMMT | MIN_PMMT > 1 ) ){ stop("ERROR: Invalid Lower Temperature Percentile Range for the Minimum Mortality Temperature !!!"); }
MAX_PMMT = 100 / 100; if( any( 0 > MAX_PMMT | MAX_PMMT > 1 ) ){ stop("ERROR: Invalid Upper Temperature Percentile Range for the Minimum Mortality Temperature !!!"); }


################################################################################
### Format Data Lists
################################################################################

DATATABLE_CALI <- DATATABLE %>%
  filter(date >= DATE1_CALI & date <= DATE2_CALI + weeks(MAX_LAG))

DATATABLE_PRED <- DATATABLE %>%
  filter(date >= DATE1_PRED & date <= DATE2_PRED + weeks(MAX_LAG))

# Filter calibration data to baseline years only
DATALIST_CALI <- lapply(vREG, function(x) {
  subset(DATATABLE_CALI, location == x & year %in% baseline_years)
}); names(DATALIST_CALI) <- vREG

DATALIST_PRED <- lapply(vREG, function(x) {
  subset(DATATABLE_PRED, location == x)
}); names(DATALIST_PRED) <- vREG

# Filter out empty counties from both lists
non_empty <- sapply(DATALIST_CALI, nrow) > 0
DATALIST_CALI <- DATALIST_CALI[non_empty]
DATALIST_PRED <- DATALIST_PRED[names(DATALIST_CALI)]
vREG <- names(DATALIST_CALI)
nREG <- length(vREG)

# Add week-of-period
for(i in 1:nREG) {
  DATALIST_CALI[[i]]$wop <- seq_len(nrow(DATALIST_CALI[[i]]))
}

################################################################################
### Model Fitting and Prediction Setup
################################################################################

# Initialize storage
COEF_MODEL <- matrix(NA, nREG, length(VAR_PRC)+1, dimnames = list(vREG))
VCOV_MODEL <- vector("list", nREG); names(VCOV_MODEL) <- vREG
MMT_REG = array( NA, dim = c( nREG ), dimnames = list( vREG ) );
CROSS_PRED_REG_META <- vector("list", nREG); names(CROSS_PRED_REG_META) <- vREG

#check:
summary(DATALIST_CALI[[i]]$temp)
MMT_REG[i]


CROSS_PRED_REG_META[[i]] <- crosspred(
  basis,
  coef = BLUP[[i]]$blup,
  vcov = BLUP[[i]]$vcov,
  model.link = "log",
  at = temp_vals,
  cen = MMT_REG[i]  # <-- MUST be included
)


#fit models
for(i in 1:nREG){
  data <- DATALIST_CALI[[i]]
  
  if(sum(!is.na(data$temp)) < 5) {
    message(paste("Skipping", vREG[i], "- not enough temp data"))
    next
  }
  
  cb <- crossbasis(data$temp, lag=c(MIN_LAG,MAX_LAG),
                   argvar=list(fun=VAR_FUN,
                               knots=quantile(data$temp, VAR_PRC, na.rm=TRUE),
                               Boundary.knots=range(data$temp, na.rm=TRUE)),
                   arglag=list(fun="integer"))
  
  model <- glm(mort ~ ns(wop, df=round(DF_SEAS*nrow(data)*7/365.25)) + cb,
               data=data, family=quasipoisson)
  
  pred <- crosspred(cb, model, at=quantile(data$temp, PRED_PRC, na.rm=TRUE))
  
  mmt <- pred$predvar[
    which.min(pred$allRRfit[which(PRED_PRC == MIN_PMMT):which(PRED_PRC == MAX_PMMT)]) +
      which(PRED_PRC == MIN_PMMT) - 1
  ]
  
  red <- crossreduce(cb, model, cen=mmt)
  
  COEF_MODEL[i,] <- coef(red)
  VCOV_MODEL[[i]] <- vcov(red)
  MMT_REG[i] <- mmt
}


# Meta-analysis
TEMP_AVG <- sapply(DATALIST_CALI, function(x) mean(x$temp, na.rm=TRUE))
TEMP_IQR <- sapply(DATALIST_CALI, function(x) IQR(x$temp, na.rm=TRUE))

# Identify valid rows: no NA coefficients
valid_rows <- rowSums(is.na(COEF_MODEL)) == 0

# Subset everything to valid counties only
COEF_MODEL <- COEF_MODEL[valid_rows, , drop = FALSE]
VCOV_MODEL <- VCOV_MODEL[valid_rows]
TEMP_AVG <- TEMP_AVG[valid_rows]
TEMP_IQR <- TEMP_IQR[valid_rows]
vREG <- rownames(COEF_MODEL)
nREG <- length(vREG)

MULTIVAR <- mixmeta(COEF_MODEL ~ TEMP_AVG + TEMP_IQR, VCOV_MODEL,
                    data = data.frame(vREG=vREG), method="reml")
BLUP <- blup(MULTIVAR, vcov=TRUE)

# Generate cumulative exposure-response predictions
for(i in 1:nREG) {
  temp_vec <- DATALIST_CALI[[i]]$temp
  
  # Check: must have at least 5 non-missing temp values
  if (sum(!is.na(temp_vec)) < 5) {
    message(paste("Skipping", vREG[i], "- not enough valid temperature data for crosspred"))
    next
  }
  
  temp_vals <- quantile(temp_vec, PRED_PRC, na.rm = TRUE)
  temp_knots <- quantile(temp_vec, VAR_PRC, na.rm = TRUE)
  temp_bounds <- range(temp_vec, na.rm = TRUE)
  
  # Skip if knots or boundaries are non-finite (e.g., all NA or constant temp)
  if (any(!is.finite(temp_knots)) || any(!is.finite(temp_bounds))) {
    message(paste("Skipping", vREG[i], "- temperature quantiles invalid"))
    next
  }
  
  basis <- onebasis(temp_vals, fun = VAR_FUN,
                    knots = temp_knots,
                    Boundary.knots = temp_bounds)
  
  CROSS_PRED_REG_META[[i]] <- crosspred(basis,
                                        coef = BLUP[[i]]$blup,
                                        vcov = BLUP[[i]]$vcov,
                                        model.link = "log",
                                        at = temp_vals,
                                        cen = MMT_REG[i])
}


################################################################################
### Plot Cumulative Exposure-Response
################################################################################

pdf("cumulative_exposure_response.pdf", width=12, height=8)
layout(matrix(seq(1,length(vREG)), ncol=2))
for(i in 1:nREG){
  plot(CROSS_PRED_REG_META[[i]]$predvar, CROSS_PRED_REG_META[[i]]$allRRfit,
       type="l", lwd=2, col="black", main=paste("County:", vREG[i]),
       xlab="Temperature (°C)", ylab="Relative Risk")
  polygon(c(CROSS_PRED_REG_META[[i]]$predvar, rev(CROSS_PRED_REG_META[[i]]$predvar)),
          c(CROSS_PRED_REG_META[[i]]$allRRlow, rev(CROSS_PRED_REG_META[[i]]$allRRhigh)),
          col=rgb(0,0,0,0.2), border=NA)
  abline(h=1, lty=2)
  abline(v=MMT_REG[i], col="red", lty=3)
}
dev.off()


################################################################################
### Estimate Excess Deaths for Summer 2023 and 2024
################################################################################

# Define summer weeks for 2023 and 2024
summer_weeks_2023 <- which(DATATABLE_PRED$date >= DATE1_SU23 & DATATABLE_PRED$date <= DATE2_SU23)
summer_weeks_2024 <- which(DATATABLE_PRED$date >= DATE1_SU24 & DATATABLE_PRED$date <= DATE2_SU24)

# Function to calculate excess deaths
estimate_excess <- function(pred_data, blup_list, mmt_val, temp_var, mort_var) {
  basis <- onebasis(pred_data[[temp_var]], fun=VAR_FUN,
                    knots=quantile(pred_data[[temp_var]], VAR_PRC, na.rm=TRUE),
                    Boundary.knots=range(pred_data[[temp_var]], na.rm=TRUE))
  center <- onebasis(mmt_val, fun=VAR_FUN,
                     knots=quantile(pred_data[[temp_var]], VAR_PRC, na.rm=TRUE),
                     Boundary.knots=range(pred_data[[temp_var]], na.rm=TRUE))
  cb_cen <- scale(basis, center=center, scale=FALSE)
  rr <- exp(cb_cen %*% blup_list$blup)
  excess <- (rr - 1) * pred_data[[mort_var]]
  return(excess)
}

# Store excess death results
excess_summer_2023 <- list()
excess_summer_2024 <- list()

for(i in 1:nREG){
  pred_data <- DATALIST_PRED[[i]]
  blup_i <- BLUP[[i]]
  mmt_i <- MMT_REG[i]
  
  summer_data_2023 <- pred_data[pred_data$date >= DATE1_SU23 & pred_data$date <= DATE2_SU23, ]
  summer_data_2024 <- pred_data[pred_data$date >= DATE1_SU24 & pred_data$date <= DATE2_SU24, ]
  
  excess_2023 <- estimate_excess(summer_data_2023, blup_i, mmt_i, "temp", "mort")
  excess_2024 <- estimate_excess(summer_data_2024, blup_i, mmt_i, "temp", "mort")
  
  excess_summer_2023[[vREG[i]]] <- data.frame(date=summer_data_2023$date,
                                              excess_deaths=as.numeric(excess_2023))
  excess_summer_2024[[vREG[i]]] <- data.frame(date=summer_data_2024$date,
                                              excess_deaths=as.numeric(excess_2024))
}

# Combine results and write CSV
df_2023 <- bind_rows(excess_summer_2023, .id="FIPS")
df_2024 <- bind_rows(excess_summer_2024, .id="FIPS")

write_csv(df_2023, "excess_deaths_summer_2023.csv")
write_csv(df_2024, "excess_deaths_summer_2024.csv")

################################################################################
### Export RR Curves for All Counties
################################################################################

rr_export <- bind_rows(lapply(seq_along(CROSS_PRED_REG_META), function(i) {
  data.frame(
    FIPS = vREG[i],
    temp = CROSS_PRED_REG_META[[i]]$predvar,
    rr = CROSS_PRED_REG_META[[i]]$allRRfit,
    rr_low = CROSS_PRED_REG_META[[i]]$allRRlow,
    rr_high = CROSS_PRED_REG_META[[i]]$allRRhigh,
    mmt = MMT_REG[i]
  )
}))

write_csv(rr_export, "RR_curves_by_county.csv")
