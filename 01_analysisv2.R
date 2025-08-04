# Required Libraries
library(dlnm)
library(mixmeta)
library(splines)
library(dplyr)
library(purrr)

# ---- SETTINGS ----
MIN_LAG <- 0
MAX_LAG <- 3
VAR_FUN <- "ns"
VAR_PRC <- c(0.1, 0.5, 0.9)
PRED_PRC <- seq(0.01, 0.99, by=0.01)
DF_SEAS <- 4

# ---- INPUT DATA ----
# Assume: DATALIST_CALI (list of county data), each with `mort`, `temp`
vREG <- names(DATALIST_CALI)
nREG <- length(vREG)

# Initialize result storage
COEF_MODEL <- matrix(NA, nREG, length(VAR_PRC) + 1)
VCOV_MODEL <- vector("list", nREG)
MMT_REG <- numeric(nREG)

# ---- MODEL FITTING PER COUNTY ----
for(i in 1:nREG){
  data <- DATALIST_CALI[[i]]
  
  if(sum(!is.na(data$temp)) < 5){
    message(paste("Skipping", vREG[i], "- not enough data"))
    next
  }
  
  # Add weekly order position for seasonal spline
  data$wop <- 1:nrow(data)
  
  cb <- crossbasis(data$temp, lag=c(MIN_LAG, MAX_LAG),
                   argvar = list(fun=VAR_FUN,
                                 knots=quantile(data$temp, VAR_PRC, na.rm=TRUE),
                                 Boundary.knots=range(data$temp, na.rm=TRUE)),
                   arglag = list(fun="integer"))
  
  model <- glm(mort ~ ns(wop, df=round(DF_SEAS * nrow(data) * 7 / 365.25)) + cb,
               data=data, family=quasipoisson)
  
  pred <- crosspred(cb, model, at=quantile(data$temp, PRED_PRC, na.rm=TRUE))
  
  # Restrict MMT search to 5–30°C
  mmt_candidates <- pred$predvar >= 5 & pred$predvar <= 30
  mmt <- pred$predvar[mmt_candidates][which.min(pred$allRRfit[mmt_candidates])]
  
  red <- crossreduce(cb, model, cen=mmt)
  
  COEF_MODEL[i,] <- coef(red)
  VCOV_MODEL[[i]] <- vcov(red)
  MMT_REG[i] <- mmt
}

# ---- META-REGRESSION ----
TEMP_AVG <- sapply(DATALIST_CALI[vREG], function(df) mean(df$temp, na.rm=TRUE))
TEMP_IQR <- sapply(DATALIST_CALI[vREG], function(df) IQR(df$temp, na.rm=TRUE))

meta_data <- data.frame(TEMP_AVG, TEMP_IQR)

MULTIVAR <- mixmeta(COEF_MODEL ~ TEMP_AVG + TEMP_IQR, VCOV_MODEL,
                    data = meta_data, method="reml")

# ---- BLUP + POOLED COEF ----
BLUP <- blup(MULTIVAR, vcov=TRUE)
coef_mat <- sapply(BLUP, function(x) x$blup)
vcov_list <- lapply(BLUP, function(x) x$vcov)

coef_avg <- rowMeans(coef_mat)
vcov_avg <- Reduce("+", vcov_list) / length(vcov_list)

# ---- POOLED PREDICTION ----
all_temp <- unlist(lapply(DATALIST_CALI[vREG], function(x) x$temp))
temp_grid <- quantile(all_temp, PRED_PRC, na.rm=TRUE)
temp_knots <- quantile(all_temp, VAR_PRC, na.rm=TRUE)
temp_bounds <- range(all_temp, na.rm=TRUE)

basis_mat <- onebasis(temp_grid, fun=VAR_FUN,
                      knots=temp_knots,
                      Boundary.knots=temp_bounds)

rr_vals <- exp(basis_mat %*% coef_avg)

# Restrict MMT search to 5–30°C
mmt_range <- which(temp_grid >= 5 & temp_grid <= 30)
mmt_overall <- temp_grid[mmt_range[which.min(rr_vals[mmt_range])]]

crosspred_all <- crosspred(basis_mat,
                           coef=coef_avg,
                           vcov=vcov_avg,
                           model.link="log",
                           at=temp_grid,
                           cen=mmt_overall)

# ---- PLOT ----
plot_range <- which(temp_grid >= -5 & temp_grid <= 40)

plot(temp_grid[plot_range], crosspred_all$allRRfit[plot_range],
     type = "l", lwd = 3, col = "black",
     xlab = "Temperature (°C)", ylab = "Relative Risk",
     main = "Pooled Exposure–Response Curve", ylim=c(0.9, max(crosspred_all$allRRhigh[plot_range])))

polygon(c(temp_grid[plot_range], rev(temp_grid[plot_range])),
        c(crosspred_all$allRRlow[plot_range], rev(crosspred_all$allRRhigh[plot_range])),
        col=rgb(0,0,0,0.2), border=NA)

abline(h=1, lty=2, col="gray")
abline(v=crosspred_all$cen, lty=3, col="red")

# ---- OPTIONAL: EXPORT CURVE ----
pooled_rr_df <- data.frame(
  temp = crosspred_all$predvar,
  RR = crosspred_all$allRRfit,
  RR_low = crosspred_all$allRRlow,
  RR_high = crosspred_all$allRRhigh
)

# write.csv(pooled_rr_df, "pooled_RR_curve.csv", row.names = FALSE)

#leave in percentiles in pooled exposure response.



# ---- Libraries ----
library(dplyr)
library(purrr)
library(dlnm)
library(readr)
library(ggplot2)

# ---- Assumptions ----
# - DATALIST_CALI: List of data frames per county with variables mort, temp, pop, week, year
# - crosspred_all: from pooled DLNM (output of crosspred)
# - region_df: data frame with columns county_id, region

# ---- Function to estimate weekly excess deaths ----
estimate_excess_deaths <- function(data, crosspred_all, mmt) {
  rr_lookup <- approxfun(crosspred_all$predvar, crosspred_all$allRRfit, rule = 2)
  
  data <- data %>%
    mutate(
      rr = rr_lookup(temp),
      af = ifelse(temp > mmt, (rr - 1) / rr, 0),  # Only apply if temp > MMT
      excess = af * mort
    )
  return(data)
}

# ---- Apply to all counties ----
excess_results <- map2_dfr(DATALIST_CALI, names(DATALIST_CALI), function(df, cname) {
  if (sum(!is.na(df$temp)) < 5) return(NULL)
  out <- estimate_excess_deaths(df, crosspred_all, crosspred_all$cen)
  out$county <- cname
  return(out)
})

# ---- Join with region info ----
region_df <- read.csv("US_States_FIPS_NCA_Region.csv")
excess_results <- left_join(excess_results, region_df, by = c("stname.x" = "StName"))

# ---- Summarise by year ----
year_summary <- excess_results %>%
  group_by(year) %>%
  summarise(
    total_excess = sum(excess, na.rm = TRUE),
    total_mort = sum(mort, na.rm = TRUE),
    pop = mean(popu, na.rm = TRUE),
    per_million = total_excess / pop * 1e5
  )

# ---- Summarise by region and year ----
region_year_summary <- excess_results %>%
  group_by(NCA4_Reg, year) %>%
  summarise(
    total_excess = sum(excess, na.rm = TRUE),
    pop = mean(popu, na.rm = TRUE),
    per_million = total_excess /pop*1e5
  )

# ---- Save summaries ----
write_csv(year_summary, "national_excess_by_year.csv")
write_csv(region_year_summary, "region_excess_by_year.csv")

# ---- Trend over time ----
trend_model <- lm(total_excess ~ year, data = year_summary)
summary(trend_model)

# ---- Plot national trend ----
ggplot(year_summary, aes(x = year, y = per_million)) +
  geom_line(color = "firebrick", size = 1.2) +
  ylim(0,10000)+
  geom_point(size = 2) +
  labs(
    title = "Heat-Attributable Excess Mortality per 100k",
    y = "Excess Deaths per 100k",
    x = "Year"
  ) +
  theme_minimal(base_size = 14)
