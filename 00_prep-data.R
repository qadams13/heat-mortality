library(dplyr)
library(lubridate)

## climate data
#read in ERA 5 population weighted data

clim <- list.files(pattern="v1.1.Rds", all.files=FALSE,
                   full.names=FALSE)
clim1<-lapply(clim, readRDS)
clim2 <-bind_rows(clim1)
clim2$month <- month(clim2$Date)
clim2$year <- year(clim2$Date)

clim3 <- clim2 %>%
  select(StCoFIPS, Date, month, year, Tmax_C)%>%
  group_by(month, year, StCoFIPS)%>%
  summarize(Tmax_C = mean(Tmax_C))

write.csv(clim3, "climatedata.csv")

##mortality data
#read in monthly data from CDC wonder
