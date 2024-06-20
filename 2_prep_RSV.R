# Housekeeping  -------------------------------------------------
# Author: Fabienne Krauer
#########################

library(tidyverse)
library(readr)
library(readxl)
library(lubridate)
library(binom)
library(ISOweek)
library(zoo)
library(mgcv)

#Sys.setlocale("LC_ALL","German")
theme_set(theme_minimal())


yearrange <- 2015:2019

# Prep ----------------------------------------

pop_25 <- readRDS("data/pop_DE_25_2015_2019.rds")

# Helper DF with all weeks from 2015-2019
time <- data.frame("date" = seq(as.Date("2015-01-01"), to=as.Date("2019-12-31"), by=1))
time$calweek <- isoweek(time$date)
#time$calweek <- as.numeric(strftime(time$date, format="%V"))  
time$year <- as.numeric(format(time$date, format="%Y"))
time$qy <- as.yearqtr(time$date, format = "%Y-%m-%d")
time$quarter <- quarters(time$date)
# Define season (=Q3,Q4 of year x and Q1,Q2 of year x+1)
time$season <- ifelse(time$quarter %in% c("Q3", "Q4"), time$year, 
                      ifelse(time$quarter=="Q2", time$year-1,NA))
time <- time[order(time$date),]
time$season[1] <- time$year[1]-1
for (i in 2:nrow(time)) {
  time$season[i] <- ifelse(is.na(time$season[i]), time$season[i-1], time$season[i])
}
rm(i)
# Restrict data to 4 complete seasons: 15/16-16/17-17/18-18/19
time <- time[time$season %in% c(2015:2018),]
time <- time[order(time$date),]
time$time <- julian(time$date, origin=as.Date(time$date[1]))


time_quarter_age <- time %>% 
                dplyr::select(season, qy, quarter) %>% unique() %>% 
                tidyr::crossing(pop_25[pop_25$year==2019, c("groupno", "age_low_y", "label", "agegrp")])

  
# TK data  ----------------------------------------

## 25 age groups -------------------------------------

TK <- data.frame(KHS_Fall=as.numeric(),
                 KW = as.character(),
                 Age_Group= as.character(),
                 Anteil=as.numeric(),
                 Anzahl=as.numeric(),
                 year=as.numeric())
for (i in yearrange) {
  
  moo <- read_delim(paste0("data/RKI/inzidenz_anteile_neu_s_", i, ".csv"), 
                    ";", escape_double = FALSE, 
                    locale = locale(decimal_mark = ","), 
                    trim_ws = TRUE)
  moo$Anzahl <- as.numeric(moo$Anzahl)
  moo$year <- i
  TK <- merge(TK, moo, by=intersect(names(TK), names(moo)),all=T)
  
}
rm(moo)

TK <- merge(time, 
            TK, 
            by.x=c("year", "calweek"), 
            by.y=c("year", "KW"), 
            all.x=T)

TK <- TK %>% filter(!is.na(Anzahl))


# Re-create the missing age group month 0-1 (currently lumped into mo0-2)
# assume month 0-1 is 68% of month 1-2, see Parukh 2017

foo <- TK %>% filter(Age_Group=="mo01") %>% 
              mutate(Age_Group="mo00",
               cases=Anzahl*0.68/1.68) 

TK$cases <- ifelse(TK$Age_Group=="mo01", TK$Anzahl/1.68, TK$Anzahl)
TK <- merge(TK, foo, by=intersect(names(TK), names(foo)), all=T)

TK$groupno <- as.numeric(factor(TK$Age_Group, levels= sort(unique(TK$Age_Group))))
TK$setting <- ifelse(TK$KHS_Fall==1, "inpatient", "outpatient")

TK <- merge(TK, 
               pop_25[pop_25$year==2019,c("groupno", "agegrp", "label", "age_low_y")], 
               by="groupno", all=T)


# calculate average proportion of hospitalised by age groups 1-12
# just for curiosity. Not needed in the model.
# prop_y1 <- TK %>% filter(grepl("^mo", Age_Group) & KHS_Fall==1) %>% 
#   dplyr::group_by(season, Age_Group) %>% 
#   dplyr::summarise(n = sum(cases, na.rm=T)) %>% 
#   dplyr::ungroup() %>% dplyr::group_by(season) %>% 
#   dplyr::mutate(N = sum(n, na.rm=T)) %>% 
#   dplyr::ungroup() %>% dplyr::group_by(Age_Group) %>% 
#   dplyr::summarise(n = sum(n, na.rm=T),
#                    N = sum(N, na.rm=T),
#                    prop = n/N)
# 
# ggplot(prop_y1) + 
#   geom_bar(aes(x=Age_Group, y=prop), stat="identity")
# 
# rm(prop_y1)

# Sum in-patients by season
# TK_25_in <- TK %>% filter(setting=="outpatient") %>% 
#   dplyr::group_by(season, label, agegrp, age_low_y, groupno) %>% 
#   dplyr::summarise(cases = sum(cases, na.rm=T)) %>% 
#   mutate(season = ifelse(season %in% c(2016, 2018),"high", "low")) %>% 
#   group_by(season, groupno, age_low_y) %>% 
#   dplyr::summarise(cases = mean(cases, na.rm=T)) %>% 
#   pivot_wider(id_cols = c(age_low_y, groupno), 
#               names_from = season,
#               values_from = cases) %>% 
#   dplyr::mutate(mean = (high + low)/2,
#                 amplitude = (high * 100 / mean - 100)/100)


# Sum in-patients by quarter, extend all quarters with missing data to all age groups
TK_25_quarter <- TK %>% 
                filter(setting=="inpatient") %>% 
                dplyr::group_by(setting, qy, season, quarter, label, agegrp, age_low_y, groupno) %>% 
                dplyr::summarise(cases = round(sum(cases, na.rm=T))) %>% 
                dplyr::ungroup() %>% 
                full_join(time_quarter_age) %>% 
                arrange(qy, groupno) %>% dplyr::mutate(setting = "inpatient") 

saveRDS(TK_25_quarter, "data/TK_25_quarter.rds")


# re-shape in-patient for fitting
TK_25_ts <- TK_25_quarter %>% 
                    arrange(qy) %>% 
                    pivot_wider(id_cols=qy, 
                                names_from = groupno, 
                                values_from = cases,
                                names_expand = TRUE)


write.csv(TK_25_ts, "data/data_TK25_hosp_quarter.csv", row.names = F, na="")


# calculate overall proportions of each age group: 25 age groups
TK_25_prop <- TK_25_quarter %>% 
              group_by(groupno) %>% 
              dplyr::summarise(num = sum(cases, na.rm=T)) %>% 
              dplyr::ungroup() %>% 
              dplyr::arrange(groupno) %>% 
              dplyr::mutate(denom = sum(num), 
                            prop = num/denom,
                            label_de = c("Monat 0", "Monat 1", "Monat 2", "Monat 3", "Monat 4", "Monat 5",
                                      "Monat 6", "Monat 7", "Monat 8", "Monat 9", "Monat 10", "Monat 11",
                                      "Jahr 1", "Jahr 2", "Jahr 3", "Jahr 4", "5-9", "10-14", "15-24", "25-34",
                                      "35-44", "45-54", "55-64", "65-74", "75+"),
                            label = c("month 0", "month 1", "month 2", "month 3", "month 4", "month 5",
                                         "month 6", "month 7", "month 8", "month 9", "month 10", "month 11",
                                         "year 1", "year 2", "year 3", "year 4", "5-9", "10-14", "15-24", "25-34",
                                         "35-44", "45-54", "55-64", "65-74", "75+"))




ggplot(TK_25_prop) + 
  geom_bar(aes(x=reorder(label, groupno), y=prop), stat="identity") +
  ylab("proportion") + xlab("age group") +
  theme(axis.text.x = element_text(angle=90))

write.csv(TK_25_prop, "data/data_TK25_prop.csv", na="", row.names = FALSE)


# AGI data --------------------------------------------

AGI <- read_excel("data/RKI/AGI_NRZ_Daten_fuer_FG33.xls")
AGI <- AGI %>% filter(AgeGroup!="unknown") %>% 
  dplyr::group_by(YearWWeek, AgeGroup) %>% 
  dplyr::summarise(ntested = n(),
                   npos = length(ResRSV[ResRSV=="P"]))
AGI$age_low_y <- as.numeric(substr(AGI$AgeGroup,1,2))
AGI$age_up_y <- as.numeric(gsub("[0-9]{2}\\.\\.", "", AGI$AgeGroup))+1
# Bin age groups 0-2, 2-5, 5-15, 15-35 and 35+ to match age structure in the model
AGI$agegrp <- cut(AGI$age_low_y, 
                  breaks = c(0,2,5,15,35,99), 
                  right = FALSE)

AGI$year <- as.numeric(substr(AGI$YearWWeek,1,4))
AGI$date <- ISOweek2date(paste(AGI$year, substr(AGI$YearWWeek,5,7), "7", sep="-"))
AGI$calweek <- isoweek(AGI$date)
AGI$quarter <- quarters(AGI$date)

AGI <- AGI %>% filter(!((quarter %in% c("Q1", "Q2") & year==2015) | (quarter %in% c("Q3", "Q4") & year==2019) | year>2019))

AGI <- AGI %>% group_by(year, calweek, quarter, date, agegrp) %>% 
  dplyr::summarise(ntested=sum(ntested, na.rm=T), npos=sum(npos, na.rm=T))

AGI <- merge(AGI, time[,c("date", "time", "season")], by="date", all.x=T)

AGI$groupno <- as.numeric(as.factor(AGI$agegrp))

saveRDS(AGI, "data/AGI_agestrat.rds")

# Calculate overall proportion of each age group
AGI_prop <- AGI %>% 
            group_by(agegrp, groupno) %>% 
            dplyr::summarise(cases = sum(npos, na.rm = T)) %>% 
            dplyr::ungroup() %>% 
            dplyr::mutate(num = cases,
                             denom = sum(cases),
                             prop = num/denom) %>% 
            dplyr::rename(label = agegrp) %>% 
            arrange(label) 


ggplot(AGI_prop) + 
  geom_bar(aes(x=label, y=prop), stat="identity") +
  ylab("proportion") + xlab("age group") +
  theme(axis.text.x = element_text(angle=90))

write.csv(AGI_prop, "data/data_AGI_prop.csv", na="", row.names = FALSE)


ggplot(AGI) + 
  geom_line(aes(x=date, y=npos, colour=as.factor(agegrp)))


ggplot(AGI) + 
  geom_line(aes(x=date, y=npos)) +
  facet_wrap(~agegrp)

# Reshape to wide for fitting
AGI_fit_age <- AGI %>% 
  arrange(time) %>% 
  select(time, npos, agegrp) %>% 
  pivot_wider(id_cols = "time", names_from = agegrp, values_from = npos, names_expand = TRUE)

write.csv(AGI_fit_age, "data/data_AGI_ts_age.csv", na="", row.names = FALSE)


# Seroprevalence data ---------------------------------


pienter <- read.csv("https://raw.githubusercontent.com/Stijn-A/RSV_serology/master/data/infection_status_csv.txt",
                    sep=",")

# Group age into intervals 
seroconv <- pienter %>% 
            mutate(agegrp= cut(pienter$age_days,
            breaks= seq(0, 365, by=30.25), 
            include.lowest = T, right=F)) %>% 
            dplyr::group_by(agegrp) %>% 
            dplyr::summarise(agemid=round(median(age_days)), # Age midpoint in age group
                             N = n(), # Total N in age group
                             noconv = N - sum(infection)) %>%  # n seroconverted in age group after infection
                            filter(!is.na(agegrp)) %>% 
            dplyr::ungroup() %>% 
            dplyr::arrange(agegrp) %>% 
            dplyr::mutate(prop_noconv = noconv/N,
                          label = c("month 0", "month 1", "month 2", "month 3", "month 4", "month 5",
                                    "month 6", "month 7", "month 8", "month 9", "month 10"),
                          groupno = 1:11)


write.csv(seroconv, "data/seroconversion.csv", row.names = FALSE)
