# Housekeeping ----------------------------------------------
# Authors: Fabienne Krauer, Felix Guenther
##########################################

library(tidyverse)
library(scales)
library(gridExtra)
library(kableExtra)

theme_set(theme_light())
source("func_helpers.R")

# fig format
figform <- "png"
dpi <- 600

TK_labels <- c("month 0", "month 1", "month 2", "month 3", "month 4", "month 5",
               "month 6", "month 7", "month 8", "month 9", "month 10", "month 11",
               "1", "2", "3", "4", "5-9", "10-14", "15-24", "25-34",
               "35-44", "45-54", "55-64", "65-74", "75+")


# Load and prep data ------------------------------------------

# population data
pop <- readRDS("data/pop_DE_25_2015_2019.rds")
pop_2019 <- pop %>%   
  filter(year==2019) %>%
  select(groupno, n_pop = n)

# AGI data
AGI <- readRDS("data/AGI_agestrat.rds")
AGI_total <- AGI %>% dplyr::group_by(year, calweek, quarter, date, time, season) %>% 
                      dplyr::summarise(npos=sum(npos, na.rm=T))
AGI_prop <- read_csv("data/data_AGI_prop.csv")
  
# TK data
TK <- readRDS("data/TK_25_quarter.rds")
TK$quarter <- as.numeric(as.factor(TK$qy))
TK_total <-  TK %>% 
            dplyr::group_by(season, qy, quarter) %>% 
            dplyr::summarise(nhosp=sum(cases, na.rm=T)) %>% 
            dplyr::ungroup() %>% arrange(qy) 

TK_prop <- read.csv("data/data_TK25_prop.csv")


# Seroconversion data
seroconv <- read.csv("data/seroconversion.csv")


# INEK data
inek_smry = read_csv2("data/INEK_derived_endpoint_rates.csv")


## fitted simulations

# Directly fitted
quantiles_sim_ts <- read_csv("output/fit_simulations/quantiles_sim_ts.csv")
quantiles_sim_ts_prop <- read_csv("output/fit_simulations/quantiles_sim_ts_prop.csv")
quantiles_sim_hosp <- read_csv("output/fit_simulations/quantiles_sim_hosp.csv")
quantiles_sim_hosp_prop <- read_csv("output/fit_simulations/quantiles_sim_hosp_prop.csv")
quantiles_sim_noconv <- read_csv("output/fit_simulations/quantiles_sim_noconv.csv")

# Other derived quantities
quantiles_sim_hl <- read_csv("output/fit_simulations/quantiles_sim_hl.csv") %>% 
                  dplyr::mutate(riskgroup = "low risk") 
quantiles_sim_hh <- read_csv("output/fit_simulations/quantiles_sim_hh.csv") %>% 
  dplyr::mutate(riskgroup = "high risk")

quantiles_sim_h_byrisk <- full_join(quantiles_sim_hl, quantiles_sim_hh)


quantiles_sim_prop_highrisk <- read_csv("output/fit_simulations/quantiles_sim_prop_highrisk.csv") %>% arrange(agegrp)
quantiles_sim_hosp_age <- read_csv("output/fit_simulations/quantiles_sim_hosp_age.csv")
quantiles_sim_ts_age <- read_csv("output/fit_simulations/quantiles_sim_ts_age.csv")

# merge all
fit_agi <- merge(AGI_total, quantiles_sim_ts, by="time", all=T)
fit_agi_prop <- merge(AGI_prop, quantiles_sim_ts_prop, by.x = "groupno", by.y="agegrp", all=T)
fit_hosp <- merge(TK_total, quantiles_sim_hosp, by="quarter", all=T)
fit_hosp_prop <- merge(TK_prop, quantiles_sim_hosp_prop, by.x = "groupno", by.y="agegrp", all=T)
fit_noconv <- merge(seroconv, quantiles_sim_noconv, by.x = "groupno", by.y="agegrp", all=T)

fit_hosp_age <- full_join(TK, quantiles_sim_hosp_age, by = c("groupno" = "agegrp", "quarter"))
fit_agi_age <- full_join(AGI, quantiles_sim_ts_age, by = c("groupno" = "agegrp", "time"))


## INEK data

inek = read_delim("data/InEK-Daten_RSV-spezifisch_2019-2023.csv", 
                  delim = ";", 
                  locale = locale(decimal_mark = ","))


# Under-ascertainment -------------------------------

inek_hosp = inek %>% filter(year==2019, Parameter == "hosp") %>%
  mutate(age_start = agegroup,
         age_end = replace_na(lead(agegroup, 1), 100),
         age_ref = (age_end+age_start)/2,
         source = "INEK") %>%
  select(year, age_start, age_end, age_ref, n_hosp=value, n_pop = population,
         inc_100k = incidence, source)

data_mod_fit = quantiles_sim_hosp_age %>%
  mutate(year = case_when(quarter %in% c(1:2) ~ 2015,
                          quarter %in% c(3:6) ~ 2016,
                          quarter %in% c(7:10) ~ 2017,
                          quarter %in% c(11:14) ~ 2018,
                          quarter %in% c(15:16) ~ 2019)) %>%
  filter(year>=2016, year <2019) %>%
  select(-quarter) %>%
  rename(groupno = agegrp) %>%
  group_by(groupno, year) %>%
  summarise(n_hosp = sum(sim_median),
            n_hosp_q025 = sum(sim_low95),
            n_hosp_q975 = sum(sim_up95)) %>%
  mutate(age_start = case_when(groupno == 1 ~ 0,
                               groupno == 2 ~ 1/12, 
                               groupno == 3~ 2/12, 
                               groupno == 4~ 3/12, 
                               groupno == 5~ 4/12, 
                               groupno == 6~ 5/12, 
                               groupno == 7~ 6/12, 
                               groupno == 8~ 7/12, 
                               groupno == 9~ 8/12, 
                               groupno == 10~9/12, 
                               groupno == 11~10/12,
                               groupno == 12~11/12, 
                               groupno == 13~1, 
                               groupno == 14~2, 
                               groupno == 15~3, 
                               groupno == 16~4, 
                               groupno == 17~5, 
                               groupno == 18~10, 
                               groupno == 19~15, 
                               groupno == 20~25, 
                               groupno == 21~35, 
                               groupno == 22~45, 
                               groupno == 23~55, 
                               groupno == 24~65,
                               groupno == 25~75)) %>%
  group_by(year) %>%
  mutate(age_end = replace_na(lead(age_start, 1), 100)) %>%
  ungroup() %>%
  mutate(age_ref = (age_start + age_end)/2) %>%
  left_join(pop_2019) %>%
  mutate(inc_100k = n_hosp/n_pop*100000,
         inc_100k_q025 = n_hosp_q025/n_pop*100000,
         inc_100k_q975 = n_hosp_q975/n_pop*100000,
         source = "model fit") %>%
  select(year, age_start, age_end, age_ref, n_hosp, n_pop,
         inc_100k, inc_100k_q025, inc_100k_q975, source)

#' External data from population-based studies/systematic reviews
ext_data = tibble(age_start=c(18, 65, 75, 85, 60, 60), 
                  age_end = c(65, 75, 85, 100, 100, 100),
                  inc_100k = c(c(0.03, 0.64, 2.13, 2.10)*100,
                               0.15/100*100000,
                               1/1000*100000),
                  study = c(rep("Osei-Yeboah\n(2023)",4),
                            "Savic (2022)",
                            "Shi (2020)"))
#' Summary of calibrated model
dat_mod_smry = data_mod_fit %>%
  group_by(age_ref, age_start, age_end) %>%
  summarise(inc_100k = mean(inc_100k),
            inc_100k_q025 = mean(inc_100k_q025),
            inc_100k_q975 = mean(inc_100k_q975)) %>%
  mutate(source="Calibrated\nmodel")

#' Derived quantities of calibrated model
mod_deriv = data_mod_fit %>% filter(age_ref>=60) %>% 
  mutate(n_pop = ifelse(age_ref==60, n_pop/2, n_pop), 
         n_hosp = ifelse(age_ref==60, n_hosp/2, n_hosp)) %>%
  group_by(year) %>%
  summarise(inc_100k = sum(n_hosp)/sum(n_pop)*100000) %>%
  summarise(inc_100k = mean(inc_100k)) %>%
  mutate(age_start=60, age_end = 100, inc_100k, 
         study = "Calibrated\nmodel 60+")

# Under-ascertainment factors
mult_data = ext_data %>% 
  left_join(mod_deriv %>% 
              rename(inc_100k_mod = inc_100k) %>%
              select(-study),  
            by = c("age_start", "age_end")) %>%
  filter(!is.na(inc_100k_mod)) %>%
  rbind(ext_data %>% 
          left_join(dat_mod_smry %>% ungroup() %>% 
                      rename(inc_100k_mod = inc_100k) %>%
                      select(-c(inc_100k_q025, inc_100k_q975, age_end, age_ref)), 
                    by = c("age_start")) %>%
          filter(!is.na(inc_100k_mod)) %>%
          select(-source)) %>%
  rbind(ext_data %>% 
          left_join(dat_mod_smry %>% ungroup() %>% 
                      rename(inc_100k_mod = inc_100k) %>%
                      select(-c(inc_100k_q025, inc_100k_q975, age_end, age_ref)) %>%
                      filter(age_start==75) %>%
                      mutate(age_start=85), 
                    by = c("age_start")) %>%
          filter(!is.na(inc_100k_mod)) %>%
          select(-source)) %>%
  mutate(mult = inc_100k / inc_100k_mod)




# Plots and results ------------------------------------------


## Proportion high risk among all hospitalised infants -----------------

quantiles_sim_prop_highrisk %>% 
  filter(agegrp <=12) %>% 
  summary()



## Fig. 2 ----------------------

# TK time series
(fig2a <- ggplot(fit_hosp) +
  labs(title = "Total RSV-specific hospitalisations (TK)", tag = "A") +
  geom_errorbar(aes(x=quarter, ymin=sim_low95, ymax = sim_up95), colour="red") +
  geom_bar(aes(x=quarter, y=sim_median), fill="red", alpha=0.5, stat="identity") +
  geom_point(aes(x=quarter, y=nhosp)) +
  ylab("quarterly hospitalisations") + xlab(NULL) +
  scale_x_continuous(breaks = seq(1,16), labels=fit_hosp$qy) +
  theme(axis.text.x = element_text(angle=90),
        panel.grid.minor.x = element_blank()))  

# TK proportions
(fig2b <- ggplot(fit_hosp_prop) +
  labs(title = "Distribution of hospitalisations by age groups (TK)", tag = "B") +
  geom_bar(aes(x=reorder(label, groupno), y=sim_median), fill="red", alpha=0.5, stat="identity") +
  geom_errorbar(aes(x=reorder(label, groupno), ymin=sim_low95, ymax=sim_up95), colour="red") +
  geom_point(aes(x=reorder(label, groupno), y=prop)) +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle=90)) +
  ylab(NULL) + xlab("age group"))

# AGI time series
(fig2c <- ggplot(fit_agi) +
  labs(title = "Total weekly outpatient cases (AGI)", tag = "C") +
  geom_point(aes(x=date, y=npos)) +
  geom_line(aes(x=date, y=sim_median), colour="red") +
  geom_ribbon(aes(x=date, ymin=sim_low95, ymax=sim_up95), fill="red", alpha=0.4) +
  ylab("cases") + xlab(NULL))


# AGI proportions
(fig2d <- ggplot(fit_agi_prop) +
  labs(title = "Distribution of cases by age groups (AGI)", tag = "D") +
  geom_bar(aes(x=reorder(label, groupno), y=sim_median), fill="red", alpha=0.5, stat="identity") +
  geom_errorbar(aes(x=reorder(label, groupno), ymin=sim_low95, ymax=sim_up95), colour="red") +  
  geom_point(aes(x=reorder(label, groupno), y=prop)) +
  scale_y_continuous(labels = scales::percent) +
  ylab(NULL) + xlab("age group (years)"))


# Seroconversion data
(fig2e <- ggplot(fit_noconv) +
  labs(title = "Percentage never infected by age group", tag = "E") +
  geom_bar(aes(x=groupno, y=sim_median), fill="red", stat="identity", alpha=0.5) +
  geom_errorbar(aes(x=groupno, ymin=sim_low95, ymax=sim_up95), colour="red") +
  geom_point(aes(x=groupno, y=prop_noconv)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_x_continuous(breaks = 0:11) +
  ylab(NULL) + xlab("age group (months)"))

# combine
fig2 <- grid.arrange(fig2a, fig2b, fig2c, fig2d, fig2e, nrow = 4, layout_matrix = rbind(c(1,1),
                                                                                        c(2,2),
                                                                                        c(3,3),
                                                                                        c(4,5)))
ggsave(paste0("output/figures/fig2.", figform) , 
       fig2, width = 8, height = 12, 
       dpi = dpi)


## Fig. S8 ------------------------------------------

(figS8a <- ggplot(fit_hosp_age) +
   labs(title = "RSV-specific hospitalisations by age group (TK)", tag = "A") +
   facet_wrap(~reorder(label, agegrp), scale="free_y", ncol=4) +
   geom_errorbar(aes(x=quarter, ymin=sim_low95, ymax = sim_up95), colour="red") +
   geom_bar(aes(x=quarter, y=sim_median), stat="identity", fill="red", alpha=0.5) +
   geom_point(aes(x=quarter, y=cases)) +
   ylab("quarterly hospitalisations") + xlab(NULL) +
   scale_x_continuous(breaks = seq(1,16), labels=unique(fit_hosp_age$qy)) +
   theme(axis.text.x = element_text(angle=90),
         panel.grid.minor.x = element_blank()))


(figS8b <- ggplot(fit_agi_age[!is.na(fit_agi_age$agegrp),]) +
    labs(title = "RSV-specific outpatients by age group (AGI)", tag = "B") +
    facet_wrap(~reorder(agegrp, groupno), scale="free_y") +
    geom_ribbon(aes(x=date, ymin=sim_low95, ymax = sim_up95), fill="red", alpha=0.5) +
    geom_line(aes(x=date, y=sim_median), colour="red") +
    geom_point(aes(x=date, y=npos)) +
    ylab("weekly cases") + xlab(NULL) +
    theme(axis.text.x = element_text(angle=90),
          panel.grid.minor.x = element_blank()))

# combine
figS8 <- grid.arrange(figS8a, figS8b, nrow = 2, layout_matrix = rbind(1,1,1,1,2,2))

ggsave(paste0("output/figures/figS8.", figform), 
       figS8, width = 8, height = 12, dpi = dpi)





## Fig. S10 ---------------------------------------------

(figS10a <- ggplot(quantiles_sim_h_byrisk[quantiles_sim_h_byrisk$agegrp<=12,]) +
   labs(title = "Probability of hospitalisation, <1-year old",
        tag = "A") +
   geom_line(aes(x=agegrp, y=sim_median, colour=riskgroup)) +
   geom_ribbon(aes(x=agegrp, ymin=sim_low95, ymax = sim_up95, fill=riskgroup), alpha=0.5) +
   scale_y_continuous(labels = scales::percent) + 
   scale_x_continuous(breaks = 1:12) +
   theme(panel.grid.minor.x = element_blank(),
         legend.position = c(0.8, 0.8)) +
   ylab(NULL) + xlab("age group (months)")) 

(figS10b <- ggplot(quantiles_sim_h_byrisk[quantiles_sim_h_byrisk$agegrp>12,]) +
    labs(title = "Probability of hospitalisation, 1+ year old",
         tag = "B") +
    geom_line(aes(x=agegrp, y=sim_median), colour="#00BFC4") +
    geom_ribbon(aes(x=agegrp, ymin=sim_low95, ymax = sim_up95), fill="#00BFC4", alpha=0.5) +
    scale_y_continuous(labels = scales::percent) + 
    scale_x_continuous(breaks = 13:25, labels = TK_labels[13:25]) +
    theme(panel.grid.minor.x = element_blank()) +
    ylab(NULL) + xlab("age group (years)")) 


figS10 <- grid.arrange(figS10a, figS10b, nrow = 1)

ggsave(paste0("output/figures/figS10.", figform), 
       figS10, width = 12, height = 6, dpi = dpi)




## Fig. S11 ------------------------------------

cols_sel = palette.colors(7)[c(2:4, 6:7)]
names(cols_sel) = c("Savic (2022)", "Calibrated\nmodel", "Osei-Yeboah\n(2023)",
                    "Calibrated\nmodel 60+", "Shi (2020)")
cols_sel = cols_sel[c("Calibrated\nmodel",
                      "Calibrated\nmodel 60+",
                      "Osei-Yeboah\n(2023)",
                      "Calibrated\nmodel 60+",
                      "Savic (2022)",
                      "Shi (2020)")]

figS11 <- inek_hosp %>% 
  filter(age_ref>=55) %>% 
  ggplot() + 
  # INEK data
  geom_col(aes(age_ref, inc_100k, fill = source), 
           width=inek_hosp %>% filter(age_ref>=55) %>% 
             mutate(age_range = age_end-age_start) %>% 
             pull(age_range), alpha=.25, col = "black") +
  scale_fill_manual(values = c("INEK: year 2019"="grey")) +
  # Model output
  geom_segment(aes(x=age_start, xend = age_end, y= inc_100k_q025, yend = inc_100k_q025,
                   col = source),
               data = dat_mod_smry %>% filter(age_ref>55)) +
  geom_segment(aes(x=age_start, xend = age_start, y= inc_100k_q025, yend = inc_100k_q975,
                   col = source),
               data = dat_mod_smry %>% filter(age_ref>55)) +
  geom_segment(aes(x=age_end, xend = age_end, y= inc_100k_q025, yend = inc_100k_q975,
                   col = source),
               data = dat_mod_smry %>% filter(age_ref>55)) +
  geom_segment(aes(x=age_start, xend = age_end, y= inc_100k_q975, 
                   yend = inc_100k_q975,
                   col = source),
               data = dat_mod_smry %>% filter(age_ref>55)) +
  geom_segment(aes(x=age_start, xend = age_end, y= inc_100k, yend = inc_100k,
                   col = source),
               data = dat_mod_smry %>%
                 select(age_ref, 
                        age_start, 
                        age_end,
                        inc_100k,
                        source) %>%
                 filter(age_ref>55), lwd = 1.5
  ) +
  geom_segment(aes(x=age_start, xend = age_end, y= inc_100k, yend = inc_100k, 
                   col = study),
               data = mod_deriv, lwd = 1.5) +
  geom_segment(aes(x=age_start, xend = age_end, y= inc_100k, yend = inc_100k, 
                   col = study),
               data = ext_data %>% filter(age_start>=60), lwd = 1.5) +
  ylab("RSV hospitalization incidence\nrate per 100.000, yearly, log-10") +
  xlab("Age") +
  coord_cartesian(xlim=c(55, 100)) +
  expand_limits(y=225) +
  scale_y_continuous(trans="log10") +
  theme(legend.title = element_blank()) +
  geom_segment(aes(x=age_start, xend = age_start, y=inc_100k_mod, yend = inc_100k),
               data = mult_data,
               lty = 2) +
  geom_label(aes(age_start, 
                 exp((log(inc_100k) + log(inc_100k_mod))/2), 
                 label = paste0("x ", round(mult,1)),
                 hjust = ifelse(study=="Savic (2022)", 1, 0),
                 col = study),
             data = mult_data, show.legend=FALSE) + 
  scale_color_manual(values = cols_sel) +
  guides(fill  = guide_legend(order = 1),
         color = guide_legend(order = 2))

ggsave(figS11, 
       filename = paste0("output/figures/figS11.", figform),
       width = 10, height = 4,
       dpi = dpi)


