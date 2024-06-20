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
library(gridExtra)
library(binom)
library(mgcv)
library(gridExtra)

#Sys.setlocale("LC_ALL","German")
theme_set(theme_minimal())

# fig format
figform <- "png"
dpi <- 600

# Load data ----------------------------------------

pop_25 <- readRDS("data/pop_DE_25_2015_2019.rds")
AGI <- readRDS("data/AGI_agestrat.rds")
TK_25_quarter <- readRDS("data/TK_25_quarter.rds")

asymptomatic <- read_excel("data/asymptomatics_raw.xlsx")
attackrates <- read_excel("data/attackrates_raw.xlsx")


# Asymptomatic data  --------------------------------

asymptomatic[,c("pa","pa_low95","pa_up95")] <- binom.confint(asymptomatic$pa_num, asymptomatic$pa_denom, method="exact")[,c("mean","lower","upper")]
asymptomatic$agemid <- asymptomatic$age_low_y + (asymptomatic$age_up_y-asymptomatic$age_low_y)/2


ggplot(asymptomatic) + 
  geom_linerange(aes(xmin=age_low_y, xmax=age_up_y, y=pa*100, colour=as.factor(ref_pa))) +
  xlab("Age (years)") + ylab("Asymptomatics (%)") 

# Interpolate data for other age groups
fit_pa <- gam(pa ~ s(agemid, bs="cr", k=5), 
              link=Gamma(type='log'),
              data=asymptomatic)


gam.check(fit_pa)

asymptomatic$pa_predicted <- predict(fit_pa)



# Fit vs. data
ggplot(asymptomatic) +
    geom_linerange(aes(xmin=age_low_y, xmax=age_up_y, y=pa*100, colour=as.factor(ref_pa))) +
    geom_point(aes(x=agemid, y=pa*100, colour=as.factor(ref_pa))) +
    geom_errorbar(aes(x=agemid, ymin=pa_low95*100, ymax=pa_up95*100, colour=as.factor(ref_pa))) +
    ylab("asymptomatics (%)") + xlab("age (years)") +
    geom_line(aes(x=agemid, y=pa_predicted*100)) +
    labs(colour="data from") +
    theme(legend.position = "bottom")


# Attack rate data  --------------------------------


# Age groups year 1 - 65 +
attackrates <- attackrates %>% 
  #filter(!(ref %in% c("Munywoki 2015", "Glezen 1986", "Glenn 2016")))  %>% 
  #filter(!(ref=="Takashima 2021" & age_low_y==0)) %>% 
  mutate(agemid = age_low_y + (age_up_y - age_low_y)/2)

ggplot(attackrates) + 
  geom_point(aes(x=agemid, y=ar*100)) +
  xlab("Age (years)") + ylab("attack rate (%)") +
  scale_x_continuous(trans="log")


# Fit exponential models 
attackrates$ar_predicted <- NA

fit_ar <- lm(ar ~ log(agemid),
              data=attackrates)
attackrates$ar_predicted <- predict(fit_ar)


ggplot(attackrates) +
  geom_point(aes(x=agemid, y=ar*100, colour=as.factor(ref))) +
  ylab("Attack rate (%)") + xlab("age (years)") +
  geom_line(aes(x=agemid, y=ar_predicted*100)) +
  scale_x_continuous(trans='log') +
  labs(colour="data from") +
    theme(legend.position = "bottom")


# Epidata --------------------------------------------


# Predict the attack rate (AR) for all the age groups in the model (based on the model fits)
epidata_25 <- pop_25 %>% filter(year==2019) %>% 
  select(groupno, agegrp, age_low_y, age_up_y, label, agemid, d_y, n, b0, b1, b2, deaths_d) %>% 
  arrange(agemid)

# Predict the Proportion asymptomatic (PA) for all the age groups in the model (based on the model fits)
epidata_25$pa <- NA
epidata_25$pa <- predict(fit_pa,
                             data.frame(agemid=epidata_25$agemid,
                                        pa=epidata_25$pa))


epidata_25 <- epidata_25 %>% arrange(groupno)

epidata_25$ar <- NA
epidata_25$ar[epidata_25$groupno>12] <- predict(fit_ar, 
                                                data.frame(agemid=epidata_25$agemid[epidata_25$groupno>12], 
                                                           ar=epidata_25$ar[epidata_25$groupno>12]))

# This is a place holder, we don't know these values
epidata_25$ar[epidata_25$groupno<=12] <- 0.5



# Hospitalisation probabilities h ---------------------------------------------

# Calculate the expected annual cumulative number of symptomatic cases by age group

# Age groups 13+
season <- TK_25_quarter %>% #filter(season > 2016) %>% 
  mutate(season = ifelse(season %in% c(2016, 2018),"high", "low")) %>% 
  group_by(season, groupno, age_low_y) %>% 
  dplyr::summarise(cases = mean(cases, na.rm=T)) %>% 
  pivot_wider(id_cols = c(age_low_y, groupno), 
              names_from = season,
              values_from = cases) %>% 
  dplyr::mutate(mean = (high + low)/2,
                amplitude = (high * 100 / mean - 100)/100)

summary(season$amplitude)
ggplot(season) + geom_point(aes(x=groupno, y=amplitude))

#fit_amp <- lm(amplitude ~ groupno, data=season)

#ggplot(season) + 
#  geom_point(aes(x=groupno, y=amplitude)) +
#  geom_line(aes(x=groupno, y=predict(fit_amp)))

#summary(fit_amp)

# --> amplitude in seasonl case load seems to be around +/- 35%


# Calculate expected number of hospitalised cases (assuming a biennial transmission pattern with amplitude of 35%)
 
underreporting <- TK_25_quarter %>%  
                dplyr::group_by(season, groupno, age_low_y) %>% 
                dplyr::summarise(cases_obs = sum(cases, na.rm=T)) %>% 
                merge(., epidata_25[,c("groupno", "age_low_y", "ar", "pa", "n", "label", "agegrp", "agemid")],
                      by=c("groupno", "age_low_y"), 
                      all=T) %>% arrange(groupno)

underreporting$cases_exp <- ifelse(underreporting$season %in% c(2016, 2018),
                          (underreporting$ar + underreporting$ar*0.35) * (1-underreporting$pa) * underreporting$n, 
                          (underreporting$ar - underreporting$ar*0.35) * (1-underreporting$pa) * underreporting$n)

underreporting$h <- underreporting$cases_obs / underreporting$cases_exp
underreporting$h <- ifelse(underreporting$h>1,1,underreporting$h)

summary(underreporting$h)


ggplot(underreporting) +
  geom_line(aes(x=season, y=h, colour=as.factor(groupno)))

# Calculate mean h
mean_h <- underreporting %>% 
  group_by(label, agegrp, groupno, agemid, age_low_y) %>% 
  dplyr::summarise(cases_exp = mean(cases_exp),
                   hmin = min(h), 
                   hmax = max(h),
                   hmean = mean(h)) %>% 
  dplyr::ungroup() %>% arrange(groupno)

epidata_25 <- merge(epidata_25, mean_h[,c("groupno", "hmean")], by="groupno", all=T)
epidata_25 <- epidata_25 %>% arrange(groupno)

write.csv(epidata_25, "data/epidata_25.csv", row.names = FALSE)


# Detection probabilities q -------------------------------------------------

# Calculate min/max estimated number of symptomatic cases

detect_obs <- AGI %>% 
  dplyr::group_by(season, agegrp) %>% 
  dplyr::summarise(cases_obs = sum(npos, na.rm=T))

detect_exp <- underreporting %>% 
  mutate(agegrp = cut(age_low_y, breaks=c(0,2,5,15,35,99), right=FALSE)) %>% 
  dplyr::group_by(agegrp, season) %>% 
  dplyr::summarise(cases_exp = sum(cases_exp)) %>% 
  merge(., detect_obs, by=c("agegrp", "season"), all=T) %>% dplyr::ungroup() %>% 
  dplyr::group_by(season) %>% dplyr::arrange(agegrp) %>% 
  dplyr::mutate(groupno = 1:n()) %>% dplyr::ungroup() %>% 
  dplyr::mutate(q = cases_obs/cases_exp)


# Calculate mean q
mean_q <- detect_exp %>% 
  group_by(agegrp, groupno) %>% 
  dplyr::summarise(cases_exp = mean(cases_exp),
                   qmin = min(q), 
                   qmax = max(q),
                   qmean = mean(q))

# save q values for use in ODE model
write.csv(mean_q[,c("groupno", "qmean")], "data/qmean_AGI.csv", row.names = FALSE)




# Fig S2 --------------------------------------------------------------------

# Detection probability of outpatient cases

(figS2a <- ggplot(detect_exp) +
   ggtitle("A") +
   geom_boxplot(aes(x=agegrp, y=cases_obs, colour="observed")) +
   geom_boxplot(aes(x=agegrp, y=cases_exp, colour="expected")) +
   scale_y_log10() + 
   xlab("age group (years)") + ylab("cases per year") +
   theme(legend.title = element_blank(),
         legend.position = c(0.2, 0.5),
         legend.background = element_rect(fill="white", colour="black")))

(figS2b <- ggplot() +
    ggtitle("B") +
    scale_y_continuous(trans = "log10") +
    geom_point(data=detect_exp, 
               aes(x=groupno, y=q, shape=as.factor(season))) +
    geom_line(data=mean_q, 
              aes(x=groupno, y=qmean, colour="mean q")) +
    geom_line(data=mean_q, 
              aes(x=groupno, y=qmean*1.5, colour="mean q + 50%")) +
    geom_line(data=mean_q, 
              aes(x=groupno, y=qmean*0.5, colour="mean q - 50%")) +
    scale_x_continuous(breaks = 1:5, labels = mean_q$agegrp) +
    ylab("proportion detected (q)") + xlab("age group (years)") +
    labs(colour=NULL) +
    guides(shape="none") +
    theme(legend.position = c(0.8, 0.8),
          #legend.text = element_text(size=5),
          #legend.title = element_text(size=8),
          legend.box.background = element_rect(fill="white"))) # "bottom", legend.box="vertical"

figS2 <- grid.arrange(figS2a, figS2b, nrow=1)


ggsave(paste0("output/figures/figS2.", figform), 
       plot=figS2, 
       width = 24, height = 12, units="cm", dpi = dpi)





# Fig S3 -------------------------------------------------------------------

helper_figS3 <- data.frame("agemid" = seq(from=1, to=81, by=5))
helper_figS3$prediction <- predict(fit_ar, newdata=helper_figS3)

(figS3a <- ggplot(attackrates) +
   ggtitle("A") +
   geom_point(aes(x=agemid, y=ar, colour=as.factor(ref))) +
   ylab("Attack rate") + xlab("age (years)") +
    scale_y_continuous(labels = scales::percent) +
   geom_line(data=helper_figS3, aes(x=agemid, y=prediction)) +
   labs(colour="data from") +
   theme(legend.position = c(0.8, 0.9),
         legend.background = element_rect(fill="white", colour="black")))

(figS3b <- ggplot(epidata_25) +
    ggtitle("B") +
    geom_bar(aes(x=reorder(label, groupno), y=ar*100), stat="identity", fill="grey") +
    ylab("attack rate (%)") + xlab("age group (years)") +
    theme(axis.text.x = element_text(angle=90, hjust=1)))

figS3 <- grid.arrange(figS3a, figS3b, nrow=1)

ggsave(paste0("output/figures/figS3.", figform), 
       plot=figS3, 
       width = 24, height = 12, units="cm", dpi = dpi)


# Fig S4 -------------------------------------------------------------------


(figS4a <- ggplot(asymptomatic) +
   ggtitle("A") +
   geom_linerange(aes(xmin=age_low_y, xmax=age_up_y, y=pa*100, colour=as.factor(ref_pa))) +
   geom_point(aes(x=agemid, y=pa*100, colour=as.factor(ref_pa))) +
   geom_errorbar(aes(x=agemid, ymin=pa_low95*100, ymax=pa_up95*100, colour=as.factor(ref_pa))) +
   ylab("asymptomatics (%)") + xlab("age (years)") +
   geom_line(aes(x=agemid, y=pa_predicted*100)) +
   labs(colour="data from") +
   #theme(legend.position = "bottom") 
   theme(legend.position = c(0.8, 0.9),
         legend.background = element_rect(fill="white", colour="black")))


(figS4b <- ggplot(epidata_25) +
    ggtitle("B") +
    geom_bar(aes(x=reorder(label, groupno), y=pa*100), stat="identity", fill="magenta") +
    ylab("asymptomatics (%)") + xlab("age group (years)") +
    theme(axis.text.x = element_text(angle=90, hjust=1)))


figS4 <- grid.arrange(figS4a, figS4b, nrow=1)

ggsave(paste0("output/figures/figS4.", figform), 
       plot=figS4, 
       width = 24, height = 12, units="cm", dpi = dpi)


# Fig S5 -----------------------------------------------------------------------

# Probability of hospitalised and reported cases

(figS5a <- ggplot(underreporting) +
   ggtitle("A") +
   geom_boxplot(aes(x=reorder(label, groupno), y=cases_obs, colour="observed")) +
   geom_boxplot(aes(x=reorder(label, groupno), y=cases_exp, colour="expected")) +
   geom_boxplot(aes(x=reorder(label, groupno), y=n, colour="population size")) +
   scale_y_log10() + 
   xlab(NULL) + ylab("cases per year") +
   theme(legend.title = element_blank(),
         legend.position = c(0.2, 0.8),
         legend.background = element_rect(fill="white", colour="black"),
         axis.text.x = element_text(angle=90, hjust=1)))

(figS5b <- ggplot() +
    ggtitle("B") +
    geom_point(data=underreporting[underreporting$groupno>12,], 
               aes(x=agemid, y=h, shape=as.factor(season))) +
    scale_y_continuous(trans = "log10") +
    geom_line(data=mean_h[mean_h$groupno>12,], 
              aes(x=agemid, y=hmean, colour="mean h")) +
    geom_line(data=mean_h[mean_h$groupno>12,], 
              aes(x=agemid, y=hmean*1.5, colour="mean h+ 50%")) +
    geom_line(data=mean_h[mean_h$groupno>12,], 
              aes(x=agemid, y=hmean*0.5, colour="mean h - 50%")) +
    #scale_y_continuous(trans = "log") +
    ylab("proportion hospitalised (h)") + xlab("age (years)") +
    labs(colour=NULL) +
    guides(shape="none") +
    theme(legend.position = c(0.3, 0.8),
          #legend.text = element_text(size=5),
          #legend.title = element_text(size=8),
          legend.box.background = element_rect(fill="white"))) # "bottom", legend.box="vertical"

figS5 <- grid.arrange(figS5a, figS5b, nrow=1)

ggsave(paste0("output/figures/figS5.", figform), 
       plot=figS5, 
       width = 24, height = 12, units="cm", dpi = dpi)



# Temperature data -----------------------------------------

stationid <- c("02315", "02201", "03545", "01207", "04036", "03591", "03730", "03811")
url <- "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/"
snippet <- "tageswerte_KL_"

filenames <- read_html(url) %>% 
            html_elements("body") %>%
            html_text2() %>% str_split(pattern = "[\r\n]") %>% 
            unlist() %>% str_split(pattern = " ") %>% unlist()
filenames <- filenames[grepl(snippet, filenames)]


data <- data.frame(MESS_DATUM = numeric(), TMK = numeric(), stationid = numeric())
for (i in 1:length(stationid)) {
  
  print(stationid[which(stationid==i)])
  foo <- filenames[grep(paste0(snippet, stationid[i]), filenames)]
  j <- which(filenames==foo)
    print(foo)
  temp <- tempfile()
  download.file(paste0(url, foo), temp)
    print(fs::path(temp))
  moo <- paste0("produkt_klima_tag_", str_split(filenames[j], "_")[[1]][4],"_", str_split(filenames[j], "_")[[1]][5], "_", stationid[i], ".txt")
    print(moo)
  out <- fread(unzip(fs::path(temp), files = moo))
  rm(temp)
  out <- out %>% dplyr::filter(MESS_DATUM >=20000101) %>% 
                  dplyr::select(MESS_DATUM, TMK) %>% dplyr::mutate(stationid = i)
  print(out)
  data <- rbind(data, out)
  
}
rm(temp)
rm(out)
data$date <- as.Date(as.character(data$MESS_DATUM), "%Y%m%d")
data <- data %>% filter(TMK!=-999)
data$doy <- yday(data$date)
data <- data %>% filter(doy!=366)

ggplot(data) +
  geom_line(aes(x=date, y=TMK, colour=as.factor(stationid)))

temp <- data %>% 
  dplyr::group_by(doy) %>% 
  dplyr::summarise(temp = mean(TMK))

# Fit sinusoidal model to find the minimum
tempfit <- lm(TMK ~ sin(2*pi*doy/365) + cos(2*pi*doy/365), data=data)
summary(tempfit)

data$fit <- predict(tempfit)

ggplot(data) +
  geom_line(aes(x=doy, y=TMK, colour=as.factor(stationid))) +
  geom_line(aes(x=doy, y=fit), colour="red") +
  geom_line(data=temp, aes(x=doy, y=temp))

# Doy when the temperature is at the minimum
data$doy[data$fit==min(data$fit)][1]
(data$doy[data$fit==min(data$fit)][1] + 180) / 364



# Adjustment of trial VE estimates for use in model with Erlang distributed waning ---------------

#' Functions
#' waning in the model expressed as 1-CDF and scaled by VE at vaccination.
#' Mean of Erlang-k distribution is k/T
#' @param VE0 protection among protected
#' @param T parameter of Erlang distribution
#' @param k parameter of Erlang distribution
erlang.decay = function(VE0, T, k=2, t=1:150) {
  res = 0
  for(n in 0:(k-1)){
    res = res + 1/factorial(n)*exp(-T*t)*(T*t)^n }
  return( VE0 * res) 
}


#' target function for optimization
#' Computes the squared difference of the assumed waning (Erlang-distribution) multiplied 
#' with endpoint-specific VE parameters and observed VE during T observation periods.
#' VE parameters are constrained in that VE against severe disease/hospitalization >= VE against symptomatic disease.
#' @param par VE parameters. vector of length 3, first entry: VE against sympt. disease on logistic scale, 
#' second entry: (additional) VE against severe disease/hospitalization on logistic scale, 
#' third entry: duration of protection
#' @param eval_t_1 Vector of length T. start of (multiple) trial data periods in days
#' @param eval_t_2 Vector of length T. end of (multiple) trial data periods in days
#' @param obs_inc_sympt Vector of length T. Observed VE (1-RR) against incident symptomatic disease in (multiple) observation periods
#' @param obs_inc_sympt Vector of length T. Observed VE (1-RR) against severe disease/hospitalization in (multiple) observation periods
#' @param k degree of Erlang distribution (assumend waning of VE), default k=2
#'    

target_fun_combined_const = function(par=c(p1 = 0,
                                           p2 = 1,
                                           dur_prot = 200), 
                                     eval_t_1,
                                     eval_t_2, 
                                     obs_inc_sympt,
                                     obs_inc_hosp,
                                     k=2) {
  T_deriv = 1/(par[3]/k)
  v0_1 = plogis(par[1])
  v0_2 = plogis(par[1] + par[2])
  sq_err = 0
  for(i in 1:length(eval_t_1)) {
    sq_err = sq_err + (mean(erlang.decay(VE0 = v0_1, 
                                         T = T_deriv, 
                                         t=eval_t_1[i]:eval_t_2[i], 
                                         k=k)) - obs_inc_sympt[i])^2 + 
      (mean(erlang.decay(VE0 = v0_2, 
                         T = T_deriv, 
                         t=eval_t_1[i]:eval_t_2[i], 
                         k=k)) - obs_inc_hosp[i])^2
  }
  sq_err
}

# Maternal vaccine
#' Maternal protection (in infants)
#' Data Kampmann et al (2023) 
#' VE on incident endpoints during 4 time periods
#' VE against medically attended RSV-associated lower respiratory tract illness
#'(Figure 2B and Table S6)
#' and hospitalization (Table S7)
data_kamp = tibble(eval_t_1 = c(0, 91, 151, 181),
                   eval_t_2 = c(90, 150, 180, 360),
                   obs_inc_sympt = pmax(c(1- (24/3495)/(56/3480), 
                                          1-((35-24)/(3495-24))/((81-56)/(3480-56)),
                                          1-((47-35)/(3495-35))/((99-81)/(3480-81)),
                                          1-((92-47)/(3495-47))/((156-99)/(3480-99))),
                                        0),
                   obs_inc_hosp = pmax(c(1- (10/3495)/(31/3480), 
                                         1-((15-10)/(3495-10))/((37-31)/(3480-31)),
                                         1-((17-15)/(3495-15))/((39-37)/(3480-37)),
                                         1-((38-17)/(3495-17))/((57-39)/(3480-39))),
                                       0))


optim_res_kamp = optim(c(0.62, 0.6, 200), 
                       function(par) target_fun_combined_const(par, 
                                                               eval_t_1 = data_kamp$eval_t_1,
                                                               eval_t_2 = data_kamp$eval_t_2,
                                                               obs_inc_sympt = data_kamp$obs_inc_sympt,
                                                               obs_inc_hosp = data_kamp$obs_inc_hosp,
                                                               k=2), 
                       lower = c(-Inf, 0, 50), 
                       upper = c(Inf, Inf, 5000),
                       method = "L-BFGS-B")

# Extract VE Parameters
(round(c(VE_0_sympt = plogis(optim_res_kamp$par[1]),
         VE_0_hosp = plogis(optim_res_kamp$par[1] + optim_res_kamp$par[2]),
         mean_dur = round(optim_res_kamp$par[3])), 2))


# Nirsevimab
# Data

data_symp_inc = tibble(obs_inc_sympt = c(0.887, 0.530, 0.486),
                       eval_t_1 = c(0, 91, 121),
                       eval_t_2 = c(90, 120, 150))

data_hosp = tibble(eval_t_1 = 0,
                   eval_t_2 = 150,
                   obs_inc_hosp = 0.81)

#' Same target function but with different observation period for VE against 
#' symptomatic infection and hospitalization
target_fun_2periods = function(par=c(ve_sympt = 0,
                                     ve_hosp = 0.2,
                                     dur_prot = 150), 
                               eval_t_1_sympt,
                               eval_t_2_sympt, 
                               eval_t_1_hosp,
                               eval_t_2_hosp, 
                               obs_inc_sympt,
                               obs_inc_hosp,
                               k=2) {
  T_deriv = 1/(par[3]/k)
  v0_1 = plogis(par[1])
  v0_2 = plogis(par[1] + par[2])
  
  sq_err = 0
  for(i in 1:length(eval_t_1_sympt)) {
    sq_err = sq_err + (mean(erlang.decay(VE0 = v0_1, 
                                         T = T_deriv, 
                                         t=eval_t_1_sympt[i]:eval_t_2_sympt[i], 
                                         k=k)) - obs_inc_sympt[i])^2
  }
  for(i in 1:length(eval_t_1_hosp)) {
    sq_err = sq_err + (mean(erlang.decay(VE0 = v0_2, 
                                         T = T_deriv, 
                                         t=eval_t_1_hosp[i]:eval_t_2_hosp[i], 
                                         k=k)) - obs_inc_hosp[i])^2
  }
  sq_err
}


optim_res_nirs = optim(c(1, 1, 200), 
                       function(par) target_fun_2periods(par, 
                                                         eval_t_1_sympt = data_symp_inc$eval_t_1,
                                                         eval_t_2_sympt = data_symp_inc$eval_t_2,
                                                         eval_t_1_hosp = data_hosp$eval_t_1,
                                                         eval_t_2_hosp = data_hosp$eval_t_2,
                                                         obs_inc_sympt = data_symp_inc$obs_inc_sympt,
                                                         obs_inc_hosp = data_hosp$obs_inc_hosp,
                                                         k=2), 
                       lower = c(-Inf, 0, 50), 
                       upper = c(Inf, Inf, 5000),
                       method = "L-BFGS-B")

(round(c(VE_0_sympt = plogis(optim_res_nirs$par[1]),
         VE_0_hosp = plogis(optim_res_nirs$par[1] + optim_res_nirs$par[2]),
         mean_dur = round(optim_res_nirs$par[3])), 2))


# Plot
params = tibble(type = c("Nirs.", 
                         "Pali.", 
                         "Maternal vacc.: infant", 
                         "Maternal vacc.: adult",
                         "Older adult vacc."),
                VE_0_sympt = c(0.96, 0.82, 0.75, 0.67, 0.77),
                VE_0_hosp = c(1.00, 0.95, .75, .92, 1.00),
                mean_dur = c(160, 40, 133, 712, 712)
)

prot_dat = expand.grid(t=0:730,
                       type = c("Nirs.", 
                                "Pali.",
                                "Maternal vacc.: infant",
                                "Maternal vacc.: adult",
                                "Older adult vacc.")) %>%
  left_join(params) %>%
  mutate(t_calc = if_else(type=="Pali.", t-rep(c(0, 30, 60, 90, 120), times=c(31, 30, 30, 30, 730-120)),
                          t)) %>%
  mutate(prot_sympt = erlang.decay(VE0 = VE_0_sympt, T=2/mean_dur, k=2, t=t_calc),
         prot_hosp = erlang.decay(VE0 = VE_0_hosp, T=2/mean_dur, k=2, t=t_calc))


# Fig. S6 -------------------------------------

(figS6 <-  prot_dat %>% pivot_longer(cols = c("prot_sympt", "prot_hosp")) %>%
  mutate(name = factor(name, levels = c("prot_sympt", "prot_hosp"),
                       labels = c("Symptomatic", "Hospitalisation"))) %>%
  mutate(type = factor(type, levels = c("Pali.",
                                        "Nirs.", 
                                        "Maternal vacc.: infant",
                                        "Maternal vacc.: adult",
                                        "Older adult vacc."))) %>%
  ggplot() + 
  theme_bw() +
  geom_line(aes(t, value, col = type, lty = type), linewidth=.8) +
  facet_grid(rows=vars(name)) +
  scale_color_brewer(type="qual", palette = 2) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  scale_linetype_manual(values = c("Pali."=1,
                                   "Nirs."=1, 
                                   "Maternal vacc.: infant"=1,
                                   "Maternal vacc.: adult"=3,
                                   "Older adult vacc."=3)) +
  ylab("Vaccine efficacy") +
  xlab("Time after immunisation (days)") +
  theme(legend.position = "right"))

ggsave(plot = figS6, 
       filename = paste0("output/figures/figS6.", figform), 
       width = 9, height = 4,
       dpi = dpi)


# Estimation of scaling factor for total hosp and ICU patients ------------------


#' Estimate the age-specific rate of ICU admissions and fatalities among
#' hospitalized based on INEK data from 2019 to 2023 (RSV-specific cases)
#' for simulation of severe disease endpoints

inek = read_delim("data/InEK-Daten_RSV-spezifisch_2019-2023.csv", 
                  delim = ";", 
                  locale = locale(decimal_mark = ","))

inek_dat = inek %>% select(year, agegroup, value, Parameter) %>%
  pivot_wider(names_from = Parameter, values_from = value) %>%
  mutate(agegroup=as.factor(agegroup))

inek_smry = inek %>% select(year, agegroup, value, Parameter) %>%
  pivot_wider(names_from = Parameter, values_from = value) %>%
  group_by(agegroup) %>%
  summarise(hosp_icu_ratio = mean(ICU/hosp),
            hosp_mort_ratio = mean(mort/hosp))

#' Estimate ICU rates per agegroup based on GLM with 
#' Poisson link and log-hospitalization counts as offset
mod_icu = glm(ICU ~ agegroup, offset = log_hosp, 
              family = poisson(link = "log"), 
              data = inek_dat %>% mutate(log_hosp = log(hosp)))
#' Estimate Mort rates per agegroup based on GLM with 
#' Poisson link and log-hospitalization counts as offset
mod_mort = glm(mort ~ agegroup, offset = log_hosp, 
               family = poisson(link = "log"), 
               data = inek_dat %>% mutate(log_hosp = log(hosp)))


# Create summary table/tibble for use in post-processing of simulation results 
inek_smry$pred_icu_link = predict(mod_icu, newdata = inek_smry %>% 
                                    select(agegroup) %>% 
                                    mutate(agegroup = as.factor(agegroup),
                                           log_hosp=0))
inek_smry$pred_icu_link_se = predict(mod_icu, newdata = inek_smry %>% 
                                       select(agegroup) %>% 
                                       mutate(agegroup = as.factor(agegroup),
                                              log_hosp=0), se.fit = T)$se.fit

inek_smry$pred_mort_link = predict(mod_mort, newdata = inek_smry %>% 
                                     select(agegroup) %>% 
                                     mutate(agegroup = as.factor(agegroup),
                                            log_hosp=0))
inek_smry$pred_mort_link_se = predict(mod_mort, newdata = inek_smry %>% 
                                        select(agegroup) %>% 
                                        mutate(agegroup = as.factor(agegroup),
                                               log_hosp=0), se.fit = T)$se.fit
inek_smry
# Join INEK data to 25 age-groups of the DTM
inek_smry_for_postproc = tibble(age = 1:25) %>%
  mutate(agegroup = case_when(age %in% 1:12 ~ 0,
                              age %in% 13:14  ~ 1,
                              age %in% 15:16 ~ 3,
                              age %in% 17 ~ 6,
                              age %in% 18 ~ 10,
                              age %in% 19 ~ 16,
                              age %in% 20 ~ 30,
                              age %in% 21 ~ 40,
                              age %in% 22 ~ 50,
                              age %in% 23 ~ 55,
                              age %in% 24 ~ 65,
                              age %in% 25 ~ 75)) %>%
  left_join(inek_smry, by = "agegroup") %>%
  select(-agegroup)

write_csv2(inek_smry_for_postproc, file = "data/INEK_derived_endpoint_rates.csv")


# Fig. S7 -----------------------------------

figS7 <- inek_smry_for_postproc %>%
  mutate(mean = exp(pred_icu_link),
         q025 = exp(qnorm(0.025, pred_icu_link, pred_icu_link_se)),
         q975 = exp(qnorm(0.975, pred_icu_link, pred_icu_link_se)),
         type = "Hosp-ICU") %>%
  rbind(inek_smry_for_postproc %>%
          mutate(mean = exp(pred_mort_link),
                 q025 = exp(qnorm(0.025, pred_mort_link, pred_mort_link_se)),
                 q975 = exp(qnorm(0.975, pred_mort_link, pred_mort_link_se)),
                 type = "Hosp-Mort")) %>%
  mutate(age_start = case_when(age == 1 ~ 0,
                               age == 2 ~ 1/12, 
                               age == 3~ 2/12, 
                               age == 4~ 3/12, 
                               age == 5~ 4/12, 
                               age == 6~ 5/12, 
                               age == 7~ 6/12, 
                               age == 8~ 7/12, 
                               age == 9~ 8/12, 
                               age == 10~9/12, 
                               age == 11~10/12,
                               age == 12~11/12, 
                               age == 13~1, 
                               age == 14~2, 
                               age == 15~3, 
                               age == 16~4, 
                               age == 17~5, 
                               age == 18~10, 
                               age == 19~15, 
                               age == 20~25, 
                               age == 21~35, 
                               age == 22~45, 
                               age == 23~55, 
                               age == 24~65,
                               age == 25~75)) %>%
  arrange(age) %>%
  group_by(type) %>%
  mutate(age_end = replace_na(lead(age_start, 1), 100)) %>%
  select(age, age_start, age_end, mean, q025, q975, type) %>%
  ggplot() + 
  geom_segment(aes(x=age_start, xend = age_end, y= mean, yend = mean,
                   col = type), lwd = 1) +
  geom_segment(aes(x=age_start, xend = age_end, y= q025, yend = q025,
                   col = type)) +
  geom_segment(aes(x=age_start, xend = age_end, y= q975, yend = q975,
                   col = type)) +
  geom_segment(aes(x=age_end, xend = age_end, y= q025, yend = q975,
                   col = type)) +
  geom_segment(aes(x=age_start, xend = age_start, y= q025, yend = q975,
                   col = type)) +
  xlab("Age") +
  scale_x_continuous(minor_breaks = seq(0,100,5)) +
  ylab("Ratio") +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(legend.title = element_blank())

ggsave(plot = figS7, 
       filename = paste0("output/figures/figS7.", figform), 
       width = 9, height = 4,
       dpi = dpi)


