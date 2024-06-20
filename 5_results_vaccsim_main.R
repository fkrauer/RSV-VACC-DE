# Housekeeping  -------------------------------------------------
# Authors: Felix Guenther, Fabienne Krauer
#---


library(tidyverse)
theme_set(theme_bw())

source("func_helpers.R")

set.seed(240605122)

# fig format
figform <- "png"
dpi <- 600

inc_plot = c("symptomatic",
             "hospitalised",
             "icu",
             "mort")


# Read data -----------------------------------

inc_orig = read_csv('output/vacc_simulations/vaccsim_inc_0_vs_[1,2,3].csv') %>% 
      full_join(read_csv('output/vacc_simulations/vaccsim_inc_2_vs_[4,5,6].csv')) %>% 
      label_inc()
 
inc_incr_orig <- read_csv('output/vacc_simulations/vaccsim_inc_4_vs_[6].csv') %>% 
  full_join(read_csv('output/vacc_simulations/vaccsim_inc_5_vs_[4].csv')) %>% 
  label_inc() %>% 
  full_join(inc_orig %>% filter(strategy == 5))


immu = read_csv('output/vacc_simulations/vaccsim_immunisations_0_vs_[1,2,3].csv') %>% 
  full_join(read_csv('output/vacc_simulations/vaccsim_immunisations_2_vs_[4,5,6].csv'))

immu_incr = read_csv('output/vacc_simulations/vaccsim_immunisations_4_vs_[6].csv') %>% 
  full_join(read_csv('output/vacc_simulations/vaccsim_immunisations_5_vs_[4].csv')) %>% 
  full_join(immu %>% filter(strategy == 2))


pop_dat = readRDS("data/pop_DE_25_2015_2019.rds") %>% tibble() %>%
  filter(year==2019) %>%
  select(age=groupno, n_pop=n)


inek_smry = read_csv2("data/INEK_derived_endpoint_rates.csv")


# Calculate ICU and mortality incidence ----------------------------------------

inc_extend <- calc_icu_deaths(inc_orig, inek_smry)

inc_incr_extend <- calc_icu_deaths(inc_incr_orig, inek_smry)

# Upscale older adults hospitalisations  ----------------------------------------

inc <- inc_extend %>% upscale_oa(scaling_fac = 8)


# Analysis -----------------------------

## Calculate Number of averted total cases (over all age groups) -------------------------------------------------

# Calculate 5-year total incidence and averted cases per strategy, endpoint and MCMC iteration
inc_per_strat = inc %>% 
  group_by(incidence, strategy, seasonal, replicate) %>% 
  summarise(n_inc = sum(n),
            n_base = sum(n_base),
            n_avert = sum(n_avert)) %>%
  ungroup()


# Calculate yearly overall incidence and averted cases per strategy, endpoint and MCMC iteration
inc_per_strat_yearly = inc %>% 
  group_by(incidence, strategy, seasonal, replicate, year) %>% 
  summarise(n_inc = sum(n),
            n_base = sum(n_base),
            n_avert = sum(n_avert),
            prop_avert = n_avert/n_base) %>% ungroup() %>%
  group_by(incidence, strategy, seasonal, replicate) %>%
  summarise(n_inc = mean(n_inc),
            n_base = mean(n_base),
            n_avert = mean(n_avert),
            prop_avert = n_avert / n_base) %>%
  ungroup()



## Calculate overall N of immunised individuals (over all age groups) ----------------------------------

# Derive total/yearly average number immunised per strategy and MCMC iteration

# Over 5 years
immu_per_strat = immu %>% 
  group_by(strategy, seasonal, replicate) %>%
  summarise(n_immu = sum(n),
            n_base = sum(n_base),
            n_immu_add = n_immu - n_base) %>% # additionally immunised
  ungroup()


# per year
immu_per_strat_yearly = immu %>% 
  group_by(strategy, seasonal, replicate, year) %>% 
  summarise(n_immu = sum(n),
            n_base = sum(n_base),
            n_immu_add = n_immu - n_base) %>% ungroup() %>%
  group_by(strategy, seasonal, replicate) %>%
  summarise(n_immu = mean(n_immu),
            n_base = mean(n_base),
            n_immu_add = mean(n_immu_add)) %>%
  ungroup()



## Infant strategies (1-3) ------------------------------------------------

inf_labels = c("Nirs. HR",
               "Nirs. all",
               "Mat.\nvacc.")


# Yearly, Absolute
case_prev_inf = inc_per_strat_yearly %>% 
  filter(seasonal == TRUE, strategy %in% 1:3) %>%
  group_by(incidence, strategy) %>%
  summarise(n_avert_med = median(n_avert),
            n_avert_q025 = quantile(n_avert, 0.025),
            n_avert_q975 = quantile(n_avert, 0.975),
            pct_avert_med = round(median(prop_avert),4)*100,
            pct_avert_q025 = round(quantile(prop_avert, 0.025),4)*100,
            pct_avert_q975 = round(quantile(prop_avert, 0.975),4)*100) %>%
  mutate(strategy = factor(strategy, levels = 3:1, labels = rev(inf_labels)))


# Yearly, absolute and percent by age groups
res_tab = do.call(rbind,
                  lapply(list(2:6, 7:12, 13, 14, 15, 1:12, 1:13), function(x) {
                    perc_avert_yearly(age_groups = x)
                  })) %>%
  mutate(age_cat = factor(age_cat, 
                          levels = unlist(lapply(list(2:6, 7:12, 13, 14, 15, 1:12, 1:13), 
                                                 function(x) paste0(x, collapse=""))),
                          labels = c("1-5m", "6-11m", "1y", "2y", "3-5y", "0-11m", "0-1y")))

# NNI in infants 

add_immu_per_strat_year_posterior_infants = immu %>%
  filter(strategy %in% 1:3, seasonal==TRUE) %>%
group_by(strategy, seasonal, replicate, year) %>%
  summarise(n_immu = sum(n),
            n_base = sum(n_base),
            n_immu_add = n_immu - n_base) %>% 
  ungroup() %>% group_by(strategy, seasonal, replicate) %>%
  summarise(n_immu = mean(n_immu),
            n_base = mean(n_base),
            n_immu_add = mean(n_immu_add)) %>%
  ungroup() %>% group_by(strategy, seasonal) %>%
  summarise(n_immu_add_year_postmed = quantile(n_immu_add, 0.5),
            n_immu_add_year_q025 = quantile(n_immu_add, 0.025),
            n_immu_add_year_q975 = quantile(n_immu_add, 0.975)) %>% ungroup()


nni_inf = inc_per_strat %>% filter(seasonal == TRUE,
                                  strategy %in% 1:3) %>%
  left_join(immu_per_strat %>% select(-n_base)) %>%
  mutate(nni = n_immu_add/n_avert) %>%
  group_by(incidence, strategy) %>%
  summarise(nni_med = median(nni),
            nni_q025 = quantile(nni, 0.025),
            nni_q975 = quantile(nni, 0.975)) %>%
  mutate(nni_med = ifelse(strategy == 1, 0, nni_med),
         nni_q025 = ifelse(strategy == 1, 0, nni_q025),
         nni_q975 = ifelse(strategy == 1, 0, nni_q975)) %>%
  mutate(strategy = factor(strategy, levels = 3:1, labels = rev(inf_labels)))

## Table S14 -------------------------------------------

# Annual incidences by age and total for all outcomes in the current strategy (0):

cases_y_all <- inc_extend %>% 
  filter(year==1 & strategy == 1) %>% 
  group_by(incidence, replicate) %>% 
  summarise(n_base = sum(n_base)) %>%
  ungroup() %>% group_by(incidence) %>% 
  summarise(inc_median = format(round(median(n_base),0), big.mark = ","),
            inc_low95PPI = format(round(quantile(n_base, 0.025),0), big.mark = ","),
            inc_up95PPI = format(round(quantile(n_base, 0.975),0), big.mark = ","),
            cases = paste0(inc_median, " [", inc_low95PPI, ", ", inc_up95PPI, "]")) %>% 
  ungroup() %>% dplyr::mutate(age = "total") %>% 
  pivot_wider(id_cols = age, names_from = incidence, values_from = cases)

cases_y_age <- inc_extend %>% 
  filter(year==1 & strategy == 1) %>% 
  dplyr::group_by(age, incidence) %>% 
  dplyr::summarise(inc_median = format(round(median(n_base),0), big.mark = ","),
                   inc_low95PPI = format(round(quantile(n_base, 0.025),0), big.mark = ","),
                   inc_up95PPI = format(round(quantile(n_base, 0.975),0), big.mark = ","),
                   cases = paste0(inc_median, " [", inc_low95PPI, ", ", inc_up95PPI, "]")) %>% 
  pivot_wider(id_cols = "age", names_from = incidence, values_from = cases) %>% 
  dplyr::mutate(age = as.character(age))

tableS14 <- rbind(cases_y_age, cases_y_all) %>% 
  relocate(c(total, symptomatic), .after=age)


tableS14


### Fig 3 -----------------------------------

fig3a = case_prev_inf %>% 
  filter(incidence %in% inc_plot) %>%
  mutate(incidence = factor(incidence, levels = inc_plot)) %>%
  mutate(inc_avert_print = round_print(n_avert_med)) %>%
  ggplot() + 
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=14),
        strip.text = element_text(size=14)) +
  labs(title = "Overall number of cases prevented", tag = "A") +
  geom_col(aes(strategy, n_avert_med), 
           fill = "lightgreen", 
           alpha = .5, col = "black") +
  geom_label(aes(strategy, 
                 n_avert_med,
                 label = label_form(n_avert_med),
                 hjust=ifelse(strategy=="Nirs. HR" |
                                 (strategy=="Mat.\nvacc." & 
                                    incidence %in% c("symptomatic", "hospitalised", "icu")), 
                               0, 1)), fill = "lightgrey",
             vjust=1, nudge_x = .45, size=5
             ) +
  geom_linerange(aes(strategy, ymin = n_avert_q025, 
                     ymax=n_avert_q975)) +
  facet_wrap(~incidence, scales = "free_x", ncol = length(inc_plot)) + 
  #ylab("Prevented cases") +
  ylab(NULL) +
  scale_x_discrete(guide = guide_axis(angle = 0)) + xlab("") +
  scale_y_continuous(n.breaks = 3) +
  scale_fill_brewer(type = "qual") +
  coord_flip() +
  theme(legend.position = "none")

fig3b = case_prev_inf %>% 
  filter(incidence %in% inc_plot) %>%
  mutate(incidence = factor(incidence, levels = inc_plot)) %>%
  ggplot() + 
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=14),
        strip.text = element_text(size=14)) +
  labs(title = "Overall percentage of cases prevented", tag = "B") +
  geom_col(aes(strategy, pct_avert_med), 
           fill = "lightpink", 
           alpha = .5, col = "black") +
  geom_label(aes(strategy, 
                 pct_avert_med,
                 label = pct_avert_med,
                 hjust=ifelse(strategy=="Nirs. HR" |
                                (strategy=="Mat.\nvacc." & 
                                   incidence %in% c("symptomatic", "hospitalised", "icu")), 
                              0, 1)), fill = "lightgrey",
             vjust=1, nudge_x = .45, size=5
  ) +
  geom_linerange(aes(strategy, ymin = pct_avert_q025, 
                     ymax=pct_avert_q975)) +
  facet_wrap(~incidence, scales = "free_x", ncol = length(inc_plot)) + 
  #ylab("Prevented cases") +
  ylab(NULL) +
  scale_x_discrete(guide = guide_axis(angle = 0)) + xlab("") +
  scale_y_continuous(n.breaks = 3) +
  scale_fill_brewer(type = "qual") +
  coord_flip() +
  theme(legend.position = "none")


fig3c = nni_inf %>%  
  mutate(nni_med_print = round_print(nni_med)) %>%
  filter(incidence %in% inc_plot) %>%
  mutate(incidence = factor(incidence, levels = inc_plot)) %>%
 ggplot() + 
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=14),
        strip.text = element_text(size=14)) +
  labs(title = "Number needed to vaccinate (NNV)", tag = "C") +
  geom_col(aes(strategy, nni_med), fill = "lightblue", alpha = .5, col = "black") +
  geom_label(aes(strategy, nni_med, 
                 label = label_form(nni_med),
                 hjust=ifelse(strategy=="Nirs. HR", 0, 1)),
             fill = "lightgrey", vjust=1, nudge_x = .45, size=5) +
  geom_linerange(aes(strategy, ymin = nni_q025, 
                     ymax=nni_q975)) +
  facet_wrap(~incidence, scales = "free_x", ncol = length(inc_plot)) + 
  #ylab("NNI") +
  ylab(NULL) +
  scale_x_discrete(guide = guide_axis(angle = 0)) + xlab("") +
  scale_fill_brewer(type = "qual") +
  scale_y_continuous(n.breaks = 3) +
  coord_flip() +
  theme(legend.position = "none")


fig3 = grid.arrange(fig3a, fig3b, fig3c, nrow = 3)

ggsave(fig3, 
       filename = paste0("output/figures/fig3.", figform), 
       width = 12, height = 7,
       dpi = dpi)

res_tab %>% filter(age_cat == "0-11m" & type == "rel_redu") %>% select(age_cat, strategy, hospitalised, icu)

case_prev_inf %>% filter(strategy == "Nirs. HR")
case_prev_inf %>% filter(strategy == "Nirs. all")
nni_inf %>% filter(strategy == "Nirs. all")

case_prev_inf %>% filter(strategy == "Mat.\nvacc.")
nni_inf %>% filter(strategy == "Mat.\nvacc.")


### Fig. S12 -------------------------------------------------------

# age-specific numbers of cases and incidences unter strategy 0

data_figS12 <- inc %>% filter(strategy==1,
                             seasonal==T,
                             incidence != "total") %>%
  full_join(pop_dat) %>%
  group_by(age, age_label, year, incidence, strategy) %>%
  summarise(n_med = quantile(n_base, 0.5),
            n_q025 = quantile(n_base, 0.025),
            n_q975 = quantile(n_base, 0.975),
            inc_med = quantile(n_base/n_pop*100000, 0.5),
            inc_q025 = quantile(n_base/n_pop*100000, 0.025),
            inc_q975 = quantile(n_base/n_pop*100000, 0.975)) %>%
  mutate(incidence = factor(incidence, levels = c("symptomatic", "hospitalised",
                                                  "icu", "mort"),
                            labels = c("sympt", "hosp", "icu", "mort")),
         year=paste0("Year: ", year))

(figS12a <- ggplot(data_figS12) +
  geom_col(aes(x=reorder(age_label, age), y=n_med)) +
  geom_linerange(aes(age, ymin = n_q025, ymax=n_q975)) +
  facet_grid(rows=vars(incidence),
             cols = vars(year), scales = "free_y") +
  ylab("N") +
  labs(title = "Strategy 0: annual cases", tag ="A") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks = age_breaks))


(figS12b <- ggplot(data_figS12) +
    geom_col(aes(x=reorder(age_label, age), y=inc_med)) +
    geom_linerange(aes(age, ymin = inc_q025, ymax=inc_q975)) +
    facet_grid(rows=vars(incidence),
               cols = vars(year), scales = "free_y") +
    ylab("N per 100,000") +
    labs(title = "Strategy 0: incidence", tag ="B") +
    xlab(NULL) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_x_discrete(breaks = age_breaks))


figS12 = grid.arrange(figS12a, figS12b, nrow = 2)

ggsave(figS12, 
       filename = paste0("output/figures/figS12.", figform), 
       width = 10, height = 8, dpi = dpi)



### Fig. S13 -------------------------------------------------------

# age-specific numbers of cases and incidences unter strategy 1

figS13 <- plot_incid_prevented(inc, strat=1, strat_name = "Nirs. HR")
ggsave(figS13, 
       filename = paste0("output/figures/figS13.", figform), 
       width = 10, height = 12.5, dpi = dpi)

### Fig. S14 -------------------------------------------------------

# age-specific numbers of cases and incidences unter strategy 2

figS14 <- plot_incid_prevented(inc, strat=2, strat_name = "Nirs. all")
ggsave(figS14, 
       filename = paste0("output/figures/figS14.", figform),
       width = 10, height = 12.5, dpi = dpi)



### Fig. S15 -------------------------------------------------------

# age-specific numbers of cases and incidences unter strategy 3

figS15 <- plot_incid_prevented(inc, strat=3, strat_name = "Mat. vacc")
ggsave(figS15, 
       filename = paste0("output/figures/figS15.", figform),
       width = 10, height = 12.5, dpi = dpi)


### Fig. S16 ------------------------------

# Proportion of prevented cases by simulation year per strategy
figS16 = inc %>%
  filter(strategy %in% c(1:6),
         seasonal==T) %>%
  group_by(year, incidence, strategy, replicate) %>%
  summarise(n=sum(n),
            n_base=sum(n_base),
            avert = n_base-n,
            rel_avert = avert/n_base) %>%
  group_by(year, incidence, strategy) %>%
  summarise(n=mean(n),
            n_base=mean(n_base),
            avert = mean(n_base-n),
            rel_avert = avert/n_base) %>%
  mutate(incidence = factor(incidence, levels = c("total", "symptomatic", "hospitalised",
                                                  "icu", "mort")),
         strategy = factor(strategy, levels = c(1:3, 5, 4, 6),
                           labels = c("Nirs HR", "Nirs all",
                                      "Mat. vacc",
                                      "OA (75+)",
                                      "OA (65+)",
                                      "OA (55+)"))) %>%
  group_by(incidence, strategy) %>%
  mutate(avert_share = avert/sum(avert)) %>%
  filter(incidence != "total") %>%
  mutate(year=factor(year, levels = 5:1)) %>%
  ggplot() + 
  geom_col(aes(strategy, avert_share, fill = year),
           col = "black") +
  facet_wrap(~incidence) +
  scale_y_continuous(breaks = seq(0,1, by=0.2)) +
  scale_fill_brewer(type="qual", palette = 1) +
  ylab("Share of prevented incidence")

ggsave(figS16, 
       filename = paste0("output/figures/figS16.", figform),
       width = 10, height = 5, dpi = dpi)




## Older adults (strategies 4-6) ------------------------------------------------


# Scale-up incidence by 1:15
inc_per_strat_oa_up <- do.call(
                      rbind,
                      lapply(c(1:15),
                             function(x) {
                               upscale_oa(inc_extend, x)
                             })) %>% 
                      group_by(incidence, strategy, seasonal, replicate, year, scaling_fac) %>% 
                      dplyr::summarise(n_inc = sum(n),
                                n_base = sum(n_base),
                                n_avert = sum(n_avert),
                                prop_avert = n_avert / n_base) %>%
                      ungroup()


case_prev_oa = inc_per_strat_oa_up %>% 
  filter(seasonal==T & strategy %in% 4:6) %>%
  group_by(incidence, strategy, scaling_fac, replicate) %>%
  summarise(n_avert = sum(n_avert)) %>%
  ungroup() %>% 
  group_by(incidence, strategy, scaling_fac) %>%
  summarise(n_avert_med = median(n_avert),
            n_avert_q025 = quantile(n_avert, 0.025),
            n_avert_q975 = quantile(n_avert, 0.975)) %>%
  dplyr::ungroup() %>% 
  mutate(strategy = factor(strategy, levels = c(5,4,6), 
                           labels = c("75+", "65+", "55+"))) %>%
  arrange(incidence, strategy, scaling_fac)

case_prev_oa_yearly = inc_per_strat_oa_up %>% 
  filter(seasonal==T & strategy %in% 4:6) %>%
  group_by(incidence, strategy, scaling_fac, replicate) %>%
  summarise(n_base = mean(n_base),
            n_avert = mean(n_avert),
            prop_avert = n_avert/n_base) %>%
  group_by(incidence, strategy, scaling_fac) %>%
  summarise(n_avert_med = median(n_avert),
            n_avert_q025 = quantile(n_avert, 0.025),
            n_avert_q975 = quantile(n_avert, 0.975),
            pct_avert_med = round(median(prop_avert),4)*100,
            pct_avert_q025 = round(quantile(prop_avert, 0.025),4)*100,
            pct_avert_q975 = round(quantile(prop_avert, 0.975),4)*100) %>%
  mutate(strategy = factor(strategy, levels = c(5,4,6), 
                           labels = c("75+", "65+", "55+"))) %>%
  arrange(incidence, strategy, scaling_fac)
 

# NNV
immu_per_strat_oa = immu_per_strat %>%
  filter(strategy %in% 4:6) %>%
  select(strategy, seasonal, replicate, n_immu) %>%
  left_join(immu_per_strat %>%
              filter(strategy %in% 2,
                     seasonal = TRUE) %>%
              select(replicate, n_base = n_immu)) %>%
  mutate(n_immu_add = n_immu-n_base)

nni_oa = inc_per_strat_oa_up %>%
  filter(seasonal==T & strategy %in% 4:6) %>%
  group_by(incidence, strategy, scaling_fac, replicate) %>%
  summarise(n_avert = sum(n_avert)) %>%
  dplyr::ungroup() %>% 
  left_join(immu_per_strat_oa %>% select(-n_base)) %>%
  mutate(nni = n_immu_add/n_avert) %>%
  group_by(incidence, strategy, scaling_fac) %>%
  summarise(nni_med = median(nni),
            nni_q025 = quantile(nni, 0.025),
            nni_q975 = quantile(nni, 0.975))  %>%
  dplyr::ungroup() %>% 
  mutate(strategy = factor(strategy, levels = c(5,4,6), 
                           labels = c("75+", "65+", "55+")))


### Fig 4 -----------------------------------

fig4a = case_prev_oa_yearly %>%
  filter(scaling_fac %in% c(8, 14)) %>%
  mutate(scaling_fac = factor(scaling_fac, levels = c(14, 8),
                              labels = paste0("Scaling: ", c(14, 8)))) %>%
  filter(incidence %in% inc_plot) %>%
  mutate(incidence = factor(incidence, levels = inc_plot)) %>%
  mutate(inc_avert_med_print = round_print(n_avert_med)) %>%
  ggplot() + 
  geom_col(aes(strategy, n_avert_med,
               fill = scaling_fac), 
           alpha = .5, col = "black", position = position_dodge2(0.5)) +
  geom_label(aes(strategy, n_avert_q025,
                 label = label_form(n_avert_med),
                 hjust=1,
                 group = scaling_fac),
             fill = "lightgrey", position = position_dodge2(0.9), size=4.5) +
  geom_linerange(aes(strategy, ymin = n_avert_q025, 
                     ymax=n_avert_q975, group = scaling_fac),
                 position = position_dodge2(0.9)) +
  facet_wrap(~incidence, scales = "free_x", ncol = length(inc_plot)) + 
  ylab(NULL) +
  labs(title = "Overall number of cases prevented") +
  scale_x_discrete(guide = guide_axis(angle = 0)) + xlab("") +
  scale_y_continuous(n.breaks = 3) +
  scale_fill_brewer(type = "qual") +
  coord_flip() +
  theme(legend.title = element_blank())



fig4b = nni_oa %>%  
  mutate(nni_med_print = round_print(nni_med)) %>%
  filter(scaling_fac %in% c(8, 14)) %>%
  mutate(scaling_fac = factor(scaling_fac, levels = c(14, 8),
                              labels = paste0("Scaling: ", c(14, 8)))) %>%
  filter(incidence %in% inc_plot) %>%
  mutate(incidence = factor(incidence, levels = inc_plot)) %>%
  ggplot() + 
  geom_col(aes(strategy, nni_med, fill = scaling_fac), 
           alpha = .5, col = "black", position = position_dodge2(0.5)) +
  geom_label(aes(strategy,
                 if_else(strategy == "75+" & incidence=="mort" & scaling_fac == "Scaling: 14", 
                         nni_q975, nni_q025),
                 label = label_form(nni_med),
                 hjust= if_else(strategy == "75+" & incidence=="mort" & scaling_fac == "Scaling: 14", 
                                0, 1),
                 group = scaling_fac),
             fill = "lightgrey", position = position_dodge2(.9), size=4.5) +
  geom_linerange(aes(strategy, ymin = nni_q025, 
                     ymax=nni_q975, group = scaling_fac),
                 position = position_dodge2(0.9)) +
  facet_wrap(~incidence, scales = "free_x", ncol = length(inc_plot)) + 
  ylab(NULL) +
  labs(title = "Number needed to vaccinate (NNV)") +
  scale_x_discrete(guide = guide_axis(angle = 0)) + xlab("") +
  scale_fill_brewer(type = "qual") +
  scale_y_continuous(n.breaks = 3) +
  coord_flip() + 
  theme(legend.title = element_blank())


fig4 = ggpubr::ggarrange(fig4a + theme(axis.text = element_text(size=13), 
                                                     axis.title = element_text(size=14),
                                                     strip.text = element_text(size=14)),
                           fig4b + theme(axis.text = element_text(size=13), 
                                              axis.title = element_text(size=14),
                                              strip.text = element_text(size=14)), 
                           nrow = 2,
                           labels = "AUTO", common.legend = T,
                           legend = "bottom")


ggsave(fig4, 
       filename = paste0("output/figures/fig4.", figform), 
       width = 12, height = 8, dpi = dpi)


case_prev_oa_yearly %>% filter(strategy == "75+" & scaling_fac == 8)


### Fig. S19 ---------------------------

data_figS19 <- inc %>% 
                filter(strategy==4,
                       seasonal==T,
                       incidence != "total") %>%
  full_join(pop_dat[,c("age", "n_pop")], join_by(age)) %>%
  group_by(age, age_label, year, incidence) %>%
  summarise(n_med = quantile(n_base, 0.5),
            n_q025 = quantile(n_base, 0.025),
            n_q975 = quantile(n_base, 0.975),
            inc_med = quantile(n_base/n_pop*100000, 0.5),
            inc_q025 = quantile(n_base/n_pop*100000, 0.025),
            inc_q975 = quantile(n_base/n_pop*100000, 0.975)) %>%
  ungroup() %>% 
  mutate(incidence = factor(incidence, levels = c("symptomatic", "hospitalised",
                                                  "icu", "mort"),
                            labels = c("sympt", "hosp", "icu", "mort")),
         year=paste0("Year: ", year))

(figS19a <- ggplot(data_figS19) +
    geom_col(aes(x=reorder(age_label, age), y=n_med)) +
    geom_linerange(aes(age, ymin = n_q025, ymax=n_q975)) +
    facet_grid(rows=vars(incidence),
               cols = vars(year), scales = "free_y") +
    ylab("N") +
    labs(title = "Strategy 2: annual cases", tag ="A") +
    xlab(NULL) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_x_discrete(breaks = age_breaks))


(figS19b <- ggplot(data_figS19) +
    geom_col(aes(x=reorder(age_label, age), y=inc_med)) +
    geom_linerange(aes(age, ymin = inc_q025, ymax=inc_q975)) +
    facet_grid(rows=vars(incidence),
               cols = vars(year), scales = "free_y") +
    ylab("N per 100,000") +
    labs(title = "Strategy 2: incidence", tag ="B") +
    xlab(NULL) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_x_discrete(breaks = age_breaks))


figS19 = grid.arrange(figS19a, figS19b, nrow = 2)

ggsave(figS19, 
       filename = paste0("output/figures/figS19.", figform), 
       width = 10, height = 8, dpi = dpi)


### Fig. S20 ----------------------------------------

figS20 <-  plot_incid_prevented(inc, strat=5, strat_name = "Older adult vacc (75+)")
ggsave(figS20, 
       filename = paste0("output/figures/figS20.", figform), 
       width = 10, height = 12.5, dpi = dpi)

### Fig. S21 ----------------------------------------

figS21 <-  plot_incid_prevented(inc, strat=4, strat_name = "Older adult vacc (65+)")
ggsave(figS21, 
       filename = paste0("output/figures/figS21.", figform), 
       width = 10, height = 12.5, dpi = dpi)

### Fig. S22 ----------------------------------------

figS22 <-  plot_incid_prevented(inc, strat=6, strat_name = "Older adult vacc (55+)")
ggsave(figS22, 
       filename = paste0("output/figures/figS22.", figform), 
       width = 10, height = 12.5, dpi = dpi)


# Older adults Incremental strategy 4-6 ----------------------------------------


# Scale-up incidence by various factors
inc_per_strat_oa_up_incr <- do.call(
  rbind,
  lapply(c(1,8,11,14),
         function(x) {
           upscale_oa(inc_incr_extend, x)
         })) %>% 
  group_by(incidence, strategy, seasonal, replicate, year, scaling_fac) %>% 
  dplyr::summarise(n_inc = sum(n),
                   n_base = sum(n_base),
                   n_avert = sum(n_avert)) %>%
  ungroup()


case_prev_oa_incr = inc_per_strat_oa_up_incr %>% 
  filter(seasonal==T) %>%
  group_by(incidence, strategy,
           scaling_fac, replicate) %>%
  mutate(n_avert = sum(n_avert)) %>%
  group_by(incidence, strategy,
           scaling_fac) %>%
  summarise(n_avert_med = median(n_avert),
            n_avert_q025 = quantile(n_avert, 0.025),
            n_avert_q975 = quantile(n_avert, 0.975)) %>%
  mutate(strategy = factor(strategy, levels = rev(c(5,4,6)), 
                           labels = rev(c("75+", "65-74", "55-64")))) %>%
  arrange(incidence, strategy, scaling_fac)

case_prev_oa_incr_yearly = inc_per_strat_oa_up_incr %>% 
  filter(seasonal==T) %>%
  group_by(incidence, strategy,
           scaling_fac, replicate) %>%
  mutate(n_avert = mean(n_avert)) %>%
  group_by(incidence, strategy,
           scaling_fac) %>%
  summarise(n_avert_med = median(n_avert),
            n_avert_q025 = quantile(n_avert, 0.025),
            n_avert_q975 = quantile(n_avert, 0.975)) %>%
  mutate(strategy = factor(strategy, levels = rev(c(5,4,6)), 
                           labels = rev(c("75+", "65-74", "55-64")))) %>%
  arrange(incidence, strategy, scaling_fac)




## Fig. S23 ------------------------------------------------

figS23a = case_prev_oa_incr_yearly %>%
  filter(scaling_fac %in% c(8, 14)) %>%
  mutate(scaling_fac = factor(scaling_fac, levels = c(14, 8),
                              labels = paste0("Scaling: ", c(14, 8)))) %>%
  filter(incidence %in% inc_plot) %>%
  mutate(incidence = factor(incidence, levels = inc_plot)) %>%
  ggplot() + 
  geom_col(aes(strategy, n_avert_med,
               fill = scaling_fac), 
           alpha = .5, col = "black", position = position_dodge2(0.5)) +
  geom_label(aes(strategy, if_else(strategy=="75+"|incidence=="symptomatic", n_avert_q025, n_avert_q975),
                 label = format(round(n_avert_med,0), big.mark = ",", trim = T),
                 hjust = if_else(strategy=="75+"|incidence=="symptomatic", 1, 0),
                 group = scaling_fac),
             fill = "lightgrey", position = position_dodge2(0.9), size=4.5) +
  geom_linerange(aes(strategy, ymin = n_avert_q025, 
                     ymax=n_avert_q975, group = scaling_fac),
                 position = position_dodge2(0.9)) +
  facet_wrap(~incidence, scales = "free_x", ncol = length(inc_plot)) + 
  ylab("Prevented cases") +
  scale_x_discrete(guide = guide_axis(angle = 0)) + xlab("") +
  scale_y_continuous(n.breaks = 3) +
  scale_fill_brewer(type = "qual") +
  coord_flip() +
  theme(legend.title = element_blank())


immu_per_strat_oa_incr = immu_per_strat %>%
  filter(strategy %in% 4:6) %>%
  select(strategy, seasonal, replicate, n_immu) %>%
  mutate(ref_strat = case_when(strategy==5 ~ 2,
                               strategy==4 ~ 5,
                               strategy==6 ~ 4)) %>%
  left_join(immu_per_strat %>%
              filter(strategy %in% c(2,5,4), seasonal = TRUE) %>%
              select(replicate, n_base = n_immu,
                     ref_strat = strategy)) %>%
  mutate(n_immu_add = n_immu-n_base)

nni_oa_incr = inc_per_strat_oa_up_incr %>%
  filter(seasonal==T) %>%
  group_by(incidence, strategy,
           scaling_fac, replicate) %>%
  summarise(n_avert = sum(n_avert)) %>%
  left_join(immu_per_strat_oa_incr %>% select(-n_base)) %>%
  mutate(nni = n_immu_add/n_avert) %>%
  group_by(incidence, strategy, scaling_fac) %>%
  summarise(nni_med = median(nni),
            nni_q025 = quantile(nni, 0.025),
            nni_q975 = quantile(nni, 0.975))  %>%
  mutate(strategy = factor(strategy, levels = rev(c(5,4,6)), 
                           labels = rev(c("75+", "65-74", "55-64"))))



figS23b = nni_oa_incr %>% 
  filter(scaling_fac %in% c(8, 14)) %>%
  mutate(scaling_fac = factor(scaling_fac, levels = c(14, 8),
                              labels = paste0("Scaling: ", c(14, 8)))) %>%
  filter(incidence %in% inc_plot) %>%
  mutate(incidence = factor(incidence, levels = inc_plot)) %>%
  ggplot() + 
  geom_col(aes(strategy, nni_med, fill = scaling_fac), 
           alpha = .5, col = "black", position = position_dodge2(0.5)) +
  geom_label(aes(strategy,
                 if_else(strategy == "55-64" | incidence=="symptomatic", 
                         nni_q025, nni_q975),
                 label = format(round(nni_med,0), big.mark = ","),
                 hjust= if_else(strategy == "55-64" | incidence=="symptomatic", 
                                1, 0),
                 group = scaling_fac),
             fill = "lightgrey", position = position_dodge2(.9), size=4.5) +
  geom_linerange(aes(strategy, ymin = nni_q025, 
                     ymax=nni_q975, group = scaling_fac),
                 position = position_dodge2(0.9)) +
  facet_wrap(~incidence, scales = "free_x", ncol = length(inc_plot)) + 
  ylab("NNV") +
  scale_x_discrete(guide = guide_axis(angle = 0)) + xlab("") +
  scale_fill_brewer(type = "qual") +
  scale_y_continuous(n.breaks = 3) +
  coord_flip() + 
  theme(legend.title = element_blank())


figS23 = ggpubr::ggarrange(figS23a + theme(axis.text = element_text(size=13), 
                                                     axis.title = element_text(size=14),
                                                     strip.text = element_text(size=14)),
                           figS23b + theme(axis.text = element_text(size=13), 
                                               axis.title = element_text(size=14),
                                               strip.text = element_text(size=14)), 
                           nrow = 2,
                           labels = "AUTO", common.legend = T,
                           legend = "bottom")

ggsave(figS23, 
       filename = paste0("output/figures/figS23.", figform), 
       width = 12, height = 8, dpi = dpi)


## Fig. S24 ------------------------------------ 

figS24a = case_prev_oa %>%
  filter(incidence != "total",
         incidence != "symptomatic") %>%
  ggplot() +
  geom_rect(data = tibble(xstart=7.5,
                          xend = 14.5), 
            aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
            fill = "lightgrey",
            alpha = 0.4) +
  geom_point(aes(scaling_fac, n_avert_med, col = strategy),
             position = position_dodge2(width=.5)) +
  geom_linerange(aes(scaling_fac, ymin = n_avert_q025, 
                     ymax=n_avert_q975,col = strategy),
                 position = position_dodge2(width=.5)) +
  facet_wrap(~incidence, ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = c(1,5, 10, 15), minor_breaks = 1:15) +
  
  xlab("Scaling factor") +
  ylab("Prevented incidence") +
  theme(legend.title = element_blank())

nni_oa = inc_per_strat_oa_up %>%
  filter(seasonal==T) %>%
  group_by(incidence, strategy,
           scaling_fac, replicate) %>%
  summarise(n_avert = sum(n_avert)) %>%
  left_join(immu_per_strat_oa %>% select(-n_base)) %>%
  mutate(nni = n_immu_add/n_avert) %>%
  group_by(incidence, strategy, scaling_fac) %>%
  summarise(nni_med = median(nni),
            nni_q025 = quantile(nni, 0.025),
            nni_q975 = quantile(nni, 0.975))  %>%
  mutate(strategy = factor(strategy, levels = c(5,4,6), 
                           labels = c("75+", "65+", "55+")))


figS24b = nni_oa %>%
  filter(incidence != "total",
         incidence != "symptomatic") %>%
  ggplot() +
  geom_rect(data = tibble(xstart=7.5,
                          xend = 14.5), 
            aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
            fill = "lightgrey",
            alpha = 0.4) +
  geom_point(aes(scaling_fac, nni_med , col = strategy),
             position = position_dodge2(width=.5)) +
  geom_linerange(aes(scaling_fac, ymin = nni_q025 , 
                     ymax=nni_q975,col = strategy),
                 position = position_dodge2(width=.5)) +
  facet_wrap(~incidence, ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = c(1,5, 10, 15), minor_breaks = 1:15) +
  xlab("Scaling factor") +
  ylab("NNV") +
  theme(legend.title = element_blank())

figS24 = ggpubr::ggarrange(figS24a,
                           figS24b, nrow = 2, labels = "AUTO", legend = "top", 
                                     common.legend = T)  

ggsave(figS24, 
       filename = paste0("output/figures/figS24.", figform), 
       width = 12, height = 6, dpi = dpi)








