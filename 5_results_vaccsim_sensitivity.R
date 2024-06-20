# Housekeeping  -------------------------------------------------
# Authors: Felix Guenther, Fabienne Krauer
##########################################


library(tidyverse)
theme_set(theme_bw())

source("func_helpers.R")

# fig format
figform <- "png"
dpi <- 600


set.seed(240605122)

# Read and prep data -----------------------------------

inc_1_8 = read_csv('output/vacc_simulations/vaccsim_inc_sensitivity_0_vs_[1, 2, 3, 4, 5, 6, 7, 8].csv')
inc_9_16 = read_csv('output/vacc_simulations/vaccsim_inc_sensitivity_0_vs_[9, 10, 11, 12, 13, 14, 15, 16].csv')
inc_17_18 = read_csv('output/vacc_simulations/vaccsim_inc_sensitivity_0_vs_[17, 18].csv')
inc_19_23 = read_csv('output/vacc_simulations/vaccsim_inc_sensitivity_0_vs_[19, 20, 21, 22, 23].csv')

inc_main = read_csv('output/vacc_simulations/vaccsim_inc_0_vs_[1,2,3].csv') %>% 
          mutate(strategy_label = factor(strategy, 
                                             levels = 1:3,
                                             labels = c("Nirs. HR",	
                                                        "Nirs. all 70%",
                                                        "Mat. vacc. 40% & Pali. HR")),
                 analysis = "main")

inc_sens = inc_1_8 %>% 
  full_join(inc_9_16) %>% 
  full_join(inc_17_18) %>% 
  full_join(inc_19_23) %>% 
  mutate(strategy_label = factor(strategy, 
                                 levels = 1:23,
                                 labels = c("Nirs. all 10%",	
                                      "Nirs. all 20%",	
                                      "Nirs. all 30%",	
                                      "Nirs. all 40%",	
                                      "Nirs. all 50%",	
                                      "Nirs. all 60%",	
                                      "Nirs. all 80%",	
                                      "Nirs. all 90%",	
                                      "Mat. vacc. 10% & Pali. HR",	
                                      "Mat. vacc. 20% & Pali. HR",	
                                      "Mat. vacc. 30% & Pali. HR",	
                                      "Mat. vacc. 50% & Pali. HR",	
                                      "Mat. vacc. 60% & Pali. HR",	
                                      "Mat. vacc. 70% & Pali. HR",	
                                      "Mat. vacc. 80% & Pali. HR",	
                                      "Mat. vacc. 90% & Pali. HR",	
                                      "Mat. vacc. seasonal & Pali. HR",	
                                      "Mat. vacc. year-round & Pali. HR",	
                                      "Nirs. all, Sep-Mar",	
                                      "Nirs. all, Okt-Mar",	
                                      "Nirs. all, Dec-Mar",	
                                      "Nirs. all, Jan-Mar",	
                                      "Nirs. all, Feb-Mar")),
         analysis = "sens")

inc_orig = inc_main %>% full_join(inc_sens)

immu_1_8 = read_csv('output/vacc_simulations/vaccsim_immunisations_sensitivity_0_vs_[1, 2, 3, 4, 5, 6, 7, 8].csv')
immu_9_16 = read_csv('output/vacc_simulations/vaccsim_immunisations_sensitivity_0_vs_[9, 10, 11, 12, 13, 14, 15, 16].csv')
immu_17_18 = read_csv('output/vacc_simulations/vaccsim_immunisations_sensitivity_0_vs_[17, 18].csv')
immu_19_23 = read_csv('output/vacc_simulations/vaccsim_immunisations_sensitivity_0_vs_[19, 20, 21, 22, 23].csv')

immu_main = read_csv('output/vacc_simulations/vaccsim_immunisations_0_vs_[1,2,3].csv') %>% 
  mutate(strategy_label = factor(strategy, levels = 1:3,
                                 labels = c("Nirs. HR",	
                                            "Nirs. all 70%",
                                            "Mat. vacc. 40% & Pali. HR")),
         analysis = "main")

immu_sens = immu_1_8 %>% 
            full_join(immu_9_16) %>% 
  full_join(immu_17_18) %>% 
  full_join(immu_19_23) %>% 
  mutate(strategy_label = factor(strategy, levels = 1:23,
                           labels = c("Nirs. all 10%",	
                                      "Nirs. all 20%",	
                                      "Nirs. all 30%",	
                                      "Nirs. all 40%",	
                                      "Nirs. all 50%",	
                                      "Nirs. all 60%",	
                                      "Nirs. all 80%",	
                                      "Nirs. all 90%",	
                                      "Mat. vacc. 10% & Pali. HR",	
                                      "Mat. vacc. 20% & Pali. HR",	
                                      "Mat. vacc. 30% & Pali. HR",	
                                      "Mat. vacc. 50% & Pali. HR",	
                                      "Mat. vacc. 60% & Pali. HR",	
                                      "Mat. vacc. 70% & Pali. HR",	
                                      "Mat. vacc. 80% & Pali. HR",	
                                      "Mat. vacc. 90% & Pali. HR",	
                                      "Mat. vacc. seasonal & Pali. HR",	
                                      "Mat. vacc. year-round & Pali. HR",	
                                      "Nirs. all, Sep-Mar",	
                                      "Nirs. all, Okt-Mar",	
                                      "Nirs. all, Dec-Mar",	
                                      "Nirs. all, Jan-Mar",	
                                      "Nirs. all, Feb-Mar")),
         analysis = "sens")

immu = immu_main %>% full_join(immu_sens)

# Calculate ICU and mortality incidence ----------------------------------

inek_smry = read_csv2("data/INEK_derived_endpoint_rates.csv")
inc_extend <- calc_icu_deaths(inc_orig, inek_smry)

# Upscale older adults hospitalisations  ----------------------------------------

inc <- inc_extend %>% upscale_oa(scaling_fac = 8)



# Fig. S17 ---------------------

dat_s17 = inc %>% filter(grepl("%", strategy_label),
               seasonal=TRUE,
               incidence %in% c("symptomatic",
                                "hospitalised",
                                "icu",
                                "mort")) %>%
  group_by(incidence, strategy_label, replicate) %>%
  summarise(n=sum(n),
            n_base = sum(n_base),
            n_avert = sum(n_avert)) %>%
  group_by(incidence, strategy_label) %>%
  summarise(post_med = median(n_avert),
            q025 = quantile(n_avert, 0.025),
            q975 = quantile(n_avert, 0.975)) %>%
  mutate(uptake = str_replace(as.character(strategy_label), "Nirs. all ", ""),
         uptake = str_replace(uptake, " & Pali. HR", ""),
         uptake = str_replace(uptake, "Mat. vacc. ", ""),
         uptake = str_replace(uptake, "%", ""),
         uptake = as.numeric(uptake)/100,
         strategy_2 = if_else(grepl("Nirs. all", strategy_label),
                              "Nirs. all",
                              "Mat. vacc."),
         incidence = factor(incidence, 
                            levels = c("symptomatic",
                                       "hospitalised",
                                       "icu",
                                       "mort")))
figS17a = dat_s17 %>%
  ggplot() + 
  geom_point(aes(uptake, post_med, col = strategy_2)) +
  geom_linerange(aes(uptake, ymin = q025, ymax=q975, col = strategy_2),
                 show.legend = FALSE) +
  geom_line(aes(uptake, post_med, col = strategy_2), lty = 2) +
  facet_wrap(~incidence, scales = "free_y", ncol = 4) +
  scale_x_continuous(labels = scales::percent, breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
  scale_color_brewer(type="qual") +
  ylab("Cases prevented\nvs. Palivizumab HR") +
  xlab("Uptake in target group") +
  expand_limits(y=0) +
  theme(legend.title = element_blank())

dat_s17_b = dat_s17 %>%
  left_join(immu %>% group_by(strategy_label, replicate) %>%
              summarise(n=sum(n),
                        n_base = sum(n_base)) %>%
              group_by(strategy_label) %>%
              summarise(add_immu_post_med = median(n-n_base)/1e6)
  )

figS17b = dat_s17_b %>%
  ggplot() +
  geom_segment(aes(x = `add_immu_post_med_Nirs. all`,
                   xend = `add_immu_post_med_Mat. vacc.`,
                   y = `post_med_Nirs. all`,
                   yend = `post_med_Mat. vacc.`),
               data = dat_s17_b %>% select(incidence, 
                                           add_immu_post_med, strategy_2, 
                                           uptake, post_med) %>% 
                 pivot_wider(values_from = c(add_immu_post_med, post_med), 
                             names_from = strategy_2, 
                             id_cols = c(incidence, uptake)),
               col = "black", lty = 2, show.legend = FALSE) +
  geom_point(aes(add_immu_post_med, post_med, col = strategy_2)) +
  geom_linerange(aes(add_immu_post_med, ymin = q025, ymax=q975, col = strategy_2),
                 show.legend = FALSE) +
  geom_line(aes(add_immu_post_med, post_med, col = strategy_2), lty = 2) +
  facet_wrap(~incidence, scales = "free_y", ncol = 4) +
  scale_x_continuous() +
  scale_color_brewer(type="qual") +
  ylab("Cases prevented\nvs. Palivizumab HR") +
  xlab("Number of individuals additionally immunised (millions)") +
  expand_limits(y=0, x=0) +
  theme(legend.title = element_blank())

figS17 = ggpubr::ggarrange(figS17a, figS17b, nrow = 2, common.legend = T, 
                        legend = "bottom",
                        labels = "AUTO")

ggsave(paste0("output/figures/figS17.", figform), 
       figS17, width = 10, height = 6, dpi = dpi)


# Fig. S18 ---------------------

s18_dat = inc %>% 
  filter(strategy_label %in% c("Mat. vacc. seasonal & Pali. HR",	
                         "Mat. vacc. year-round & Pali. HR"),
         incidence %in% c("symptomatic",
                          "hospitalised",
                          "icu",
                          "mort")) %>%
  mutate(incidence = factor(incidence, 
                            levels = c("symptomatic",
                                       "hospitalised",
                                       "icu",
                                       "mort"))) %>%
  left_join(immu %>% 
              mutate(add_immu = n-n_base) %>%
              select(age, year, strategy_label, replicate, add_immu) %>%
              group_by(age, year, strategy_label, replicate) %>%
              summarise(add_immu = sum(add_immu))) %>%
  group_by(incidence, strategy_label, replicate) %>%
  summarise(n=sum(n),
            n_base = sum(n_base),
            n_avert = sum(n_avert),
            add_immu = sum(add_immu)) %>%
  group_by(incidence, strategy_label) %>%
  summarise(post_med = median(n_avert/5),
            q025 = quantile(n_avert/5, 0.025),
            q975 = quantile(n_avert/5, 0.975),
            nnv_med = median(add_immu/n_avert),
            nnv_q025 = quantile(add_immu/n_avert, 0.025),
            nnv_q975 = quantile(add_immu/n_avert, 0.975)) %>%
  mutate(strategy_label = factor(strategy_label,
                           levels = c("Mat. vacc. seasonal & Pali. HR",	
                                      "Mat. vacc. year-round & Pali. HR"),
                           labels = c("Mat vacc. seasonal",
                                      "Mat vacc. year-round")))


figS18a = s18_dat  %>%
  mutate(inc_avert_print = round_print(post_med)) %>%
  ggplot() + 
  geom_col(aes(strategy_label, post_med), 
           fill = "lightgreen", 
           alpha = .5, col = "black") +
  geom_label(aes(strategy_label, 
                 post_med,
                 label = label_form(post_med),
                 hjust=1), fill = "lightgrey",
             vjust=1, nudge_x = .45, size=5
  ) +
  geom_linerange(aes(strategy_label, ymin = q025, 
                     ymax=q975)) +
  facet_wrap(~incidence, scales = "free_x", ncol = 4) + 
  ylab("Prevented cases") +
  scale_x_discrete(guide = guide_axis(angle = 0)) + xlab("") +
  scale_y_continuous(n.breaks = 3) +
  scale_fill_brewer(type = "qual") +
  coord_flip() +
  theme(legend.position = "none")

figS18b = s18_dat  %>%
  mutate(nnv_print = round_print(nnv_med)) %>%
  ggplot() + 
  geom_col(aes(strategy_label, nnv_med  ), 
           fill = "lightgreen", 
           alpha = .5, col = "black") +
  geom_label(aes(strategy_label, 
                 nnv_med  ,
                 label = label_form(nnv_med),
                 hjust=1), fill = "lightgrey",
             vjust=1, nudge_x = .45, size=5
  ) +
  geom_linerange(aes(strategy_label, ymin = nnv_q025, 
                     ymax=nnv_q975)) +
  facet_wrap(~incidence, scales = "free_x", ncol = 4) + 
  ylab("NNV") +
  scale_x_discrete(guide = guide_axis(angle = 0)) + xlab("") +
  scale_y_continuous(n.breaks = 3) +
  scale_fill_brewer(type = "qual") +
  coord_flip() +
  theme(legend.position = "none")

figS18 = ggpubr::ggarrange(figS18a, figS18b, nrow = 2, common.legend = T, 
                        legend = "bottom",
                        labels = "AUTO")

ggsave(paste0("output/figures/figS18.", figform), 
       figS18, width = 10, height = 6, dpi = dpi)

