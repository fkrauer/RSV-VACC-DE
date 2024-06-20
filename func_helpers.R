# Helper functions
# Author: Felix Guenther

# labelling
label_form = function(x) {
  sapply(x,
         function(y) {
           format(y, big.mark = ",", trim = TRUE, scientific = FALSE, 
                  digits=ifelse(y>1, 2, 1), nsmall=ifelse(y<10, 1, 0))
         }
  )
}

# Rounding
round_print = function(x) {
  case_when(x>1000000 ~ round(x, -4),
            x>100000 ~ round(x, -3),
            x>1000 ~ round(x, -2),
            x>100 ~ round(x, -1),
            x>10 ~ round(x, -0),
            x<=10 ~ round(x, 1))
}

# Age breaks
age_breaks = c(paste0(0:11, "m"),
               paste0(1:4, "y"),
               "5-9y",
               "10-14y",
               "15-24y",
               "25-34y",
               "35-44y",
               "45-54y",
               "55-64y",
               "65-74y",
               "75+")[seq(1,25, 2)]

label_inc <- function(inc) {
  
  out <- inc %>% mutate(strategy_label = factor(strategy, 
                                 levels = 1:6,
                                 labels = c("Nirs. HR",	
                                            "Nirs. all 70%",
                                            "Mat. vacc. 40% & Pali. HR",
                                            "OA vacc. 65+ & Nirs. all",	
                                            "OA vacc. 75+ & Nirs. all",	
                                            "OA vacc. 55+ & Nirs. all")),
         age_label = factor(age, 
                            levels = 1:25, 
                            labels = c(paste0(0:11, "m"),
                                       paste0(1:4, "y"),
                                       "5-9y",
                                       "10-14y",
                                       "15-24y",
                                       "25-34y",
                                       "35-44y",
                                       "45-54y",
                                       "55-64y",
                                       "65-74y",
                                       "75+")))
  
  return(out)
  
}


# Calculation of yearly percentage of averted cases
perc_avert_yearly = function(age_groups = 1:6,
                             sel_incidence = c("symptomatic", "hospitalised", "icu")) {
  inc %>% 
    filter(age %in% age_groups,
           strategy %in% 1:3,
           incidence %in% sel_incidence,
           seasonal==TRUE) %>%
    select(-seasonal) %>%
    group_by(incidence, strategy, replicate, year) %>% 
    summarise(n_inc = sum(n),
              n_base = sum(n_base),
              n_avert = sum(n_avert)) %>%
    group_by(incidence, strategy,  replicate) %>%
    summarise(mean_inc = mean(n_inc),
              mean_base = mean(n_base),
              mean_avert = mean(n_avert)) %>%
    ungroup() %>%
    group_by(incidence, strategy) %>%
    summarise(av_perc_avert_yearly_med = round(median(mean_avert/mean_base)*100,1),
              av_perc_avert_yearly_q025 = round(quantile(mean_avert/mean_base, 0.025)*100,1),
              av_perc_avert_yearly_q975 = round(quantile(mean_avert/mean_base, 0.975)*100,1),
              mean_avert_med = round(median(mean_avert)),
              mean_avert_q025 = round(quantile(mean_avert, 0.025)),
              mean_avert_q975 = round(quantile(mean_avert, 0.975))
    ) %>%
    mutate(strategy = factor(strategy, levels = 3:1, labels = rev(str_replace(inf_labels, "\n", " ")))) %>%
    mutate(rel_redu = paste0(av_perc_avert_yearly_med, "% (",av_perc_avert_yearly_q025,"-",av_perc_avert_yearly_q975, ")"),
           abs_redu = paste0(mean_avert_med, " (", mean_avert_q025,"-",mean_avert_q975, ")")
    ) %>%
    select(incidence, strategy, abs_redu, rel_redu) %>%
    pivot_longer(cols = c(abs_redu, rel_redu)) %>%
    pivot_wider(names_from=incidence, values_from = c(value)) %>%
    mutate(age_cat = paste0(age_groups, collapse = "")) %>%
    select(all_of(c("age_cat", "strategy", type = "name", sel_incidence)))
}


# scaling up of older adult hospitalisations
upscale_oa = function(inc, scaling_fac = 1) {
  inc_up = inc %>% 
    mutate(n_base = if_else(incidence %in% c("hospitalised", "icu", "mort") & age >= 20,
                         n_base*scaling_fac,
                         n_base),
           n = if_else(incidence %in% c("hospitalised", "icu", "mort") & age >= 20,
                       n*scaling_fac,
                       n),
           n_avert = n_base - n) %>% 
    mutate(scaling_fac = scaling_fac)
  
  return(inc_up)
  
}



# Plotting cases and incidence prevented
plot_incid_prevented = function(inc, strat=1, strat_name = "Nirs. HR") {
  
  
  df = inc %>% filter(strategy %in% strat,
                      seasonal==T,
                      incidence != "total") %>%
    full_join(pop_dat, join_by(age)) %>%
    group_by(age, age_label, year, incidence, strategy, strategy_label) %>%
    summarise(n_avert_med = quantile(n_avert, 0.5),
              n_avert_q025 = quantile(n_avert, 0.025),
              n_avert_q975 = quantile(n_avert, 0.975),
              inc_avert_med = quantile(n_avert/n_pop*100000, 0.5),
              inc_avert_q025 = quantile(n_avert/n_pop*100000, 0.025),
              inc_avert_q975 = quantile(n_avert/n_pop*100000, 0.975),
              prop_avert_med = quantile(prop_avert, 0.5),
              prop_avert_q025 = quantile(prop_avert, 0.025),
              prop_avert_q975 = quantile(prop_avert, 0.975)) %>%
    ungroup() %>% 
    mutate(incidence = factor(incidence, levels = c("symptomatic", "hospitalised",
                                                    "icu", "mort"),
                              labels = c("sympt", "hosp", "icu", "mort")),
           year=paste0("Year: ", year)) 
  
  
  
    figa <- ggplot(df) +
    geom_col(aes(reorder(age_label, age), n_avert_med)) +
    geom_linerange(aes(reorder(age_label, age), ymin = n_avert_q025, ymax=n_avert_q975)) +
    facet_grid(rows=vars(incidence),
               cols = vars(year), scales = "free_y") +
    ylab("Cases prevented, absolute") +
    labs(title = strat_name, tag = "A") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

 
    figb <- ggplot(df) +
    ggtitle(strat_name) +
    geom_col(aes(reorder(age_label, age), inc_avert_med)) +
    geom_linerange(aes(reorder(age_label, age), ymin = inc_avert_q025, ymax=inc_avert_q975)) +
    facet_grid(rows=vars(incidence),
               cols = vars(year), scales = "free_y") +
    ylab("Incidence prevented, rate per 100k") +
    labs(title = strat_name, tag = "B") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
    figc <- ggplot(df) +
    ggtitle(strat_name) +
    geom_col(aes(reorder(age_label, age), prop_avert_med)) +
    geom_linerange(aes(reorder(age_label, age), ymin = prop_avert_q025, ymax=prop_avert_q975)) +
    facet_grid(rows=vars(incidence),
               cols = vars(year), scales = "free_y") +
    ylab("Percentage prevented compared to base strategy") + 
    labs(title = strat_name, tag = "C") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_y_continuous(labels = scales::percent)
  
  fig = grid.arrange(figa, figb, figc, nrow = 3)  
  
  return(fig)
  
}

calc_icu_deaths <- function(inc, inek_smry) {
  
  icu_mort_rate = expand_grid(replicate = 1:100,
                              inek_smry) %>%
    mutate(hosp_icu_ratio = exp(rnorm(n(), pred_icu_link, pred_icu_link_se)),
           hosp_mort_ratio = exp(rnorm(n(), pred_mort_link, pred_mort_link_se))) %>%
    select(replicate, age, hosp_icu_ratio, hosp_mort_ratio)
  
  inc_icu = inc %>% filter(incidence=="hospitalised") %>%
    left_join(icu_mort_rate, by = c("age", "replicate")) %>%
    mutate(n=n*hosp_icu_ratio,
           n_base = n_base*hosp_icu_ratio,
           n_avert = n_base - n,
           prop_avert = n_avert/n_base,
           incidence = "icu") %>%
    select(-hosp_icu_ratio, -hosp_mort_ratio)
  
  inc_death = inc %>% filter(incidence=="hospitalised") %>%
    left_join(icu_mort_rate, by = c("age", "replicate")) %>%
    mutate(n=n*hosp_mort_ratio,
           n_base = n_base*hosp_mort_ratio,
           n_avert = n_base - n,
           prop_avert = n_avert/n_base,
           incidence = "mort") %>%
    select(-hosp_icu_ratio, -hosp_mort_ratio)
  
  
  # Combine incidence of different endpoints
  out = inc %>% 
    full_join(inc_icu) %>%
    full_join(inc_death)
  
  return(out)
  
  
}

