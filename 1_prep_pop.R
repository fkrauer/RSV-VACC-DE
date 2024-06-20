# Housekeeping  -------------------------------------------------
# Author: Fabienne Krauer
#########################

library(tidyverse)
library(readr)
library(readxl)
library(lubridate)
library(socialmixr)
library(zoo)


#Sys.setlocale("LC_ALL","German")
theme_set(theme_minimal())

# fig format
figform <- "png"
dpi <- 600

# Demographics  -------------------------------------------------

# Function to generate a dataframe with the age groups
make_agegroups <- function(n=25) {
  
  ngroups <- 25
  agegroups <- data.frame(groupno=1:ngroups,
                          unit=c(rep("month", 12), rep("year", 13)),
                          d=c(rep(1,12), rep(1,4), rep(5,2), rep(10,6), 24))
  agegroups$d_y <- ifelse(agegroups$unit=="month", agegroups$d/12, agegroups$d)
  agegroups$age_low_y <- c(seq(0,1,by=1/12)[1:12], seq(1,4, by=1), seq(5,10, by=5), seq(15,65, by=10), 75)
  
  agegroups$age_up_y <- agegroups$age_low_y + agegroups$d_y
  agegroups$label <- ifelse(agegroups$unit=="month", agegroups$age_low_y*12, agegroups$age_low_y)
  agegroups$label <- paste(agegroups$unit, agegroups$label, sep="_")
  agegroups$agegrp <- cut(c(agegroups$age_low_y, 99), 
                              breaks =c(agegroups$age_low_y, 99), 
                              right=FALSE)[1:ngroups]
  agegroups$agemid <- agegroups$age_low_y + (agegroups$age_up_y-agegroups$age_low_y)/2
  
  return(agegroups)
  
}

agegroups_25 <- make_agegroups(25)

# Generates a dataframe with the average births (b0), the amplitude (b1) and the phase shift (b2)
# for a give year or yerrange, assuming a cosine function
make_births <- function(yearrange) {
  
  births <- readr::read_delim("data/destatis/destatis_births_monthly_1950_2022.csv",
                                ";", col_names = TRUE,
                                locale = locale(encoding = 'latin1'),
                                skip=5, trim_ws=TRUE)
    

  # Calculate average daily number of births and amplitude for each year
  colnames(births) <- c("year", "month", "male", "female", "total")
  
  births <- births[grepl("^[0-9]{4}", births$year),] # Remove some tail info
  births$total <- as.numeric(births$total)
  births <- births[!is.na(births$total),]
  births$year <- as.numeric(births$year)
  births <- births %>% filter(year %in% yearrange)
  
  # Sum by year and calculate average daily births
  births$monthno <- match(births$month, format(ISOdate(2004,1:12,1),"%B"))
  births$ym <- as.yearmon(paste0(births$year, "-", births$monthno), format="%Y-%m")
  births$days <- days_in_month(births$ym)
  births$date <- as.Date(births$ym, frac=0.5)
  births <- births %>% group_by(year) %>% 
    dplyr::mutate(doy = julian(date, origin=as.Date(paste0(year[1], "-01-01")))+1)
  births$doy <- ifelse(births$doy==366, 365, births$doy)
  births <- births %>% group_by(year, month) %>% 
    dplyr::mutate(births_d = sum(total, na.rm=T)/sum(days)) %>% ungroup()
  
  # Estimate amplitude and average of daily births assuming a cosine function
  births <- as_tibble(births) %>% group_by(year) %>% 
    do(., fit = nls(births_d ~ b0 * (1 + b1 * cos(2*pi*(doy-b2)/365)), 
                    start=list(b0 = 2000, b1 = 0.1, b2=200), 
                    data=.))
  births <- cbind(year=births[[1]], plyr::ldply(births$fit, coef))
  
  births$groupno <- 1   
  
  return(births)
  
  
}

births_DE_2015_2019 <- make_births(c(2015:2019))


# Generates a dataframe with the population sizes, births and deaths for each age group
# For a given year or yearrange
make_pop_static <- function(agegroups, yearrange, births) {
  
  pop <- readr::read_delim("data/destatis/destatis_popsizes_yearly_age_1970_2020.csv",
                           ";", col_names = TRUE,
                           locale = locale(encoding = 'latin1'),
                           skip=6, trim_ws=TRUE)
  colnames(pop)[1] <- "age_string"
  pop <- pop[1:(nrow(pop)-5),] # Remove some tail info
  pop <- pivot_longer(pop, cols=2:ncol(pop), names_to = "date", values_to = "n_pop")
  pop$date <- as.Date(pop$date)
  pop <- pop %>% filter(age_string!="Total") %>% 
    dplyr::mutate(year = as.numeric(format(date, format="%Y"))) %>% 
    filter(year %in% yearrange) # filter year of interest

  # split under 1-years into 12 bins (monthly age groups)
  foo <- pop[pop$age_string=="under 1 year",]
  foo$n_pop <- foo$n_pop/12
  foo <- foo[rep(seq_len(nrow(foo)), each = 12),]
  foo$age_low_y <- rep(seq(0,1,by=1/12)[1:12], length(yearrange))

  # Define lower age limits
  pop <- pop %>% filter(age_string!="under 1 year") %>%
                  dplyr::mutate(age_low_y = as.numeric(gsub("(^[0-9]{1,2})\\s.+", "\\1", age_string)))
  pop <- merge(pop, foo, by=intersect(names(pop), names(foo)), all=T)
  
  # Bin according to desired population age groups
  pop$agegrp = cut(c(pop$age_low_y, 99),
                   breaks = c(agegroups$age_low_y, 99), right = FALSE)[1:nrow(pop)]
  
  # Group according to age groups and year and sum population
  pop <- pop %>% dplyr::group_by(agegrp, year) %>% 
                  dplyr::summarise(n = sum(n_pop, na.rm=T))
  
  pop <- pop %>% dplyr::group_by(year) %>% 
                  dplyr::mutate(prop = n/sum(n)) %>% arrange(agegrp)
                
  pop <- merge(pop, agegroups, by="agegrp") %>% arrange(groupno)


  # Calculate deaths based on transitions from lower age group and births
  deaths <- pop %>% 
        select(year, groupno, n, d_y) %>% 
        merge(., births, by=c("year", "groupno"), all=T) %>% 
        group_by(year) %>% 
        arrange(groupno) %>% 
        dplyr::mutate(b0 = ifelse(is.na(b0), 0, b0),
                      n_out = n/(365*d_y),
                      n_in = ifelse(b0==0, lag(n_out,1), b0),
                      deaths_d = n_in - n_out) %>% 
        select(-n_in, -n_out)
  
  # Combine all data
  out <- merge(pop, deaths, by=intersect(names(pop), names(deaths)), all=T)
  out <- out[order(out$year, out$groupno),]
  out$b1 <- ifelse(is.na(out$b1),0,out$b1)
  out$b2 <- ifelse(is.na(out$b2),0,out$b2)
  out <- out %>% arrange(year, groupno)

  return(out)
  
}


pop_DE_25_2015_2019 <- make_pop_static(agegroups_25, 
                                       c(2015:2019), 
                                       births_DE_2015_2019)


#write.csv(pop_DE_25_2015_2019, file=paste0("data/pop_DE_25_2015_2019.csv"), 
#            row.names = FALSE)
saveRDS(pop_DE_25_2015_2019, "data/pop_DE_25_2015_2019.rds")


# Contact matrix data -------------------------------------------------

make_cmat <- function(pop, popyear, countries, correct=FALSE, reciprocal=FALSE) {
  
  
  pop <- pop %>% filter(year==popyear) %>% 
        arrange(groupno) 
  popn <- pop$n
  labels <- pop$agegrp
  ngroups <- length(popn)
  
  library(socialmixr)
  data(polymod)

    # from row to col:
    contacts_raw <- contact_matrix(polymod, 
                                 countries = countries, 
                                 age.limits = c(0,1,2,3,4,5,10,15,25,35,45,55,65,75), 
                                 symmetric = FALSE,
                                 counts = FALSE,
                                 weigh.dayofweek = TRUE,
                                 split = FALSE)$matrix
  
  
  if (correct) {
    # replace 0.0 contacts with the minimum
    contacts_raw[contacts_raw==0.0] <- min(contacts_raw[contacts_raw!=0.0])
  }


  # Expand contact rows for age groups <1 year old (no specific data available): contacts / 12
  contacts_raw <- cbind(matrix(rep(contacts_raw[,1] / 12, 12), nrow = 14, ncol = 12), 
                        contacts_raw[,2:14])
  
  contacts_raw <- rbind(matrix(rep(contacts_raw[1,] / 12, 12), nrow = 12, ncol = 25, 
                               byrow=T), 
                        contacts_raw[2:14,])
  
  colnames(contacts_raw) <- rownames(contacts_raw) <- labels
  

  if (reciprocal) {
    # Rescale (make reciprocal) according to population numbers in the model
    out <- matrix(rep(NA, ngroups^2), ngroups, ngroups, byrow=T)
    
    for (i in 1:ngroups) {
      for (j in 1:ngroups) {
        out[i,j] <- (contacts_raw[i,j] * popn[i] + contacts_raw[j,i] * popn[j])/(2*popn[i])
      }
    }
    
    colnames(out) <- rownames(out) <- labels


    return(out)
    
  } else {
    
    return(contacts_raw)
    
  }

}
      

contacts_25_raw <- make_cmat(pop_DE_25_2015_2019, 2019, c("Germany"), correct=TRUE, reciprocal=FALSE)

write.table(contacts_25_raw, file=paste0("data/contacts_25_raw.csv"), 
            row.names = FALSE, col.names = FALSE, sep=", ")



make_cdf <- function(cmat, labels) {
  
  cdf <- data.frame(cmat)
  colnames(cdf) <- labels #rownames(cdf)
  cdf$from <- labels #rownames(cdf)
  cdf$order_from <- 1:nrow(cdf)
  cdf <- pivot_longer(cdf, cols=1:nrow(cdf), 
                      names_to = "to", 
                      values_to = "n_contacts")
  
  levels <- unique(cdf[,c("from", "order_from")])
  levels <- levels[order(levels$order_from),]
  cdf$from <- factor(cdf$from, levels=levels$from)
  cdf$to <- factor(cdf$to, levels=levels$from)
  
  return(cdf)

}


labels <- pop_DE_25_2015_2019 %>% 
          dplyr::filter(year==2019) %>% 
          dplyr::mutate(labels = ifelse(groupno<=12, label, as.character(agegrp))) %>% 
          select(labels) %>% 
          dplyr::pull(labels)

labels <- gsub("_", " ", labels)



cdf_25 <- make_cdf(contacts_25_raw, labels)


# Fig. S1 -----------------------------------

(figS1 <- ggplot(cdf_25) +
        geom_tile(aes(x=from, 
                      y=to, 
                      fill=n_contacts)) + 
        scale_fill_viridis_c(option="magma", direction=-1) +
        geom_tile(data=cdf_25[cdf_25$n_contacts==0.0,], aes(x=from, 
                                                            y=to), fill="grey") +
        labs(fill="N contacts") +
    theme(axis.text.x = element_text(angle=90, hjust=1)))


ggsave(paste0("output/figures/figS1.", figform), 
       plot=figS1, 
       width = 16, height = 14, units="cm", dpi = dpi)

