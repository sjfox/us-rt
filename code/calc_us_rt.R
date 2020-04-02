library(incidence)
library(EpiEstim)
library(ggplot2)
library(pracma)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

nyc_fips <- "000001"
kc_fips <- "000002"
county_data <- read_csv('data/us-counties.csv') %>% 
  mutate(fips = ifelse(county == "New York City", nyc_fips, fips),
         fips = ifelse(county == "Kansas City", kc_fips, fips)) 

counties_to_estimate <- county_data %>% 
  arrange(desc(date)) %>% 
  group_by(county) %>% 
  slice(1) %>% 
  ungroup() %>% 
  filter(cases > 20, !is.na(fips)) %>% pull(fips)

county_incidence <- county_data %>% 
  arrange(fips, date) %>% 
  group_by(fips) %>% 
  mutate(new_cases = diff(c(0,cases))) %>%
  mutate(new_cases = ifelse(new_cases<0, 0, new_cases)) %>% 
  filter(fips %in% counties_to_estimate)


get_rt <- function(df, county){
  ## Convert cases into an incidence object (f)
  cty_incidence <- incidence(dates = rep(df$date, times = df$new_cases))
  
  ## Disease parameters
  shape=3.63
  rate = 0.71
  mean = shape/rate
  std = (shape/rate^2)^(0.5)
  
  ## Epidemic timing information
  if(dim(cty_incidence)[1] < 9) return(NULL)
  t_start <- seq(2, dim(cty_incidence)[1] - (7-1))   
  t_end <- t_start + 7-1
  
  
  ## Run estimation
  res_with_imports <- try(estimate_R(cty_incidence, method = "parametric_si", config = make_config(list( mean_si = mean, std_si = std, t_start = t_start, 
                                                                                                         t_end = t_end))))
  lengR = dim(res_with_imports$R)[1]
  lengT = length(res_with_imports$dates) 
  
  return(tibble(date = res_with_imports$dates[(lengT-lengR+1):lengT],
                median_rt = res_with_imports$R$`Median(R)`,
                lb_rt = res_with_imports$R$`Quantile.0.025(R)`,
                ub_rt = res_with_imports$R$`Quantile.0.975(R)`))
  # return(res_with_imports)
}

county_rt <- county_incidence %>% 
  select(-deaths, -cases) %>% 
  nest(date, new_cases) %>% 
  mutate(rt = map2(data, fips, get_rt))


plot_eff_r0 <- function(df, pl_title) {
  df %>% 
    ggplot(aes(date, median_rt)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = lb_rt, ymax = ub_rt), alpha = .2) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_hline(yintercept = 0, lty = 2, alpha =0) +
    background_grid(major = "yx", minor = 'y') +
    labs(x = "", y = "Effective R0", title = pl_title) 
}

save_plots_pdf = function(list, filename) {
  #start pdf
  pdf(filename)
  
  #loop
  for (p in list) {
    print(p)
  }
  
  #end pdf
  dev.off()
  
  invisible(NULL)
}

plot_df <- county_rt %>% 
  filter(! (map(rt, is.null) %>% unlist())) %>% 
  mutate(plot = map2(rt,  paste0(county, ", ", state), plot_eff_r0)) 

if(!dir.exists("figs")){
  dir.create("figs")
}
save_plots_pdf(plot_df %>% pull(plot), 'figs/us-city-rt.pdf')
