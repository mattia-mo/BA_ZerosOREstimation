# load required packages
library(tidyverse)
library(ggplot2)
if (!requireNamespace("naniar")) install.packages("naniar")
library(naniar)
library(ggpubr)

# read data
estimates <- readRDS("Data_simulated/1estimates_big.rds")

est_CI <- estimates %>% 
  filter(or >=1) %>% 
  mutate_at(vars(starts_with("upper") | starts_with("lower")),
            ~ ifelse(is.na(.) | is.infinite(.), NA, .)) %>% 
  pivot_longer(cols = starts_with("or_"), names_to = "method", values_to = "estimation_or") %>%
  mutate(method = as.factor(str_extract(method, "[^_]+$"))) %>% 
  group_by(rep,scenario) %>% 
  #das geht noch schöner
  mutate(lower = ifelse(method == "wald", lower_wald, ifelse(method == "fisher", lower_fisher, 
                        ifelse(method == "small", lower_small, ifelse(method == "midp", lower_midp, NA)))),
         upper = ifelse(method == "wald", upper_wald, ifelse(method == "fisher", upper_fisher, 
                        ifelse(method == "small", upper_small, ifelse(method == "midp", upper_midp, NA))))) %>% 
  select(!(starts_with("upper_") | starts_with("lower_"))) %>% 
  mutate(across(starts_with("estimation_"), ~na_if(., Inf)),
         log_or = round(log(or), 2), log_est = round(log(estimation_or), 2),
         log_lower = round(log(lower), 2), log_upper = round(log(upper), 2),
         across(starts_with("log_"), ~na_if(., Inf))) 

# zero handling
# complete case analysis
est_CCA <- est_CI %>% 
  group_by(table) %>% 
  mutate(estimation_or = if(any(is.na(estimation_or))) NA else estimation_or,
         lower = if(any(is.na(estimation_or))) NA else lower,
         upper = if(any(is.na(estimation_or))) NA else upper,
         log_est = if(any(is.na(estimation_or))) NA else log_est,
         log_lower = if(any(is.na(estimation_or))) NA else log_lower,
         log_upper = if(any(is.na(estimation_or))) NA else log_upper)

est_CCA$handling <- "CCA"

# Available case analysis
est_ACA <- est_CI
est_ACA$handling <- "ACA"
# unavailable cases are already NA

est_handled <- rbind(est_CCA, est_ACA)

#getting an overview about zeros in tables
zero_perc <- est_handled %>% 
  group_by(scenario, rep) %>% 
  summarize(zero_in_table = any(zero), .groups = "drop") %>% 
  group_by(scenario) %>% 
  summarize(zero_perc = sum(zero_in_table)/max(rep))

scenarios_overview <- estimates %>%
  select(scenario, or, n_obs, prop) %>%
  unique() %>% 
  arrange(scenario) %>% 
  as.data.frame() %>% 
  merge(zero_perc, by = "scenario", all.x = T)

scen_zeros <- zero_perc %>% filter(zero_perc !=0) %>% 
  select(scenario) %>% 
  unlist() %>% 
  unname()
scen_zeros_ordered <- zero_perc %>% filter(zero_perc !=0) %>% 
  arrange(zero_perc) %>% 
  select(scenario, zero_perc)

# visualization of estimates #####################################################
# EDA
EDA_Plots <- est_CI %>% 
  group_by(scenario) %>% 
  group_split() %>% 
  map( ~ ggplot(.x) +
         geom_histogram(aes(estimation_or)) +
         geom_vline(aes(xintercept = or)) +
         facet_wrap(~method) + ggtitle(unique(.x$scenario)))
# all estimates
swarmplot <- est_handled %>%
  filter(handling == "ACA") %>% 
  mutate(n_and_prop = interaction(n_obs, prop)) %>% 
  ggplot() +
  geom_point(aes(x=rep,y=log_est, color = method), alpha = 0.3) +
  geom_hline(aes(yintercept = log_or), color = "red") +
  facet_grid(log_or ~ n_and_prop, scales = "fixed") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.text = element_text(), legend.title = element_blank())+
  labs(y="")
swarmplot

# estimates of zero tables emphasized
est_handled %>% 
  filter(handling == "ACA") %>% 
  mutate(zero_and_performed = ifelse(zero, paste0("TRUE_", method), "FALSE")) %>% 
  ggplot(aes(rep, log_est, color = zero_and_performed)) +
  geom_point(data = . %>% filter(!startsWith(zero_and_performed, "TRUE_")), alpha = 0.7, position = "identity") +
  geom_point(data = . %>% filter(startsWith(zero_and_performed, "TRUE_")), alpha = 0.5) +
  geom_hline(aes(yintercept = log_or)) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE_fisher" = "green", "TRUE_small" = "red", "TRUE_wald" = "blue", "TRUE_midp" = "orange")) +
  facet_grid(or + log_or~n_obs + prop, labeller = partial(label_both, sep = " = ")) +
  labs(title = "Log-Estimation of OR with performances despite zero in table emphasized",
       x = "Tables", y = expression(log(hat(theta))),
       color = "Table containing Zero & \n method that performed") +
  theme(legend.title = element_text(size = 10),
        axis.ticks.x = element_blank(), axis.text.x = element_blank())

# NAs plotted
est_handled %>%
  filter(zero) %>%
  ggplot() +
  geom_miss_point(aes(rep, estimation_or, color = method)) +
  geom_hline(aes(yintercept = or)) +
  facet_wrap(scenario ~ handling, ncol = 4, labeller = partial(label_both, sep = " = ", 
                                                               multi_line = T)) +
  labs(title = "tables containing zeros: difference of handling them - missing values included",
       y = "estimated odds ratio", x = "") +
  theme(strip.text = element_text(
    size = 8, 
    lineheight = 0.6),
  axis.ticks.x = element_blank(), axis.text.x = element_blank())

# performance measurements #######################################################

# n after handling and NA removed
est_handled %>% 
  filter(!is.na(estimation_or)) %>% 
  group_by(scenario, method, handling) %>%
  mutate(n_after_handling = n()) %>%
  ungroup() %>% 
  mutate(size = if_else(handling == "ACA", 2, 1),
         size = as.factor(size)) %>% 
  ggplot() +
  geom_line(aes(scenario, n_after_handling, color = handling, linewidth = size)) +
  facet_wrap(~method, ncol = 2) +
  scale_linewidth_manual(values = c(0.7, 0.9))

# coverage
est_handled %>% 
  mutate(inside_CI = lower<or & upper>or,
         inside_CI_log = log_lower<log_or & log_upper>log_or) %>% 
  group_by(scenario, method, handling) %>%
  summarize(n_performed = sum(!is.na(estimation_or)),
            sumxx = sum(inside_CI, na.rm = TRUE),
            coverage = sum(inside_CI, na.rm = TRUE)/n_performed,
            coverage_log = sum(inside_CI_log, na.rm = TRUE)/n_performed) %>% 
  merge(zero_perc, by = "scenario", all.x = T) %>% 
  ggplot() +
  geom_point(aes(zero_perc, coverage, color = handling)) +
  geom_line(aes(zero_perc, coverage, color = handling)) +
  facet_wrap(~method, ncol = 2)

# closest estimation
find_closest_estimation <- function(x, log_or, method) {
  a <- x - log_or
  if(all(is.na(a))){
    closest_method <- NA
  } else {
    closest_index <- which.min(abs(a))
    closest_method <- as.character(method[closest_index])
  }
  return(closest_method)
}
est_handled %>% 
  filter(scenario %in% scen_zeros) %>% 
  group_by(scenario, rep, zero, handling) %>% 
  summarize(best = find_closest_estimation(log_est, log_or, method), .groups = "drop") %>% 
  filter(!is.na(best)) %>%
  mutate(best=factor(best, levels = c("fisher", "midp", "wald", "small"))) %>% 
  merge(scen_zeros_ordered, by = "scenario", all.x = T) %>% 
  mutate(scenario = as.character(scenario),
         scenario = factor(scenario, levels = as.character(scen_zeros_ordered$scenario), 
                           labels = paste(scen_zeros_ordered$scenario, scen_zeros_ordered$zero_perc, sep = ": "))) %>%
  ggplot() +
  geom_bar(aes(x =  scenario, fill = best), position = "fill") +
  facet_wrap(~handling) +
  labs(title = expression(paste("Percentage of clostest log(", hat(theta), ") of the methods")), x = "scenario", y = "") +
  theme(axis.text.x = element_text(size = 6.5, angle = 45))

# bias
est_handled %>% 
  mutate(ratio_ors = log(estimation_or/or),
         across(starts_with("ratio"), ~replace(., is.infinite(.), NA))) %>% 
  group_by(scenario, method, handling) %>%
  summarize(n_performed = sum(!is.na(estimation_or)),
            me = sum(ratio_ors, na.rm = TRUE)/n_performed, .groups = "drop") %>% 
  merge(zero_perc, by = "scenario", all.x = T) %>% 
  ggplot() +
  geom_point(aes(zero_perc, me, color = handling)) +
  facet_wrap(~method, ncol = 2) +
  labs(y = expression(paste("mean of  ", log(frac(hat(theta),theta)))))
ggsave("plots/PM_bias.png", width = 10, height = 7) 

# MSE
est_handled %>% 
  mutate(dist_squared = (estimation_or - or)^2,
         dist_squared_log = (log_est - log_or)^2) %>% 
  group_by(scenario, method, handling) %>%
  summarize(n_performed = sum(!is.na(estimation_or)),
            mse = sum(dist_squared, na.rm = TRUE)/n_performed,
            mse_log = sum(dist_squared_log, na.rm = TRUE)/n_performed, .groups = "drop") %>% 
  merge(zero_perc, by = "scenario", all.x = T) %>% 
  ggplot() +
  geom_line(aes(zero_perc, mse_log, color = handling), size = 0.7) +
  geom_point(aes(zero_perc, mse_log, color = handling), size = 1.4) +
  facet_wrap(~method, ncol = 2) +
  labs(y = expression(paste("mean of MSE with log(", hat(theta), ")")))
