library(dplyr)

## Summary stats table for all your vars
my_vars <- your_data %>% select(., covar1:covarN) # select all the metric variables you want to include

my_summary <- my_vars %>% summarise_each(funs(
  n=sum(!is.na(.)),
  Mean = mean(.,na.rm = TRUE),
  Median = median(.,na.rm = TRUE),
  SD = sd(.,na.rm = TRUE),
  Min = min(.,na.rm = TRUE),
  Q25 = quantile(., 0.25, na.rm = T),
  Q75 = quantile(., 0.75, na.rm = T),
  Max = max(.,na.rm = TRUE)
))

my_table <- my_summary %>% tidyr::gather(stat, val) %>%
  tidyr::separate(stat, into = c("var", "stat"), sep = "_") %>%
  tidyr::spread(stat, val)


## Count observations
ggplot(your_data) +
  geom_bar(aes(factor(manatee)))


## Count observations separating by another covariate
ggplot(your_data) +
  geom_bar(aes(factor(manatee))) +
  facet_wrap(~covariate)


## Violin/Boxplots (distribution of covariates by manatee/no manatee)
ggplot(EMA, aes(factor(manatee), covariate)) +
  geom_violin() +
  geom_boxplot() +
  geom_point(alpha=.4)


## Multiple covariates (with similar scale) in the same plot
long_data <- your_data %>% gather(var,value, covar1:covarN) # gather covariates into long format

ggplot(long_data, aes(var, value, fill=factor(manatee))) +
  geom_violin() +
  geom_boxplot() +
  geom_point(alpha=.4)

