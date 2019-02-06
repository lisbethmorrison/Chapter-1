library(dplyr)
library(ggplot2)

ds <- data.frame(y = runif(10), x = runif(10)) ## 10 random values 2 columns (y and x)
ds

# Using spearman correlation to be consistent with the next example
cor(ds$y, ds$x, method = "spearman")
## -0.2

# Resampled the y only (reshuffled the order)
ds$resample_y <- sample(ds$y)
ds ## y is reshuffled (without replacement)

# Correlation again, with the resampled y
cor(ds$resample_y, ds$x, method = "spearman")
## 0.054

# 1000 resamples of the correlation of each.
data_frame(num = 1:1000) %>% 
  group_by(num) %>% 
  mutate(corr = cor(sample(ds$y), ds$x, method = "spearman")) %>% 
  ggplot(aes(x = corr)) +
  geom_freqpoly()
