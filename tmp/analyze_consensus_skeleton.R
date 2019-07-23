# Infering
library(miic)

data(hematoData)
miic.res = miic(inputData = hematoData, latent = TRUE, confidenceShuffle = 10,
                confidenceThreshold = 0.001, doConsensus = 80, nSkeletons = 10)

# Analyzing
library(dplyr)
as_tibble(skeletons) %>% group_by(x, y) %>% mutate(count = n()*100/10) %>% View

as_tibble(skeletons) %>% group_by(x, y) %>% summarise(i = max(I)) %>% View

as_tibble(skeletons) %>% group_by(x, y) %>% summarise(i = mean(as.numeric(I)))

library(ggplot2)
as_tibble(skeletons) %>%
  group_by(x, y) %>%
  summarise(count = n(), I = mean(as.numeric(I))) %>%
  ggplot(aes(x=count, y=I)) + geom_point() + geom_jitter(height = 5) +
  geom_smooth()

library(ggplot2)
as_tibble(skeletons) %>%
  group_by(x, y) %>%
  summarise(count = n(),
            ai_vect = mean(as.numeric(ai_vect_n)),
            I = mean(as.numeric(I))) %>%
  ggplot(aes(x=count, y=ai_vect, color = log(I))) + geom_jitter(height = 0.01) + geom_point() +
  geom_smooth()
