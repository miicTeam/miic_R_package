# Example of command line using bootstraping
library(miic)

data(hematoData)
miic.res = miic(inputData = hematoData, latent = TRUE, confidenceShuffle = 10,
                confidenceThreshold = 0.001, doConsensus = 80, nSkeletons = 10)

# Some plots to evaluate the robustness of MIIC
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

# Consensus
# The code below does the same thing that the base R snippet within miic to filte out edges
# less common than the consensus
library(dplyr)
as_tibble(skeletons) %>%
  group_by(x, y) %>%
  mutate(count = n()) %>%
  filter(count >= 8) %>%
  select(x,y) %>%
  distinct(x,y) %>%
  View

# Comparing final network to consensus skeleton
final_network <- as.data.frame(cbind(miic.res$retained.edges.summary$x,
                                     miic.res$retained.edges.summary$y))
colnames(final_network) <- colnames(consensus_table) <- c('x', 'y')

# Proportion of final network which is not in the consensus table
n_edges_not_in_consensus <- nrow(dplyr::setdiff(final_network, consensus_table))
n_total_final_network <- nrow(final_network)
prop_edges_not_in_consensus <- n_edges_not_in_consensus*100/n_total_final_network
prop_edges_in_consensus <- 100-prop_edges_not_in_consensus
