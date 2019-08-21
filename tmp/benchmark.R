library(dplyr)
library(ggplot2)
library(optparse)
library(miic)
library(bnlearn)

# Running benchmark evaluation algorithm per se
benchmark <- function(nSkeletons, consensus, fraction, seed=2019,
                      edgeFiltering=0, dataset, nSamples, outputPath) {
  # Which dataset?
  if (dataset == 'insurance') {
    setwd('/home/mribeirodantas/dev/miic_r_package/tmp/INSURANCE/')
    myBnNet = read.bif(file='insurance.bif', debug = FALSE)
    input = try(rbn(myBnNet, n=nSamples, debug = F))
  }

  # Running MIIC
  setwd(outputPath)
  miic.res = miic(inputData = input, latent = TRUE,
                  confidenceShuffle = 10, confidenceThreshold = 0.001,
                  doConsensus = consensus,
                  nSkeletons = nSkeletons,
                  proportionToUndersample = fraction,
                  whereToEdgeFilter=edgeFiltering)
  # Plotting average CI by edge
  as_tibble(miic.res$skeletons) %>%
    group_by(x, y) %>%
    summarise(count = n(), I = mean(as.numeric(I))) %>%
    ggplot(aes(x=count*100/nSkeletons, y=I)) + geom_point() + geom_jitter(height = 5) +
    geom_smooth() +
    labs(title="Edges in skeletons inferred from resamplings",
         x = "How often this edge was inferred in the skeletons",
         y = "Average Conditional Mutual Information for this edge")
  ggsave(paste("CI", consensus, nSkeletons, fraction, ".jpg", sep = ' '))
  # Plotting number of nodes in d-separation set
  as_tibble(miic.res$skeletons) %>%
    group_by(x, y) %>%
    summarise(count = n(),
              ai_vect = mean(as.numeric(ai_vect_n)),
              I = mean(as.numeric(I))) %>%
    ggplot(aes(x=count*100/nSkeletons, y=ai_vect, color = log(I))) + geom_jitter(height = 0.01) + geom_point() +
    geom_smooth() +
    labs(title="Edges in skeletons inferred from resamplings",
         x = "How often this edge was inferred in the skeletons",
         y = "Number of nodes in the d-separation set",
         fill = "Log of average conditional mutual information")
  ggsave(paste("d-sep", consensus, nSkeletons, fraction, ".jpg", sep = ' '))
  # Proportion of similarity
  final_network <- as.data.frame(cbind(miic.res$retained.edges.summary$x,
                                       miic.res$retained.edges.summary$y))
  colnames(final_network) <- colnames(miic.res$consensus_table) <- c('x', 'y')
  # Proportion of final network which is not in the consensus table
  n_edges_not_in_consensus <- nrow(dplyr::setdiff(final_network,
                                                  miic.res$consensus_table))
  n_total_final_network <- nrow(final_network)
  prop_edges_not_in_consensus <- n_edges_not_in_consensus*100/n_total_final_network
  prop_edges_in_consensus <- 100-prop_edges_not_in_consensus
  write(paste0(prop_edges_in_consensus, ',',
               consensus, ',',
               nSkeletons, ',',
               fraction),
        file="percent_similarity.txt",append=TRUE)
}

# Handling arguments
option_list = list(
  make_option(c("-k", "--nSkeletons"), type="integer", default=10,
              help=paste("number of skeletons inferred for consensus skeleton",
                         "[default=%default]."),
              metavar="integer"),
  make_option(c("-c", "--consensus"), type="integer", default=80,
              help=paste("frequency of occurrence for adding edge to consensus",
                         "network [default=%default]."),
              metavar="integer"),
  make_option(c("-f", "--fraction"), type="integer", default=90,
              help=paste("fraction that will be randomly sampled from original",
                         "data [default=%default]."),
              metavar="integer"),
  make_option(c("-s", "--seed"), type="integer", default=2019,
              help="seed for reproduction of results [default=%default].",
              metavar="integer"),
  make_option(c("-o", "--output_path"),
              help="output destination [default=%default].",
              metavar="string"),
  make_option(c("-e", "--edge_filtering"), type="integer",
              help="0: filter edges only from full dataset
                1: filter no edges
                [default=%default].",
              metavar="integer"),
  make_option(c("-n", "--n_samples"), type="integer",
              help="number of samples from benchmarking dataset
                [default=%default].",
              metavar="integer"),
  make_option(c("-d", "--dataset"),
              help="dataset to be used for benchmarking
                [default=%default].",
              metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Calling code with handled arguments
benchmark(nSkeletons = opt$nSkeletons,
          consensus = opt$consensus,
          fraction = opt$fraction,
          seed = opt$seed,
          edgeFiltering = opt$edge_filtering,
          dataset=opt$dataset,
          nSamples=opt$n_samples,
          outputPath=opt$output_path)
