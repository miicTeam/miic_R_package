library(dplyr)
library(ggplot2)
library(purrr)
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
  final_network <- data.frame(cbind(miic.res$retained.edges.summary$x,
                                    miic.res$retained.edges.summary$y),
                              stringsAsFactors = FALSE)
  colnames(final_network) <- colnames(miic.res$consensus_table) <- c('x', 'y')

  # Calculating F-score
  # MIIC edge filtering network: miic_ef
  # Consensus:                   miic_cons
  # Real network:                real_network
  miic_cons <- data.frame(miic.res$consensus_table, stringsAsFactors = FALSE)
  miic_cons[] <- lapply(miic_cons, as.character)
  miic_ef <- data.frame(final_network, stringsAsFactors = FALSE)
  # Obtaining real network
  real_network <- bn.net(myBnNet)$arcs 
  real_network <- as.data.frame(real_network)
  colnames(real_network) <- c('x', 'y')
  real_network[] <- lapply(real_network, as.character)
  # TP are edges that are in the inferred network and in the real network
  # FP are edges that are in the inferred network and not in the real network
  # FN are edges that are not in the inferred network and are in the real network
  # precision = TP/(TP+FP)
  # recall = TP/(TP+FN)
  # F-score = 2*(Prec.Rec)/(Prec+Rec)

  is_in <- function(x, y, target) {
    any(
      (target[[1]] == x & target[[2]] == y) | (target[[2]] == x & target[[1]] == y)
    )
  }

  # TP of miic_ef vs real_network
  miic_ef %>%
    mutate(flag = map2_lgl(x, y, is_in, target = real_network)) %>%
    pluck(., "flag") %>%
    sum -> TP_ef

  # TP of miic_cons vs real_network
  miic_cons %>%
    mutate(flag = map2_lgl(x, y, is_in, target = real_network)) %>%
    pluck(., "flag") %>%
    sum -> TP_cons

  # FP of miic_ef vs real_network
  miic_ef %>%
    mutate(flag = !map2_lgl(x, y, is_in, target = real_network)) %>%
    pluck(., "flag") %>%
    sum -> FP_ef

  # FP of miic_cons vs real_network 
  miic_cons %>%
    mutate(flag = !map2_lgl(x, y, is_in, target = real_network)) %>%
    pluck(., "flag") %>%
    sum -> FP_cons

  # FN of miic_ef vs real_network
  real_network %>%
    mutate(flag = !map2_lgl(x, y, is_in, target = miic_ef)) %>%
    pluck(., "flag") %>%
    sum -> FN_ef

  # FN of miic_cons vs real_network
  real_network %>%
    mutate(flag = !map2_lgl(x, y, is_in, target = miic_cons)) %>%
    pluck(., "flag") %>%
    sum -> FN_cons

  # Precision of miic_ef
  prec_ef <- TP_ef/(TP_ef+FP_ef)
  # Precision of miic_cons
  prec_cons <- TP_cons/(TP_cons+FP_cons)
  # Recall of miic_ef
  rec_ef <- TP_ef/(TP_ef+FN_ef)
  # Recall of miic_cons
  rec_cons <- TP_cons/(TP_cons+FN_cons)

  # Calculate F-score miic_ef
  Fscore_ef <- 2*(prec_ef*rec_ef)/(prec_ef+rec_ef)
  # Calculate F-score miic_cons
  Fscore_cons <- 2*(prec_cons*rec_cons)/(prec_cons+rec_cons)
  # Calculate F-scores of miic_cons minus miic_ef
  fscore_diff <- Fscore_cons - Fscore_ef

  if (!file.exists('percent_similarity.txt')) {
    write(paste0('True Positive EF', ',', 'True Positive Consensus', ',',
                 'False Positive EF', ',', 'False Positive Consensus', ',',
                 'False Negative EF', ',', 'False Negative Consensus', ',',
                 'Recall EF', ',', 'Recall Consensus', ',',
                 'Precision EF', ',', 'Precision Consensus', ',',
                 'Fscore EF',  ',', 'Fscore Consensus', ',', 'Fscore Cons-EF',
                 ',', 'Number of samples', ',','Consensus %', ',',
                 'Number of Skeletons', ',', 'Fraction'), file="percent_similarity.txt")
  }

  write(paste0(TP_ef, ',',
               TP_cons, ',',
               FP_ef, ',',
               FP_cons, ',',
               FN_ef, ',',
               FN_cons, ',',
               rec_ef, ',',
               rec_cons, ',',
               prec_ef, ',',
               prec_cons, ',',
               Fscore_ef, ',',
               Fscore_cons, ',',
               fscore_diff, ',',
               nSamples, ',',
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
