# deconvolute DNA methylation profile

library(io)
library(ggplot2)
library(tidyr)
library(stringr)
library(deconvMe)
library(MeDeCom)
library(RefFreeEWAS)
library(clue)

#' function for filtering for N most variably-methylated probes
#'
#' @param beta   input beta matrix to be filtered
#' @param n      number of most variably-methylated probes to filter for
filter_variable_probes <- function(beta, n = 10000){
	beta.probevar <- apply(beta, 1, var)
	var.idx <- head(order(beta.probevar, decreasing = TRUE), n = n)
	beta.v <- beta[var.idx,]
}

#' fix label swapping in outputs of reference-free deconvolution algorithms
#' using correlation between predicted and ground truth as cost matrix 
#' to maximize via Hungarian algorithm
#'
#' @param predicted  N x K matrix of predicted composition OR J x K matrix of predicted profiles
#' @param actual     corresponding ground truth composition or compartment-specific profiles with
#'		     same dimensions as predicted matrix
fix_label_swap <- function(predicted, actual){
        cor.mat <- cor(actual, predicted)
        cor.mat.pos <- cor.mat + 1
        swap.idx <- as.integer(solve_LSAP(cor.mat.pos, maximum = TRUE))
        predicted.aligned <- predicted[, swap.idx, drop = FALSE]
}

nrmse <- function(pred, obs){
  sqrt(mean((pred - obs)^2))/mean(obs)
}

# function for plotting distribution of composition estimates as boxplots for each cell type
plot_composition <- function(w, method, mean_threshold = NULL){
	if(!is.null(mean_threshold)){
		w <- w[, colMeans(w, na.rm = TRUE) >= mean_threshold, drop = FALSE]
	}
	w.long <- stack(as.data.frame(w))
	ggplot(w.long, aes(x = ind, y = values, fill = ind)) +
	  geom_boxplot(width = 0.5, color = "black", alpha = 0.5, outlier.shape = NA) +
	  stat_summary(fun = median, geom = "point", shape = 21, size = 3, fill = "white") +
	  theme_classic() +
	  labs(
	    title = str_wrap(paste0(method, "-deconvoluted composition estimates for TCGA breast cancer data"), 40),
	    x = "Cell types",
	    y = "Composition estimates"
	  ) +
	  theme(
		legend.position = "none",
		axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)
	  )
}


# TODO: functions for running deconvolution methods

# methylresolver
run_mr <- function(beta, dataset, n_probes = NULL){
	out.dir.methylresolver <- file.path("output", dataset, "methylresolver")
	res.mr <- run_methylresolver(
	  beta_matrix = beta,
	  doPar = F,
	  numCores = 8,
	  alpha = seq(0.5, 0.9, by = 0.05),
	  absolute = TRUE,
	  purityModel = MethylResolver::RFmodel,
	  seed = 1,
	  cpg_subset = NULL
	)
	if(!is.null(n_probes)){
		probe_tag <- paste0("_n", n_probes)
	} else {
		probe_tag <- NULL
	}
	qwrite(res.mr, file.path(out.dir.methylresolver, paste0(dataset, probe_tag, "_methylresolver_composition.tsv")))
	res.mr
}

# RefFreeCellMix
run_rfcm <- function(beta, K){
	beta <- as.matrix(beta)
	mu <- RefFreeCellMixInitialize(Y = beta, K = K) 
	res.rfcm <- RefFreeCellMix(Y = beta, mu = mu, iter = 10)
}


