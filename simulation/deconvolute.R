# perform deconvolution on simulated data
library(dplyr)
library(io)
library(ggplot2)
library(ggsci)

X11.options(type = "cairo")

source("../R/deconvolution.R")

#---
# load data

simdat.list <- list.files("data")

for(simdat in simdat.list){
	J <- sub(".*_j([0-9]+).*", "\\1", simdat)
	K <- sub(".*-k([0-9]+).*", "\\1", simdat)
	N <- sub(".*-n([0-9]+).*", "\\1", simdat)
	P <- sub(".*-p([0-9]+).*", "\\1", simdat)
	assign(paste0("simdat_j", J, "_k", K, "_n", N, "_p", P), readRDS(file.path("data", simdat)))
}

#---
# algorithm running functions

# run RefFreeCellMix, write output and plots
run_rfcm <- function(data, K, write = TRUE, plot = TRUE, seed = 123){
        beta <- as.matrix(t(data$mixed))
	set.seed(seed)
        mu <- RefFreeCellMixInitialize(Y = beta, K = K)
        res.rfcm <- RefFreeCellMix(Y = beta, mu = mu, iter = 10)
	if(write){
		J <- data$params$J
		K <- data$params$K
		N <- data$params$N
		p <- data$params$tumour_ploidy
		data_name <- paste0("sim-j", J, "-k", K, "-n", N, "-p", p)
		qwrite(res.rfcm, file.path("output", data_name, 
					    paste0("meth-decon_reffreecellmix_", data_name, ".rds")))
	}
	if(plot){
		evaluate_meth_decon(res.rfcm$Omega, data$composition, type = "composition", algo_name = "RefFreeCellMix", plot = TRUE, params = data$params)
	}
	res.rfcm
}

# run MeDeCom, write output and plots
run_mdc <- function(data, K, lambda = 1e-3, write = TRUE, plot = TRUE, seed = 123){
        beta <- as.matrix(t(data$mixed))
	set.seed(seed)
	res.medecom <- runMeDeCom(beta, K, c(0,10^(-5:-1)),
                           NINIT = 10, NFOLDS = 10, ITERMAX = 300, NCORES = 10)
	res.medecom.lmcs <- getLMCs(res.medecom, K = K, lambda = lambda)
	res.medecom.prop <- t(getProportions(res.medecom, K = K, lambda = lambda))
	if(write){
		J <- data$params$J
		K <- data$params$K
		N <- data$params$N
		p <- data$params$tumour_ploidy
		data_name <- paste0("sim-j", J, "-k", K, "-n", N, "-p", p)
		qwrite(res.medecom, file.path("output", data_name, 
					    paste0("meth-decon_medecom_", data_name, ".rds")))
	}
	if(plot){
		evaluate_meth_decon(res.medecom.lmcs, data$profiles, type = "profile", algo_name = "MeDeCom", plot = TRUE, params = data$params)
		evaluate_meth_decon(res.medecom.prop, data$composition, type = "composition", algo_name = "MeDeCom", plot = TRUE, params = data$params)
	}
	list(res = res.medecom,
	     lmc = res.medecom.lmcs,
	     prop = res.medecom.prop
	     )
}

# plotting function
est_plot <- function(predicted, actual, type = "composition", method, params){
	J <- params$J
        K <- params$K
        N <- params$N
        P <- params$tumour_ploidy
	prop.df0 <- pivot_longer(as.data.frame(predicted),
		     cols = everything(),
		     names_to = "compartment",
		     values_to = "estimate"
	)
	if(type == "composition"){
		title <- "Composition"
		labels <- "proportions"
	} else if(type == "profile"){
		title <- "Methylation beta profile"
		labels <- "beta values"
	}
	prop.df <- data.frame(prop.df0, actual = c(t(actual)))
	ggplot(prop.df, aes(x = estimate, y = actual, color = compartment)) +
	    geom_point(show.legend = TRUE) +
	    geom_abline() +
	    scale_color_manual(values = pal_npg("nrc", alpha = 0.7)(5)) +
	    ggtitle(str_wrap(paste0(title, " estimation using ", method, " on simulated data
	                            with J=", J, ", K=", K, ", N=", N, ", ploidy=", P), width = 50)) +
	    ylab(paste0("Actual ", labels)) +
	    xlab(paste0("Estimated ", labels)) +
	    theme_classic() +
	    scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
	    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
	    coord_fixed(ratio = 1)
}


# evaluate NRMSE and R2
evaluate_meth_decon <- function(predicted, actual, type = "composition", algo_name, plot = FALSE, params){
	predicted <- fix_label_swap(predicted, actual)
	colnames(predicted) <- paste0("compartment", 1:ncol(predicted))
	nrmse <- nrmse(c(predicted), c(actual))
	cor <- cor(c(predicted), c(actual))^2
	if(plot){
		J <- params$J
                K <- params$K
                N <- params$N
                p <- params$tumour_ploidy
                data_name <- paste0("sim-j", J, "-k", K, "-n", N, "-p", p)
		qdraw(est_plot(predicted, actual, type, algo_name, params),
		      file.path("output", data_name, paste0("meth-decon_", tolower(algo_name), "_", type, "-eval_", data_name, ".pdf")))
	}
	list(nrmse = nrmse,
	     cor = cor
	     )
}

#---

# run RefFreeCellMix
res.rfcm.sim1 <- run_rfcm(simdat_j1000_k5_n100_p4, K = 5)
res.rfcm.sim2 <- run_rfcm(simdat_j5000_k5_n100_p4, K = 5)
res.rfcm.sim3 <- run_rfcm(simdat_j1000_k5_n100_p2, K = 5)
res.rfcm.sim4 <- run_rfcm(simdat_j1000_k5_n100_p3, K = 5)

#---
# MeDeCom

res.medecom.sim1 <- run_mdc(simdat_j1000_k5_n100_p4, K = 5)
res.medecom.sim2 <- run_mdc(simdat_j5000_k5_n100_p4, K = 5)
res.medecom.sim3 <- run_mdc(simdat_j1000_k5_n100_p2, K = 5)
res.medecom.sim4 <- run_mdc(simdat_j1000_k5_n100_p3, K = 5)

