library(devtools)

load_all("~/projects/medotant")

dir.create("data")

simulate_meth_mixtures <- function(J, K, N, tumour_ploidy, seed = NULL, write = FALSE){
	if(!is.null(seed)){
		set.seed(seed)
	}
	profiles <- simulate_meth_profiles(K, J)
	composition <- matrix(0, nrow=N, ncol=K);
	composition[, K] <- runif(N, 0.5, 1.0);
	composition[, -K] <- runif(N*(K-1), 0, 0.2);
	composition <- composition / rowSums(composition);
	mixed <- simulate_mixed_meth_profiles(composition, tumour_ploidy, profiles);
	data <- list(profiles = profiles,
		     composition = composition,
		     mixed = mixed,
		     params = list(J = J,
				   K = K,
				   N = N,
				   tumour_ploidy = tumour_ploidy,
				   seed = seed
				   )
	)
	if(write){
		saveRDS(data, paste0("data/sim-tumour-beta_j", J, "-k", K, "-n", N, "-p", tumour_ploidy, ".rds"))
	}
	data
}


simulate_meth_mixtures(J = 1000, K = 5, N = 50, tumour_ploidy = 4, seed = 2345, write = TRUE)
simulate_meth_mixtures(J = 5000, K = 5, N = 50, tumour_ploidy = 4, seed = 2345, write = TRUE)

simulate_meth_mixtures(J = 1000, K = 5, N = 100, tumour_ploidy = 4, seed = 2345, write = TRUE)
simulate_meth_mixtures(J = 5000, K = 5, N = 100, tumour_ploidy = 4, seed = 2345, write = TRUE)

simulate_meth_mixtures(J = 1000, K = 2, N = 100, tumour_ploidy = 4, seed = 2345, write = TRUE)
simulate_meth_mixtures(J = 5000, K = 2, N = 100, tumour_ploidy = 4, seed = 2345, write = TRUE)

simulate_meth_mixtures(J = 1000, K = 5, N = 100, tumour_ploidy = 3, seed = 2345, write = TRUE)
simulate_meth_mixtures(J = 5000, K = 5, N = 100, tumour_ploidy = 3, seed = 2345, write = TRUE)

simulate_meth_mixtures(J = 1000, K = 5, N = 100, tumour_ploidy = 2, seed = 2345, write = TRUE)
simulate_meth_mixtures(J = 5000, K = 5, N = 100, tumour_ploidy = 2, seed = 2345, write = TRUE)


