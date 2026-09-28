# Load packages
library(here)
library(R2jags)
library(MCMCvis)
library(coda)
library(ggplot2)
library(patchwork)

# ============================================================
# Load datasets
# ============================================================

### Degree of frugivory
rho_df <- read.table(
  here("data", "bird_rho.txt"),
  header = TRUE
)

rho <- rho_df$rho
names(rho) <- rho_df$code

length(rho)

### Regional abundance / species pool
pool_df <- read.table(
  here("data", "pops_pool.txt"),
  header = TRUE
)

log_pool <- t(log(pool_df[, -c(1, 2)] + 1))

dim(log_pool)


### Binary mask for the distribution of species across populations
presence_by_pop <- read.table(
  here("data", "presence_by_pop.txt"),
  header = TRUE
)

presence_by_pop <- t(presence_by_pop[, -c(1:2)])

dim(presence_by_pop)


# Check that species in log_pool are consistent with
# the geographic presence mask
diff_mat <- presence_by_pop -
  1 * t((pool_df[, -c(1, 2)] > 0))

any(diff_mat == -1)

# ============================================================
# Visits per population
# ============================================================

files <- paste(
  paste("pop", 1:8, sep = "_"),
  ".txt",
  sep = ""
)

list_counts <- lapply(files, function(f) {
  
  dat <- read.table(
    here("data", "visits x population", f),
    header = TRUE
  )
  
  mat <- as.matrix(dat[, -1])
  rownames(mat) <- dat[, 1]
  
  mat
})


# ============================================================
# Dimensions and species augmentation
# ============================================================

n_pops <- length(list_counts)

n_species_aug <- length(rho)

n_plants <- sapply(list_counts, ncol)

max_reps <- max(n_plants)

species_names <- names(rho)


# ============================================================
# Check that all observed species are included in rho
# ============================================================

all_observed_species <- unique(
  unlist(
    lapply(list_counts, rownames)
  )
)

missing_from_rho <- setdiff(
  all_observed_species,
  species_names
)

if (length(missing_from_rho) > 0) {
  
  stop(
    paste(
      "These observed species are missing from rho:",
      paste(missing_from_rho, collapse = ", ")
    )
  )
}


# ============================================================
# Augment count matrices
# species × plants
# ============================================================

list_counts_aug <- lapply(list_counts, function(mat) {
  
  out <- matrix(
    0,
    nrow = n_species_aug,
    ncol = ncol(mat),
    dimnames = list(
      species_names,
      colnames(mat)
    )
  )
  
  common_species <- intersect(
    rownames(mat),
    species_names
  )
  
  out[common_species, ] <- mat[common_species, ]
  
  out
})


# ============================================================
# Check dimensions
# ============================================================

sapply(
  list_counts_aug,
  dim
)


# ============================================================
# Create 3D observation array
# species × plants × populations
# ============================================================

C_array <- array(
  0,
  dim = c(
    n_species_aug,
    max(n_plants),
    n_pops
  )
)

for (k in 1:n_pops) {
  
  C_array[
    , 1:n_plants[k],
    k
  ] <- list_counts_aug[[k]]
}

C <- C_array


# Missing plant-level observations are NA
# beyond the actual number of plants in each population

C_full <- array(
  NA,
  dim = c(
    n_species_aug,
    max(n_plants),
    n_pops
  )
)

for (k in 1:n_pops) {
  
  C_full[
    , 1:n_plants[k],
    k
  ] <- C[
    , 1:n_plants[k],
    k
  ]
}

C <- C_full

# Species order in all datasets
head(rownames(C[, 1, 1]))
head(rownames(presence_by_pop))
head(names(rho))
head(rownames(log_pool))

identical(
  rownames(presence_by_pop),
  names(rho)
)

identical(
  rownames(presence_by_pop),
  rownames(log_pool)
)

identical(
  rownames(presence_by_pop),
  rownames(list_counts_aug[[1]])
)


# ============================================================
# Check consistency between observations and geographic mask
# ============================================================

for (i in 1:n_species_aug) {
  
  for (k in 1:n_pops) {
    
    obs <- C[i, , k]
    obs <- obs[!is.na(obs)]
    
    if (
      length(obs) > 0 &&
      any(obs > 0) &&
      presence_by_pop[i, k] == 0
    ) {
      
      stop(
        paste(
          "Observed visits for species",
          rownames(presence_by_pop)[i],
          "in population",
          k,
          "where presence_by_pop = 0"
        )
      )
    }
  }
}


# ============================================================
# JAGS model
# ============================================================

cat("
model {

  ############################################################
  # Detection model
  ############################################################

  for (i in 1:n_species) {

    logit(p[i]) <- eta[i]

    eta[i] ~ dnorm(mu_p, tau_p)

  }


  ############################################################
  # Frugivory-based occurrence constraint
  ############################################################

  for (i in 1:n_species) {

    M[i] <- step(rho[i] - 1.0E-10)

  }


  ############################################################
  # Ecological process
  ############################################################

  for (k in 1:n_pops) {

    for (i in 1:n_species) {

      lambda[i,k] <- presence_by_pop[i,k] * M[i] * exp(
        alpha +
        beta * log_pool_c[i,k] +
        gamma * rho_c[i]
      )


      ########################################################
      # Latent local abundance
      ########################################################

      N[i,k] ~ dpois(lambda[i,k])


      ########################################################
      # Observation model
      ########################################################

      for (j in 1:n_plants[k]) {

        C[i,j,k] ~ dbin(
          p[i],
          N[i,k]
        )

      }
    }
  }


  ############################################################
  # Global ecological parameters
  ############################################################

  alpha ~ dnorm(0, 0.01)

  beta ~ dnorm(0, 0.01)

  gamma ~ dnorm(0, 0.01)


  ############################################################
  # Hierarchical detection model
  ############################################################

  mu_p ~ dnorm(0, 0.1)

  tau_p <- pow(sigma_p, -2)

  sigma_p ~ dunif(0, 5)

}
", file = "Nmix_model.txt")


# ============================================================
# Data list for JAGS
# ============================================================

# Center covariates
mean_log_pool <- mean(log_pool)
mean_rho <- mean(rho)

log_pool_c <- log_pool - mean_log_pool
rho_c <- rho - mean_rho

data_list <- list(
  C = C,
  log_pool_c = log_pool_c,
  rho = rho, # original rho for M
  rho_c = rho_c, # centered rho for gamma
  n_species = n_species_aug,
  n_pops = n_pops,
  n_plants = n_plants,
  presence_by_pop = presence_by_pop
)


# ============================================================
# Initial values
# ============================================================

set.seed(2026)

make_inits <- function() {
  
  eta_init <- rnorm(
    n_species_aug,
    0,
    1
  )
  
  N_init <- matrix(
    1,
    nrow = n_species_aug,
    ncol = n_pops
  )
  
  for (i in 1:n_species_aug) {
    
    for (k in 1:n_pops) {
      
      obs <- C[i, , k]
      obs <- obs[!is.na(obs)]
      
      maxC <- ifelse(
        length(obs) == 0,
        0,
        max(obs)
      )
      
      if (
        presence_by_pop[i, k] == 1 &&
        rho[i] > 0
      ) {
        
        N_init[i, k] <- max(
          maxC + 1,
          1
        )
        
      } else {
        
        N_init[i, k] <- 0
        
      }
    }
  }
  
  list(
    N = N_init,
    eta = eta_init,
    alpha = rnorm(1, 0, 1),
    beta = rnorm(1, 0, 1),
    gamma = rnorm(1, 0, 1),
    mu_p = 0,
    sigma_p = 1
  )
}


# ============================================================
# Parameters monitored
# ============================================================

params <- c(
  
  # Ecological parameters
  "alpha",
  "beta",
  "gamma",
  
  # Expected local abundance
  "lambda",
  
  # Latent local abundance
  "N",
  
  # Detection
  "p",
  "eta",
  "mu_p",
  "sigma_p"
  
)


# ============================================================
# Model fit
# ============================================================

fit <- jags(
  
  data = data_list,
  
  inits = make_inits,
  
  parameters.to.save = params,
  
  model.file = "Nmix_model.txt",
  
  n.chains = 4,
  
  n.iter = 120000,
  
  n.burnin = 20000,
  
  n.thin = 10
  
)


# ============================================================
# Save / load
# ============================================================

 write.csv(
   fit$BUGSoutput$summary,
   "mcmcb.csv"
 )

 save(
   fit,
   file = "fit_R2bjags.RData"
 )

# load("fit_R2jags.RData")


# ============================================================
# MCMC diagnostics
# ============================================================

summary_fit <- fit$BUGSoutput$summary

rhat <- summary_fit[, "Rhat"]

hist(
  rhat,
  main = "Rhat",
  xlab = "Rhat"
)

max(
  rhat,
  na.rm = TRUE
)

Neff <- summary_fit[, "n.eff"]

hist(
  Neff,
  main = "Neff",
  xlab = "Neff"
)

max(
  Neff,
  na.rm = TRUE
)

mcmc <- as.mcmc(fit)

traceplot(
  mcmc[, c("alpha","beta","gamma")],
  smooth = FALSE
)

# ============================================================
# Posterior estimates of ecological parameters
# ============================================================

summary_fit[
  c(
    "alpha",
    "beta",
    "gamma"
  ),
]


# ============================================================
# Posterior latent abundance
# ============================================================

N_post <- MCMCchains(
  fit,
  params = "N"
)

N_mean <- apply(
  N_post,
  2,
  mean
)

N_latent <- matrix(
  N_mean,
  nrow = n_species_aug,
  ncol = n_pops
)

rownames(N_latent) <- rownames(list_counts_aug[[1]])

colnames(N_latent) <- pool_df$poblacion.ID

N_latent

# Compare with observed visits
N_med <- matrix(N_latent, nrow = n_species_aug, ncol = n_pops)

list_Nmed <- vector("list", n_pops)

for (k in 1:n_pops) {
  n_pl <- n_plants[k]           
  mat_k <- matrix(
    N_med[, k],                 
    nrow = n_species_aug,
    ncol = n_plants[k]
  )
  
  rownames(mat_k) <- rownames(list_counts_aug[[k]])
  colnames(mat_k) <- colnames(list_counts_aug[[k]])
  
  list_Nmed[[k]] <- mat_k
}

obs_all <- unlist(list_counts_aug)
pred_all <- unlist(list_Nmed)

plot(obs_all, pred_all) # should be above the 1:1 line
abline(a = 0, b = 1)
