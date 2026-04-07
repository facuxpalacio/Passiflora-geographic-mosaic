# Load datasets
### Degree of frugivory
rho_df <- read.table("bird_rho.txt", head = T)
rho <- rho_df$rho
names(rho) <- rho_df$code
length(rho)

### Local abundance
pool_df <- read.table("pops_pool.txt", head = T)
log_pool <- t(log(pool_df[,-c(1,2)] + 1))
dim(log_pool)

### Binary mask for the distribution of species across populations
presence_by_pop <- read.table("presence_by_pop.txt", head = T)
presence_by_pop <- t(presence_by_pop[, -c(1:2)])
dim(presence_by_pop)

# Check that species in log_pool are present in presence_by_pop
diff_mat <- presence_by_pop - 1*t((pool_df[, -c(1:2)]>0))
any(diff_mat == -1)

### Visits per population
setwd("visits x population")
files <- paste(paste("pop", 1:8, sep = "_"), ".txt", sep = "")
list_counts <- lapply(files, function(f) {
  dat <- read.table(f, header = TRUE)
  mat <- as.matrix(dat[,-1])
  rownames(mat) <- dat[,1]
  mat
})

# Model parameters
n_pops <- length(list_counts)
n_species_obs <- nrow(list_counts[[1]])
n_species_aug <- length(rho)
n_augmented <- n_species_aug - n_species_obs
n_plants <- sapply(list_counts, ncol)
max_reps <- max(n_plants)
species_names <- rho_df$code
names_sp_obs <- rownames(list_counts[[1]]) # observed species
names_sp_nonobs <- species_names[!species_names %in% names_sp_obs] # non-observed species

list_counts_aug <- lapply(list_counts, function(mat) {
  
  # pad with zeros for unobserved species
  # reorder species to match the other datasets
  padded <- rbind(mat,
                  matrix(0, nrow = n_augmented, ncol = ncol(mat)))
  rownames(padded) <- c(names_sp_obs, names_sp_nonobs)
  padded_ordered <- padded[order(rownames(padded)), ]
  
  return(padded_ordered)
})

# Check dimensions
sapply(list_counts_aug, dim)

# Create 3D array
C_array <- array(0, c(n_species_aug, max(n_plants), n_pops))

for(k in 1:n_pops){
  C_array[, 1:n_plants[k], k] <- list_counts_aug[[k]]
}

C <- C_array

# Make sure that C has dimensions [species, max(n_plants), pops]
max_plants <- max(n_plants)
C_full <- array(NA, dim = c(n_species_aug, max_plants, n_pops))

for(k in 1:n_pops){
  C_full[, 1:n_plants[k], k] <- C[, 1:n_plants[k], k]
}
C <- C_full

# Hierarchical model
cat("
model {

  #### Detection model (species-specific only)
  for (i in 1:n_species_aug) {
    logit(p[i]) <- eta[i]
    eta[i] ~ dnorm(mu_eta, tau_eta)
  }

  #### Ecological process
  for (k in 1:n_pops) {          # populations
    for (i in 1:n_species_aug) { # species

      # Regional pool abundance
      log(lambda_potential[i,k]) <- alpha[k] + beta[k]*log_pool[i,k]

      # Latent inclusion
      w_raw[i,k] ~ dbern(omega[k]*rho[i])

      # Enforced biogeographic constraint
      w[i,k] <- w_raw[i,k]*presence_by_pop[i,k]

      # Expected local abundance
      lambda_actual[i,k] <- w[i,k]*lambda_potential[i,k]

      # Realized abundance (Poisson)
      N_actual[i,k] ~ dpois(lambda_actual[i,k])

      # Observations on plants
      for (j in 1:n_plants[k]) {
        C[i,j,k] ~ dbin(p[i], N_actual[i,k])
      }
    }

    # Population-level parameters
    alpha[k] ~ dnorm(mu_alpha, tau_alpha)
    beta[k]  ~ dnorm(mu_beta,  tau_beta)
    omega[k] ~ dbeta(2, 2)
  }

  #### Hyperpriors
  mu_alpha ~ dnorm(0, 0.01)
  tau_alpha <- pow(sigma_alpha, -2)
  sigma_alpha ~ dunif(0, 10)

  mu_beta ~ dnorm(0, 0.01)
  tau_beta <- pow(sigma_beta, -2)
  sigma_beta ~ dunif(0, 10)

  mu_eta ~ dnorm(0, 0.1)
  tau_eta <- pow(sigma_eta, -2)
  sigma_eta ~ dunif(0, 5)

}
", file = "model.txt")

# No observed species is forced to be absent
for (i in 1:n_species_aug) {
  for (k in 1:n_pops) {
    
    obs <- C[i, , k]
    obs <- obs[!is.na(obs)]
    
    if (length(obs) > 0 && any(obs > 0)) {
      presence_by_pop[i, k] <- 1
    }
  }
}

# Data list for JAGS
data_list <- list(
  C = C, # 3D array: species × plants × populations
  log_pool = log_pool, # matrix species × populations
  rho = rho, # species vector
  n_species_aug = n_species_aug,
  n_pops = n_pops,
  n_plants = n_plants, # vector of length n_pops
  presence_by_pop = presence_by_pop
)

# Initial values
set.seed(99)
make_inits <- function() {
  
  omega_init <- runif(n_pops, 0.1, 0.9)
  
  # w_raw: MUST be 0 if presence_by_pop == 0
  w_raw_init <- matrix(0, n_species_aug, n_pops)
  
  for (i in 1:n_species_aug) {
    for (k in 1:n_pops) {
      
      obs <- C[i, , k]
      obs <- obs[!is.na(obs)]
      maxC <- ifelse(length(obs) == 0, 0, max(obs))
      
      if (maxC > 0) {
        w_raw_init[i,k] <- 1
      } else if (presence_by_pop[i,k] == 1) {
        w_raw_init[i,k] <- rbinom(1, 1, omega_init[k] * rho[i])
      } else {
        w_raw_init[i,k] <- 0
      }
    }
  }
  
  # N_actual: MUST be 0 if presence_by_pop == 0 OR w_raw == 0
  # N_actual: STRICTLY consistent with observed counts
  N_init <- matrix(0, n_species_aug, n_pops)
  
  for (i in 1:n_species_aug) {
    for (k in 1:n_pops) {
      
      obs <- C[i, , k]
      obs <- obs[!is.na(obs)]
      
      maxC <- ifelse(length(obs) == 0, 0, max(obs))
      
      if (presence_by_pop[i,k] == 1 && w_raw_init[i,k] == 1) {
        
        # must be >= max observed
        N_init[i,k] <- max(maxC, 1) + 50
        
      } else {
        
        # must be zero if not present
        if (maxC > 0) {
          stop(paste("Observed count > 0 where presence = 0 at", i, k))
        }
        
        N_init[i,k] <- 0
      }
    }
  }
  
  
  list(
    w_raw      = w_raw_init,
    N_actual   = N_init,
    eta        = rnorm(n_species_aug, 0, 1),
    alpha      = rnorm(n_pops, 0, 1),
    beta       = rnorm(n_pops, 0, 1),
    omega      = omega_init,
    mu_alpha   = 0,
    sigma_alpha = 1,
    mu_beta    = 0,
    sigma_beta  = 1,
    mu_eta     = 0,
    sigma_eta   = 1
  )
}

# Parameters monitored
params <- c(
  "w",
  "N_actual",
  "lambda_potential",
  "lambda_actual",
  "p",        # now p[i] only
  "eta",
  "mu_eta",
  "sigma_eta",
  "alpha",
  "beta",
  "omega",
  "mu_alpha",
  "sigma_alpha",
  "mu_beta",
  "sigma_beta"
)


# Model fit
library(R2jags)
fit <- jags(
  data = data_list,
  inits = make_inits,
  parameters.to.save = params,
  model.file = "model.txt",
  n.chains = 3,
  n.iter = 30000,
  n.burnin = 8000,
  n.thin = 5
)

# write.csv(fit$BUGSoutput$summary, "mcmc.csv")
# save(fit, file = "fit_R2jags.RData")
# load("fit_R2jags.RData")

# Model summary
# Gelman-Rubin for all parameters
rhat <- fit$BUGSoutput$summary[, "Rhat"]
hist(rhat)
max(rhat)

# Compare observations with true abundances
library(MCMCvis)
# Create list of matrices with true abundances, analogous to list_counts_aug
N_post <- MCMCchains(fit, params = "N_actual")
N_med <- apply(N_post, 2, mean)   
N_med <- matrix(N_med, nrow = n_species_aug, ncol = n_pops)

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

plot(obs_all, pred_all)
abline(a = 0, b = 1)


