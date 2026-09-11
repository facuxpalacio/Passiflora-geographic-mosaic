##################################################################
### SIMULATIONS TO TEST SIMPLIFIED N-MIXTURE MODEL
##################################################################

set.seed(2026)

library(R2jags)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)


# ---- Study design ----

n_species <- 83
n_pops <- 10
n_plants <- rep(20, n_pops)


# ---- True parameters ----

# Global intercept
alpha_true <- -0.3

# Global slope for regional abundance
beta_true <- 1

# Effect of degree of frugivory
gamma_true <- 1

# Degree of frugivory
rho <- rbeta(
  n_species,
  5,
  2
)

# Species-specific detection probability
mu_p_true <- 0.2
sigma_p_true <- 0.5

p_true <- plogis(
  rnorm(
    n_species,
    mean = mu_p_true,
    sd = sigma_p_true
  )
)


# ---- Regional abundance ----

log_pool <- matrix(
  rnorm(
    n_species * n_pops,
    mean = 1.5,
    sd = 0.8
  ),
  nrow = n_species,
  ncol = n_pops
)

log_pool_c <- sweep(
  log_pool,
  2,
  colMeans(log_pool)
)


# ---- Presence/absence mask ----
# Known regional presence/absence

presence_by_pop_true <- matrix(
  rbinom(
    n_species * n_pops,
    size = 1,
    prob = 0.85
  ),
  nrow = n_species,
  ncol = n_pops
)


# ---- Simulate latent abundance ----

lambda_true <- matrix(
  0,
  nrow = n_species,
  ncol = n_pops
)

N_true <- matrix(
  0,
  nrow = n_species,
  ncol = n_pops
)


for (i in 1:n_species) {
  
  for (k in 1:n_pops) {
    
    if (presence_by_pop_true[i, k] == 1) {
      
      lambda_true[i, k] <-
        exp(
          alpha_true +
            beta_true * log_pool_c[i, k] +
            gamma_true * rho[i]
        )
      
      N_true[i, k] <- rpois(
        1,
        lambda_true[i, k]
      )
      
    } else {
      
      lambda_true[i, k] <- 0
      N_true[i, k] <- 0
      
    }
  }
}


# ---- Simulate observed visits ----

C_sim <- array(
  NA,
  dim = c(
    n_species,
    max(n_plants),
    n_pops
  )
)


for (k in 1:n_pops) {
  
  for (i in 1:n_species) {
    
    for (j in 1:n_plants[k]) {
      
      if (presence_by_pop_true[i, k] == 1) {
        
        C_sim[i, j, k] <-
          rbinom(
            1,
            size = N_true[i, k],
            prob = p_true[i]
          )
        
      } else {
        
        C_sim[i, j, k] <- 0
        
      }
    }
  }
}


# ---- Data list ----

data_list_sim <- list(
  
  C = C_sim,
  log_pool_c = log_pool_c,
  rho = rho,
  n_species = n_species,
  n_pops = n_pops,
  n_plants = n_plants,
  presence_by_pop = presence_by_pop_true
)


# ---- Initial values ----

make_inits_sim <- function(
    C,
    presence,
    n_species,
    n_pops
) {
  
  N_init <- matrix(
    0,
    nrow = n_species,
    ncol = n_pops
  )
  
  
  for (i in 1:n_species) {
    
    for (k in 1:n_pops) {
      
      obs <- C[i, , k]
      obs <- obs[!is.na(obs)]
      
      maxC <- ifelse(
        length(obs) == 0,
        0,
        max(obs)
      )
      
      
      if (presence[i, k] == 1) {
        
        N_init[i, k] <-
          max(
            maxC,
            1
          ) +
          rpois(
            1,
            2
          )
        
      } else {
        
        N_init[i, k] <- 0
      }
    }
  }
  
  
  list(
    
    N = N_init,
    
    eta = rnorm(
      n_species,
      0,
      1
    ),
    
    alpha = rnorm(
      1,
      0,
      0.5
    ),
    
    beta = rnorm(
      1,
      1,
      0.3
    ),
    
    gamma = 0,
    
    mu_p = 0,
    
    sigma_p = 1
  )
}


inits_sim <- list(
  
  make_inits_sim(
    C_sim,
    presence_by_pop_true,
    n_species,
    n_pops
  ),
  
  make_inits_sim(
    C_sim,
    presence_by_pop_true,
    n_species,
    n_pops
  ),
  
  make_inits_sim(
    C_sim,
    presence_by_pop_true,
    n_species,
    n_pops
  )
)


###############################
### JAGS MODEL ################
###############################

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
  # Ecological process
  ############################################################

  for (k in 1:n_pops) {

    for (i in 1:n_species) {

      ########################################################
      # Regional abundance + frugivory
      ########################################################

      lambda[i,k] <- presence_by_pop[i,k] * exp(
        alpha +
        beta * log_pool_c[i,k] +
        gamma * rho[i]
      )


      ########################################################
      # Latent abundance
      ########################################################

      N[i,k] ~ dpois(lambda[i,k])


      ########################################################
      # Observations
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
  # Hyperpriors: detection
  ############################################################

  mu_p ~ dnorm(0, 0.1)

  tau_p <- pow(sigma_p, -2)

  sigma_p ~ dunif(0, 5)

}
", file = "Nmix_simple1.txt")


# ---- Parameters to monitor ----

params <- c(
  "alpha",
  "beta",
  "gamma",
  "p",
  "mu_p",
  "sigma_p",
  "N"
)


fit_sim <- jags(
  
  data = data_list_sim,
  
  inits = inits_sim,
  
  parameters.to.save = params,
  
  model.file = "Nmix_simple1.txt",
  
  n.chains = 3,
  n.iter = 20000,
  n.burnin = 5000,
  n.thin = 5
)

write.csv(
  fit_sim$BUGSoutput$summary,
  "mcmc_sim_simple_Nmix1.csv"
)

save(
  fit_sim,
  file = "fit_sim_simple_Nmix1.RData"
)

load("fit_sim_simple_Nmix1.RData")

# Convergence diagnostics -------------------------------------------------
# -------------------------------------------------------------------------
# Rhat and effective sample size: all monitored parameters
# -------------------------------------------------------------------------

metrics <- as.data.frame(fit_sim$BUGSoutput$summary)

metrics2 <- metrics %>%
  select(Rhat, n.eff) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Metrica",
    values_to = "Valor"
  )

# Diagnóstico específico para N[i,k]

metrics_N <- metrics %>%
  filter(grepl("^N\\[\\d+,\\d+\\]$", rownames(.)))


pct_neff <- mean(
  metrics_N$n.eff < 100,
  na.rm = TRUE
) * 100

pct_rhat <- mean(
  metrics_N$Rhat > 1.1,
  na.rm = TRUE
) * 100


# Histograma de Neff

p1 = ggplot(metrics_N, aes(x = n.eff)) +
  geom_histogram(
    bins = 30,
    fill = "steelblue",
    color = "black"
  ) +
  labs(
    title = "Neff for estimated total number of visits",
    subtitle = sprintf(
      "%.1f%% of posteriors with Neff < 100",
      pct_neff
    ),
    x = "Neff",
    y = "Freq"
  )

# Histograma de Rhat

p2 = ggplot(metrics_N, aes(x = Rhat)) +
  geom_histogram(
    bins = 30,
    fill = "tomato",
    color = "black"
  ) +
  labs(
    title = "Rhat estimated for total number of visits",
    subtitle = sprintf(
      "%.1f%% parameters Rhat > 1.1",
      pct_rhat
    ),
    x = "Rhat",
    y = "Freq"
  )


p1 + p2

# Revisar los Neff de N
summary(metrics_N$n.eff)

quantile(
  metrics_N$n.eff,
  probs = c(.01, .05, .10, .25, .50, .75, .90, .95, .99),
  na.rm = TRUE
)

mean(metrics_N$n.eff < 50, na.rm = TRUE) * 100
mean(metrics_N$n.eff < 100, na.rm = TRUE) * 100
mean(metrics_N$n.eff < 200, na.rm = TRUE) * 100
mean(metrics_N$n.eff < 500, na.rm = TRUE) * 100

metrics_N %>%
  arrange(n.eff) %>%
  head(30)

# cuantos N tienen Neff=0
metrics_N %>%
  summarise(
    n = n(),
    zero_sd = sum(sd == 0, na.rm = TRUE),
    prop_zero_sd = mean(sd == 0, na.rm = TRUE),
    n_eff_1 = sum(n.eff == 1, na.rm = TRUE),
    prop_neff_1 = mean(n.eff == 1, na.rm = TRUE)
  )
table(metrics_N$mean[metrics_N$sd == 0])

# Para cada especie × población:
# ¿cuántas detecciones hubo entre todas las plantas?

C_total <- matrix(0, n_species, n_pops)

for (i in 1:n_species) {
  for (k in 1:n_pops) {
    C_total[i, k] <- sum(
      C_sim[i, 1:n_plants[k], k],
      na.rm = TRUE
    )
  }
}

diag_N <- metrics_N %>%
  mutate(
    i = as.integer(sub("N\\[(\\d+),.*", "\\1", rownames(metrics_N))),
    k = as.integer(sub("N\\[\\d+,(\\d+)\\]", "\\1", rownames(metrics_N)))
  ) %>%
  mutate(
    C_total = mapply(
      function(i, k) C_total[i, k],
      i, k
    )
  )

table(
  N_zero_SD = diag_N$sd == 0,
  C_total_zero = diag_N$C_total == 0
)

diag_N %>%
  filter(sd == 0, C_total > 0) %>%
  select(i, k, mean, sd, C_total, Rhat, n.eff)

diag_N$C_max <- mapply(
  function(i, k) {
    max(
      C_sim[i, 1:n_plants[k], k],
      na.rm = TRUE
    )
  },
  diag_N$i,
  diag_N$k
)

table(
  diag_N$sd == 0 & diag_N$C_total > 0,
  diag_N$mean == diag_N$C_max
)

# Diagnostico posterior excluyendo los casos degenerados
bad_rhat %>%
  summarise(
    min_Rhat = min(Rhat),
    median_Rhat = median(Rhat),
    max_Rhat = max(Rhat),
    mean_neff = mean(n.eff),
    min_neff = min(n.eff)
  )


# PARAMETER RECOVERY
post <- fit_sim$BUGSoutput$sims.list

# -------------------------------------------------------------------------
# Function to calculate posterior mean and 95% CI
# -------------------------------------------------------------------------

build_recovery_df <- function(
    posterior,
    true_vals,
    param_name
) {
  
  posterior <- as.matrix(posterior)
  
  n <- ncol(posterior)
  
  mean_est <- apply(
    posterior,
    2,
    mean
  )
  
  lo <- apply(
    posterior,
    2,
    quantile,
    probs = 0.025
  )
  
  hi <- apply(
    posterior,
    2,
    quantile,
    probs = 0.975
  )
  
  data.frame(
    group = param_name,
    index = 1:n,
    true = as.vector(true_vals),
    mean = mean_est,
    lo = lo,
    hi = hi,
    covered = as.vector(true_vals) >= lo &
      as.vector(true_vals) <= hi
  )
}


# -------------------------------------------------------------------------
# alpha
# -------------------------------------------------------------------------

recovery_alpha <- build_recovery_df(
  posterior = post$alpha,
  true_vals = alpha_true,
  param_name = "alpha"
)


# -------------------------------------------------------------------------
# beta
# -------------------------------------------------------------------------

recovery_beta <- build_recovery_df(
  posterior = post$beta,
  true_vals = beta_true,
  param_name = "beta"
)


# -------------------------------------------------------------------------
# gamma
# -------------------------------------------------------------------------

gamma_mean <- mean(post$gamma)

gamma_lo <- quantile(
  post$gamma,
  0.025
)

gamma_hi <- quantile(
  post$gamma,
  0.975
)

recovery_gamma <- data.frame(
  group = "gamma",
  index = 1,
  true = gamma_true,
  mean = gamma_mean,
  lo = gamma_lo,
  hi = gamma_hi,
  covered =
    gamma_true >= gamma_lo &
    gamma_true <= gamma_hi
)


# -------------------------------------------------------------------------
# p
# -------------------------------------------------------------------------

recovery_p <- build_recovery_df(
  posterior = post$p,
  true_vals = p_true,
  param_name = "p"
)


# -------------------------------------------------------------------------
# N
# -------------------------------------------------------------------------

# post$N has dimensions:
# iterations x species x populations

N_est <- apply(
  post$N,
  c(2, 3),
  mean
)

N_lo <- apply(
  post$N,
  c(2, 3),
  quantile,
  probs = 0.025
)

N_hi <- apply(
  post$N,
  c(2, 3),
  quantile,
  probs = 0.975
)


recovery_N <- data.frame(
  
  group = "N",
  
  index = 1:length(N_true),
  
  true = as.vector(N_true),
  
  mean = as.vector(N_est),
  
  lo = as.vector(N_lo),
  
  hi = as.vector(N_hi),
  
  covered =
    as.vector(N_true) >= as.vector(N_lo) &
    as.vector(N_true) <= as.vector(N_hi)
)

# -------------------------------------------------------------------------
# mu_p
# -------------------------------------------------------------------------

mu_p_mean <- mean(post$mu_p)

mu_p_lo <- quantile(
  post$mu_p,
  0.025
)

mu_p_hi <- quantile(
  post$mu_p,
  0.975
)

recovery_mu_p <- data.frame(
  group = "mu_p",
  index = 1,
  true = mu_p_true,
  mean = mu_p_mean,
  lo = mu_p_lo,
  hi = mu_p_hi,
  covered =
    mu_p_true >= mu_p_lo &
    mu_p_true <= mu_p_hi
)


# -------------------------------------------------------------------------
# sigma_p
# -------------------------------------------------------------------------

sigma_p_mean <- mean(post$sigma_p)

sigma_p_lo <- quantile(
  post$sigma_p,
  0.025
)

sigma_p_hi <- quantile(
  post$sigma_p,
  0.975
)

recovery_sigma_p <- data.frame(
  group = "sigma_p",
  index = 1,
  true = sigma_p_true,
  mean = sigma_p_mean,
  lo = sigma_p_lo,
  hi = sigma_p_hi,
  covered =
    sigma_p_true >= sigma_p_lo &
    sigma_p_true <= sigma_p_hi
)

# -------------------------------------------------------------------------
# Combine
# -------------------------------------------------------------------------

recovery_all <- bind_rows(
  recovery_alpha,
  recovery_beta,
  recovery_gamma,
  recovery_mu_p,
  recovery_sigma_p,
  recovery_p,
  recovery_N
)


recovery_all$group <- factor(
  recovery_all$group,
  levels = c(
    "alpha",
    "beta",
    "gamma",
    "mu_p",
    "sigma_p",
    "p",
    "N"
  )
)


# -------------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------------

p_recovery <- ggplot(
  recovery_all,
  aes(
    x = true,
    y = mean,
    ymin = lo,
    ymax = hi,
    color = covered
  )
) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey50",
    linewidth = 1
  ) +
  
  geom_pointrange(
    alpha = 0.7,
    fatten = 1.5,
    linewidth = 0.7
  ) +
  
  facet_wrap(
    ~group,
    scales = "free",
    ncol = 3
  ) +
  
  scale_color_manual(
    values = c(
      "TRUE" = "steelblue",
      "FALSE" = "firebrick"
    )
  ) +
  
  labs(
    title = "Parameter recovery",
    x = "True value",
    y = "Posterior mean (95% CI)",
    color = "Within CI 95%"
  ) +
  
  theme_minimal() +
  
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold"
    ),
    axis.title = element_text(
      size = 15
    ),
    axis.text = element_text(
      size = 12
    ),
    legend.title = element_text(
      size = 14
    ),
    legend.text = element_text(
      size = 12
    ),
    strip.text = element_text(
      size = 14,
      face = "bold"
    )
  )

p_recovery

# Check posterior correlations --------------------------------------------

safe_cor <- function(x, y) {
  
  if (
    sd(x) == 0 ||
    sd(y) == 0
  ) {
    return(NA_real_)
  }
  
  cor(x, y)
}


# -------------------------------------------------------------------------
# N vs p, by species x population
# -------------------------------------------------------------------------

corr_N_p <- matrix(
  NA,
  nrow = n_species,
  ncol = n_pops
)


for (i in 1:n_species) {
  
  for (k in 1:n_pops) {
    
    corr_N_p[i, k] <- safe_cor(
      post$N[, i, k],
      post$p[, i]
    )
  }
}


df_N_p <- data.frame(
  index = 1:length(corr_N_p),
  species = rep(1:n_species, times = n_pops),
  population = rep(1:n_pops, each = n_species),
  correlation = as.vector(corr_N_p),
  pair = "N vs p"
)

# -------------------------------------------------------------------------
# Combine
# -------------------------------------------------------------------------

corr_summary <- df_N_p


corr_summary$pair <- factor(
  corr_summary$pair,
  levels = c("N vs p")
)

# -------------------------------------------------------------------------
# Numerical summary
# -------------------------------------------------------------------------

cat(
  "Mean absolute posterior correlation N-p:",
  mean(
    abs(df_N_p$correlation),
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Maximum absolute posterior correlation N-p:",
  max(
    abs(df_N_p$correlation),
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Percentage with |r| > 0.5:",
  mean(
    abs(df_N_p$correlation) > 0.5,
    na.rm = TRUE
  ) * 100,
  "%\n"
)

aggregate(
  correlation ~ pair,
  data = corr_summary,
  FUN = function(x) {
    
    c(
      mean_abs = mean(
        abs(x),
        na.rm = TRUE
      ),
      
      max_abs = max(
        abs(x),
        na.rm = TRUE
      ),
      
      pct_strong = mean(
        abs(x) > 0.5,
        na.rm = TRUE
      ) * 100
    )
  }
)


# -------------------------------------------------------------------------
# Boxplot + individual correlations
# -------------------------------------------------------------------------

p_corr <- ggplot(
  subset(
    df_N_p,
    !is.na(correlation)
  ),
  aes(
    x = "",
    y = correlation
  )
) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.4,
    size = 2,
    color = "steelblue"
  ) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.7
  ) +
  
  geom_hline(
    yintercept = c(-0.5, 0.5),
    linetype = "dotted",
    linewidth = 0.7
  ) +
  
  geom_boxplot(
    width = 0.35, alpha = 0,
    outlier.shape = NA
  ) +
  
  labs(
    title = "Posterior correlation between N and p",
    subtitle = "Dotted lines indicate |r| = 0.5",
    x = NULL,
    y = "Posterior correlation (Pearson's r)"
  ) +
  
  theme_minimal() +
  
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )

p_corr

##################################################
### INDEPENDENT DATASET VALIDATION
celtis <- read.csv("C:/RD/SXY.csv")
plants <- celtis[, c(1, 14:27)]

# Matriz de conteos
C_int <- as.matrix(plants[, -1])

# Número de especies y réplicas
n_species <- ncol(plants) - 1
n_pops <- 1
n_plants <- 158

# Matriz de conteos: especies x árboles
C_mat <- t(as.matrix(plants[, -1]))

# Agregar dimensión de población
C_intensive <- array(
  C_mat,
  dim = c(n_species, n_plants, n_pops)
)

dim(C_intensive)

# Geographic mask
presence_by_pop <- matrix(
  1,
  nrow = 14,
  ncol = 1
)

set.seed(0)

# -----------------------------
# Datos
# -----------------------------

C_all <- as.matrix(plants[, -1])

n_species <- ncol(C_all)
n_pops <- 1
n_plants <- 20

# Seleccionar 20 plantas al azar
idx <- sample(1:nrow(C_all), n_plants)

C <- array(
  C_all[idx, ],
  dim = c(n_species, n_plants, n_pops)
)

# Todas las especies fueron observadas en esta población
presence_by_pop <- matrix(
  1,
  nrow = n_species,
  ncol = n_pops
)


# Frugivory degree and abundance pool
df <- read.table("C:/RD/traits_abu_pool.txt",head=T)
rho <- df$frugivory_degree
log_pool <- log(df$abu_pool)
log_pool_c <- log_pool - mean(log_pool, na.rm = TRUE)

log_pool_c_20 <- matrix(
  log_pool_c,
  nrow = n_species,
  ncol = 1
)

rho_20 <- rho

# Model
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
  # Ecological process
  ############################################################

  for (k in 1:n_pops) {

    for (i in 1:n_species) {

      lambda[i,k] <- presence_by_pop[i,k] * exp(
        alpha +
        beta * log_pool_c[i,k] +
        gamma * rho[i]
      )

      N[i,k] ~ dpois(lambda[i,k])

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
  # Hyperpriors: detection
  ############################################################

  mu_p ~ dnorm(0, 0.1)

  tau_p <- pow(sigma_p, -2)

  sigma_p ~ dunif(0, 5)

}
", file = "Nmix_simple1.txt")

# Make inits
# Máximo observado por especie en las 20 plantas
N_init <- apply(C[, , 1], 1, max) + 50

inits <- function() {
  
  list(
    N = matrix(
      apply(C[, , 1], 1, max) +
        sample(10:30, n_species, replace = TRUE),
      nrow = n_species,
      ncol = n_pops
    )
  )
}

# Model fit
data_jags <- list(
  C = C,
  n_species = n_species,
  n_pops = n_pops,
  n_plants = n_plants,
  presence_by_pop = presence_by_pop,
  log_pool_c = log_pool_c_20,
  rho = rho_20
)

params <- c(
  "alpha",
  "beta",
  "gamma",
  "p",
  "mu_p",
  "sigma_p",
  "N"
)

fit_20 <- jags(
  data = data_jags,
  parameters.to.save = params,
  model.file = "Nmix_simple1.txt",
  n.chains = 3,
  n.iter = 30000,
  n.burnin = 10000,
  n.thin = 2
)

N_summary <- data.frame(
  species = colnames(C_all),
  N_20 = apply(N_post[, , 1], 2, mean),
  lower_20 = apply(N_post[, , 1], 2, quantile, probs = 0.025),
  upper_20 = apply(N_post[, , 1], 2, quantile, probs = 0.975),
  N_158 = N_158
)

N_summary

ggplot(N_summary, aes(x = N_158, y = N_20)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_errorbar(
    aes(ymin = lower_20, ymax = upper_20),
    width = 0
  ) +
  geom_point(size = 3) +
  geom_text(
    aes(label = species),
    hjust = -0.05,
    vjust = 0.5,
    size = 3
  ) +
  labs(
    x = "Observed abundance across 158 plants",
    y = "Estimated abundance from 20 plants"
  ) +
  theme_classic()