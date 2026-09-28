# ============================================================
# Visits per population
# ============================================================

files <- paste(
  paste("pop", 1:10, sep = "_"),
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

# One value per site and species
comm_matrix <- do.call(
  cbind,
  lapply(list_counts, rowSums, na.rm = T)
)

colnames(comm_matrix) <- tools::file_path_sans_ext(basename(files))
comm_matrix <- t(comm_matrix)
X <- c("SM", "CP", "LZ", "MG", "LP", "EC", "ER", "TU", "SF", "RA")
comm_matrix <- data.frame(
  X = X,
  comm_matrix,
  check.names = FALSE
)

# Check distance between model-based estimates and observed visits
comm_latent <- read.csv(here("data", "comm_passionflower.csv"), head = T)
d_obs <- vegdist(comm_matrix[,-1], method = "bray")
d_latent <- vegdist(comm_latent[,-1], method = "bray")

summary(d_obs)
summary(d_latent)

cor(
  as.numeric(d_obs),
  as.numeric(d_latent),
  method = "pearson"
)

# Beta functional diversity (mean pairwise functional distances between assemblages)
birdtraits <- read.table(here("data", "frugivore_traits.txt"), head = T)
comm <- comm_matrix[,-1]
rownames(comm) <- comm_matrix$X

# Match species names
birdtraits <- birdtraits[birdtraits$code %in% colnames(comm), ]
rownames(birdtraits) <- birdtraits$sp

#### Body mass
Body_dist <- dist(scale(birdtraits[, "BodyMass.Value"]))
attr(Body_dist, "Labels") <- birdtraits$sp
hc <- hclust(Body_dist, method = "average")

colnames(comm) <- birdtraits$sp
beta_FD_Body <- comdist(comm, as.matrix(Body_dist),
                        abundance.weighted = TRUE)

vec_body <- as.vector(beta_FD_Body[lower.tri(beta_FD_Body)])
mean(vec_body,na.rm=T)
sd(vec_body,na.rm=T)
100*sd(vec_body,na.rm=T)/mean(vec_body,na.rm=T)

#### Foraging stratum
ForStrat_dist <- dist(scale(birdtraits[, c("ForStrat.ground",
                                           "ForStrat.understory",
                                           "ForStrat.midhigh",
                                           "ForStrat.canopy")]))
hc <- hclust(ForStrat_dist, method = "average")

colnames(comm) <- birdtraits$sp
beta_FD_ForStrat <- comdist(comm, as.matrix(ForStrat_dist),
                            abundance.weighted = TRUE)

vec_ForStrat <- as.vector(beta_FD_ForStrat[lower.tri(beta_FD_ForStrat)])
mean(vec_ForStrat,na.rm=T)
sd(vec_ForStrat,na.rm=T)
100*sd(vec_ForStrat,na.rm=T)/mean(vec_ForStrat,na.rm=T)

#### Fruit acquisition
fruit_acq_dist <- dist(scale(birdtraits[, c("fruit_aerial",
                                            "fruit_glean",
                                            "fruit_ground")]))
hc <- hclust(fruit_acq_dist, method = "average")

beta_FD_fruit_acq <- comdist(comm, as.matrix(fruit_acq_dist),
                             abundance.weighted = TRUE)

vec_acq <- as.vector(beta_FD_fruit_acq[lower.tri(beta_FD_fruit_acq)])
mean(vec_acq,na.rm=T)
sd(vec_acq,na.rm=T)
100*sd(vec_acq,na.rm=T)/mean(vec_acq,na.rm=T)


# Mean pairwise phylogenetic distances between assemblages
phylo_birdTree <- read.nexus(here("data", "birdTree_phylogeny_Hackett.nex"))
code_sp <- read.table(here("data", "code_birdTree_sp.txt"), head = T)

# Create consensus tree
cons_tree <- consensus(phylo_birdTree, p = 0.5)
cons_tree <- consensus.edges(phylo_birdTree, cons_tree,
                             method = "least.squares")

# Check and match name codes
sum(cons_tree$tip.label %in% code_sp$birdTree)
setdiff(cons_tree$tip.label, code_sp$birdTree)

lookup <- setNames(code_sp$code, code_sp$birdTree)
cons_tree$tip.label <- lookup[cons_tree$tip.label]

# Compute mean parwise phylogenetic distances between assemblage pairs
# change codes by names
sp_tree <- names(cons_tree$tip.label)
codes <- unname(cons_tree$tip.label)
cons_tree$tip.label <- sp_tree

colnames(comm)[colnames(comm) == "Pipraeidea_bonariensis"] <- 
  "Thraupis_bonariensis"
phyD <- comdist(comm, cophenetic(cons_tree), abundance.weighted = T)

vec_phy <- as.vector(phyD[lower.tri(phyD)])
mean(vec_phy,na.rm=T)
sd(vec_phy,na.rm=T)
100*sd(vec_phy,na.rm=T)/mean(vec_phy,na.rm=T)


# Generalized dissimilarity modeling
## Prepare response and predictors
##### Linear gradient matrix (response 1)
lin_grad <- read.csv(here("data", "beta_sim.csv"))
lin_grad[, -1] <- 1 - abs(lin_grad[, -1]) # Absolute values and dissimilarity
colnames(lin_grad) <- c("site", 1:nrow(lin_grad))
lin_grad[,1] <- 1:nrow(lin_grad)
dim(lin_grad) # first column = site

##### Nonlinear gradient matrix (response 2) 
nonlin_grad <- read.csv(here("data", "gamma_sim.csv"))
nonlin_grad[, -1] <- 1 - nonlin_grad[, -1]
colnames(nonlin_grad) <- c("site", 1:nrow(nonlin_grad))
nonlin_grad[,1] <- 1:nrow(nonlin_grad)
dim(nonlin_grad) # first column = site

##### Beta-diversity (Bray-Curtis) (predictor 1) 
prop_mat <- t(apply(comm_matrix[,-1], 1, function(x) x/sum(x)))
beta_div2 <- as.matrix(vegdist(prop_mat, method = "bray"))
rownames(beta_div2) <- colnames(beta_div2) <- comm_matrix$X

vec_comp <- as.vector(beta_div2[lower.tri(beta_div2)])
mean(vec_comp,na.rm=T)
sd(vec_comp,na.rm=T)
100*sd(vec_comp,na.rm=T)/mean(vec_comp,na.rm=T)

rownames(beta_div2) <- colnames(beta_div2)
dim(beta_div2)

# Order rows/columns alphabetically
ord <- order(rownames(beta_div2))
beta_div2 <- beta_div2[ord, ord]

# Create numeric site IDs
sites_char <- rownames(beta_div2)
sites_num  <- 1:length(sites_char)

# Create numeric-labeled distance matrix for GDM
beta_div2_numeric <- beta_div2
rownames(beta_div2_numeric) <- sites_num
colnames(beta_div2_numeric) <- sites_num

# Add "site" column (REQUIRED)
beta_div2_numeric <- data.frame(site = sites_num, 
                               beta_div2_numeric)

#### Functional diversity (mean body mass) (predictor 2)
FD_BM2 <- read.csv(here("data", "observed visits", "beta_FD_meanBodyMass2.csv"))[,-1]
rownames(FD_BM2) <- colnames(FD_BM2)

# Order rows/columns alphabetically
ord <- order(rownames(FD_BM2))
FD_BM2 <- FD_BM2[ord, ord]

# Create numeric site IDs
sites_char <- rownames(FD_BM2)
sites_num  <- 1:nrow(FD_BM2)

# Create numeric-labeled distance matrix for GDM
FD_BM2_numeric <- FD_BM2
rownames(FD_BM2_numeric) <- sites_num
colnames(FD_BM2_numeric) <- sites_num

# Add "site" column (REQUIRED)
FD_BM2_numeric <- data.frame(site = sites_num, FD_BM2_numeric)
dim(FD_BM2_numeric)

#### Functional diversity (foraging stratum) (predictor 3)
FD_strat2 <- read.csv(here("data", "observed visits", "beta_FD_ForStrat2.csv"))[,-1]
rownames(FD_strat2) <- colnames(FD_strat2)

# Order rows/columns alphabetically
ord <- order(rownames(FD_strat2))
FD_strat2 <- FD_strat2[ord, ord]

# Create numeric site IDs
sites_char <- rownames(FD_strat2)
sites_num  <- 1:nrow(FD_strat2)

# Create numeric-labeled distance matrix for GDM
FD_strat2_numeric <- FD_strat2
rownames(FD_strat2_numeric) <- sites_num
colnames(FD_strat2_numeric) <- sites_num

# Add "site" column (REQUIRED)
FD_strat2_numeric <- data.frame(site = sites_num, FD_strat2_numeric)
dim(FD_strat2_numeric)

#### Functional diversity (fruit harvesting) (predictor 4)
FD_harv2 <- read.csv(here("data", "observed visits", "beta_FD_fruit_acq2.csv"))[,-1]
rownames(FD_harv2) <- colnames(FD_harv2)

# Order rows/columns alphabetically
ord <- order(rownames(FD_harv2))
FD_harv2 <- FD_harv2[ord, ord]

# Create numeric site IDs
sites_char <- rownames(FD_harv2)
sites_num  <- 1:nrow(FD_harv2)

# Create numeric-labeled distance matrix for GDM
FD_harv2_numeric <- FD_harv2
rownames(FD_harv2_numeric) <- sites_num
colnames(FD_harv2_numeric) <- sites_num

# Add "site" column (REQUIRED)
FD_harv2_numeric <- data.frame(site = sites_num, FD_harv2_numeric)
dim(FD_harv2_numeric)

##### Population density (predictor 5)
populationID <- c("CP", "EC", "ER", "LZ", "LP", "MG", "RA",
                  "SF", "SM", "TU")
pop_density <- c(760, 134, 910, 4115, 239, 
                 0, 874, 1507, 695, 5445)
pop_df <- data.frame(populationID, pop_density)
log_pop <- log(pop_df$pop_density + 1)

# pairwise absolute log-differences
pop_dist <- as.data.frame(
  as.matrix(dist(log_pop, method = "euclidean"))
)
rownames(pop_dist) <- populationID
colnames(pop_dist) <- populationID
pop_dist <- data.frame(populationID, pop_dist)
pop_dist$populationID <- 1:nrow(pop_dist)
colnames(pop_dist) <- c("site", pop_dist$populationID)
dim(pop_dist) # first column = site

##### Geographic coordinates (predictor 6)
plant_data <- read.csv(here("data", "Passionflower_pops.csv"))

coords <- plant_data %>% 
  group_by(populationID) %>%
  summarise(
    Long = mean(as.numeric(longitude), na.rm = TRUE),
    Lat  = mean(as.numeric(latitude), na.rm = TRUE)
  )
coords$populationID <- c("CP", "EC", "ER", "LZ", "LP", "MG", "RA",
                         "SF", "SM", "TU")
coords[10, 2] <- -65.26979233 # Longitude Tucuman
coords[10, 3] <- -26.80470395 # Latitude Tucuman
coords[c(4,5), ] <- coords[c(5,4), ] # switch LZ and LP

###############################
### Prediction data for GDM ###
predData <- data.frame(
  site = sites_num,
  Long = coords$Long,
  Lat  = coords$Lat
)

# Coordinate projection
pts <- st_as_sf(predData, coords = c("Long", "Lat"), crs = 4326)

# Project to Argentina TM (EPSG:5347)
pts_proj <- st_transform(pts, 5347)

# Extract projected coordinates back to a data frame
coords <- st_coordinates(pts_proj)

predData_proj <- data.frame(
  site = predData$site,
  Long = coords[,1],
  Lat  = coords[,2]
)
geo_dist <- as.matrix(dist(coords))
rownames(geo_dist) <- populationID
colnames(geo_dist) <- populationID
geo_dist <- data.frame(populationID, geo_dist)
geo_dist$populationID <- 1:nrow(geo_dist)
colnames(geo_dist) <- c("site", geo_dist$populationID)
dim(geo_dist) # first column = site

##### Phylogenetic distance (predictor 7)
phyD2 <- read.csv(here("data", "observed visits", "phylo_dist2.csv"))[,-1]

# Order rows/columns alphabetically
ord <- order(rownames(phyD2))
phyD2 <- phyD2[ord, ord]

# Create numeric site IDs
sites_char <- rownames(phyD2)
sites_num  <- 1:length(sites_char)

# Create numeric-labeled distance matrix for GDM
phyD2_numeric <- phyD2
rownames(phyD2_numeric) <- sites_num
colnames(phyD2_numeric) <- sites_num

# Add "site" column (REQUIRED)
phyD2_numeric <- data.frame(site = sites_num, 
                           phyD2_numeric)

## Collinearity between predictors
matrices <- list(
  beta_div2  = beta_div2,
  FD_BM2     = FD_BM2,
  FD_strat2  = FD_strat2,
  FD_harv2   = FD_harv2,
  pop_dist  = pop_dist[,-1],
  geo_dist  = geo_dist[,-1],
  phyD2      = phyD2
)

matrices <- lapply(matrices, as.matrix)

n <- length(matrices)

corr_mantel <- matrix(
  NA,
  nrow = n,
  ncol = n,
  dimnames = list(names(matrices), names(matrices))
)

for (i in 1:n) {
  for (j in 1:n) {
    
    if (i == j) {
      corr_mantel[i, j] <- 1
      
    } else {
      
      res <- mantel(
        matrices[[i]],
        matrices[[j]],
        method = "pearson",
        permutations = 9999
      )
      
      corr_mantel[i, j] <- res$statistic
    }
  }
}

round(corr_mantel, 2)

## Combine predictors into gdm format
# Create initial site-pair table for linear gradients (bioFormat = 3)
gdmTab1 <- formatsitepair(
  bioData = lin_grad,     # biological distance matrix + site column
  bioFormat = 3,
  XColumn = "Long",
  YColumn = "Lat",
  siteColumn = "site",
  predData  = predData_proj      # coordinates table
)

# bioFormat = 4
## adds a predictor matrix to an existing site-pair table, in this case,
## predData needs to be provided, but is not actually used

# Combine predictors into a list
distPreds_list_obs <- list(geo_dist = geo_dist, 
                       beta_div2 = beta_div2_numeric,
                       FD_BM2 = FD_BM2_numeric, 
                       FD_strat2 = FD_strat2_numeric,
                       FD_harv2 = FD_harv2_numeric,
                       phylog2 = phyD2_numeric,
                       population = pop_dist) 

# Build the GDM site-pair table
gdmTab2_obs <- formatsitepair(
  bioData = gdmTab1,     # dissimilarity matrix with 'site' column
  bioFormat = 4,            # dissimilarity matrix
  siteColumn = "site",
  predData  = predData_proj,    # REQUIRED for format 4
  distPreds = distPreds_list_obs
)

## Linear selection: GDM fit and plotting
##### Fit GDM
gdmTabMod_lin_obs <- gdm(gdmTab2_obs, geo = F)

# Percent Deviance Explained: goodness-of-fit
# Intercept: expected dissimilarity between sites that do not differ in the predictors
# Summary of the fitted I-splines for each predictor, including the values of the 
# coefficients and their sum. 
# The sum indicates the amount of compositional turnover associated with that variable, 
# holding all other variables constant. I-spline summaries are order by coefficient sum. 
# Variables with all coefficients=0 have no relationship with the modeled biological pattern.
summary(gdmTabMod_lin_obs)

# Create initial site-pair table for nonlinear gradients (bioFormat = 3)
gdmTab3_obs <- formatsitepair(
  bioData = nonlin_grad,     # biological distance matrix + site column
  bioFormat = 3,
  XColumn = "Long",
  YColumn = "Lat",
  siteColumn = "site",
  predData  = predData_proj      # coordinates table
)

# bioFormat = 4
## adds a predictor matrix to an existing site-pair table, in this case,
## predData needs to be provided, but is not actually used

# Geographic distance matrix
# geo_dist <- as.matrix(data.frame(site = 1:nrow(coords), 
#                                 as.matrix(dist(coords))))

# build the GDM site-pair table
gdmTab4_obs <- formatsitepair(
  bioData = gdmTab3_obs,     # dissimilarity matrix with 'site' column
  bioFormat = 4,            # dissimilarity matrix
  siteColumn = "site",
  predData  = predData_proj,    # REQUIRED for format 4
  distPreds = distPreds_list_obs
)

## Nonlinear selection: GDM fit and plotting
##### Fit GDM
gdmTabMod_nonlin_obs <- gdm(gdmTab4_obs, geo = F)
summary(gdmTabMod_nonlin_obs)
