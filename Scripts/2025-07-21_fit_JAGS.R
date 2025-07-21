source("./Scripts/TurkeyMultiOccu_package_load.R")


## Load in habitat covariates ----------------------------------------------

#scaled_covs <- read.csv("./Data/2024-10-17_scaled_covariates.csv")
#covs <- read.csv("./Data/2024-10-17_covariates.csv")

## Process observation data ----------------------------------------------

# Turkey Detections -------------------------------------------------------

# load in the detection histories
det <- read.csv("./Data/TurkeyOccupancyReport.csv")

# site names
sites <- det$Site[1:77]

# keep only detection. For some reasonit comes with a weird last column
y_mat <- det[ ,c(3,6:(ncol(det)-1))]

# make this into weekly detection histories
n_split <- 7
n_start <- 2

# split seasons apart
season_split <- seq(n_split, ncol(y_mat), n_split) |>   
  map(~ select(y_mat[,-(1:(n_start-1))],(.-(n_split-1)):.))

# sum across rows to get total number of detections at each site for each season
weekly_det_list <- lapply(season_split, rowSums, na.rm = TRUE)
weekly_det <- do.call("cbind", weekly_det_list)

# find NA's
no_sample <- lapply(season_split, function(x){
  rowSums(is.na(x))
})
no_samples_mat <- do.call("cbind", no_sample)

# add to detections to identify sites that were not sampled during a week
weekly_det[no_samples_mat == 7] <- NA

# change values greater than 1 to 1 so we just have 1's and 0's
weekly_det[weekly_det > 1] <- 1

weekly_det_2 <- data.frame(season = det$Season, weekly_det)

# split into chunks based on season
tmp_list <- weekly_det_2 |>
  group_by(season) |>
  group_split()

# turn them into matrices and remove first column (season)
mat_list <- lapply(tmp_list, function(x) {
  as.matrix(x[,-1])
})

# turn this into an array for modeling in JAGS
y_array <- array(unlist(mat_list), dim = c(length(sites), 
                                           6, 
                                           n_distinct(det$Season)))


# Set up observation covariates -------------------------------------------

# create a detection covariate of number of days a camera was active per week

# our `no_sample_mat` was the number of NA's so 7 minus that would be the number
# of days a camera was active
cam_days <- 7 - no_samples_mat

cam_days_2 <- data.frame(season = det$Season, cam_days)

# split into chunks based on season
cam_list <- cam_days_2 |>
  group_by(season) |>
  group_split()

# turn them into matrices and remove first column (season)
cam_mats <- lapply(cam_list, function(x) {
  as.matrix(x[,-1])
})

# turn this into an array for modeling in JAGS
cam_array <- array(unlist(cam_mats), dim = c(length(sites), 
                                             6, 
                                             n_distinct(det$Season)))


# Prepare data for JAGS model ---------------------------------------------


# data list for model
data_list <- list(
  noccasion = dim(y_array)[2],
  nsite = dim(y_array)[1],
  nseason = dim(y_array)[3],
  ncovar = ncol(scaled_covs),
  X = as.matrix(scaled_covs),
  A_a = cam_array,
  y = y_array
)

# initial values
my_inits <- function(chain){
  gen_list <- function(chain = chain){
    list(
      z = matrix(1, ncol = data_list$nseason, nrow = data_list$nsite),
      psi0 = rnorm(1),
      #ridge = rnorm(data_list$ncovar),
      laplace = runif(data_list$ncovar, -3, 3),
      lambda = rgamma(1,0.1,0.01),
      pi = round(rbeta(data_list$ncovar,1,1)),
      pp=rbeta(1,1,1),
      a0 = rnorm(1),
      a1 = rnorm(1),
      theta = rnorm(1),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Wichmann-Hill",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}

# fit the model
my_mod <- runjags::run.jags(
  model = "./Scripts/2025-07-21-lasso_model_jags.R",
  monitor = c("theta", "psi0", "psi_beta", "a0", "a1", "pi", "lambda", "pp"),
  data = data_list,
  inits = my_inits,
  n.chains = 4,
  adapt = 1000,
  burnin = 25000,
  sample = 35000,
  thin = 4,
  method = "parallel"
)

# model diagnostics
my_sum <- summary(my_mod)
MCMCtrace(as.mcmc.list(my_mod))
plotAuto(as.matrix(as.mcmc.list(my_mod)), thin =8)

caterplot(
  my_mod,
  collapse = TRUE,
  reorder = FALSE,
  quantiles = list(outer=c(0.025,0.975),inner=c(0.10,0.90))
)


#saveRDS(my_mod, "./Results/2024-10-18_auto_occu_model_results.rds")


# Run Reduced Model -------------------------------------------------------

reduced_covs_names <- c("dist2road_scaled", 
                        "canopy_4000_scaled",
                        "dist2water_scaled",
                        "habitat_1000_scaled",
                        "log_dist_trail_scaled",
                        "imp_1000_scaled",
                        "imp_2000_scaled",
                        "imp_4000_scaled",
                        "dist2stream_scaled",
                        "rugged_2000_scaled",
                        "rugged_4000_scaled"
                        )

# data list for model
data_list <- list(
  noccasion = dim(y_array)[2],
  nsite = dim(y_array)[1],
  nseason = dim(y_array)[3],
  ncovar = ncol(scaled_covs[,reduced_covs_names]),
  X = as.matrix(scaled_covs[,reduced_covs_names]),
  A_a = cam_array,
  y = y_array
)

# initial values
my_inits <- function(chain){
  gen_list <- function(chain = chain){
    list(
      z = matrix(1, ncol = data_list$nseason, nrow = data_list$nsite),
      psi0 = rnorm(1),
      #ridge = rnorm(data_list$ncovar),
      psi_beta = runif(data_list$ncovar, -3, 3),
      lambda = rgamma(1,0.1,0.01),
      a0 = rnorm(1),
      a1 = rnorm(1),
      theta = rnorm(1),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Wichmann-Hill",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}

# fit the model
my_mod_reduced <- runjags::run.jags(
  model = "./Scripts/2025-07-21_reduced_lasso_model_jags.R",
  monitor = c("theta", "psi0", "psi_beta", "a0", "a1", "lambda"),
  data = data_list,
  inits = my_inits,
  n.chains = 4,
  adapt = 1000,
  burnin = 25000,
  sample = 35000,
  thin = 4,
  method = "parallel"
)

# model diagnostics
my_sum <- summary(my_mod_reduced)
#write.csv(my_sum, "./Results/2025-07-18_reduced_model_summary.csv")
MCMCtrace(as.mcmc.list(my_mod_reduced))
plotAuto(as.matrix(as.mcmc.list(my_mod_reduced)), thin =8)