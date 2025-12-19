library(SpatialDownscaling)
library(ggplot2)
library(reshape2)
library(abind)
library(gridExtra)
library(abind)
library(terra)
library(reshape2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# HELPER METHODS

plot_variable <- function(data, time_idx = 1, title = "",
                          var_range = NULL) {
  data_subset <- data[, , time_idx]
  data_melted <- melt(data_subset)
  colnames(data_melted) <- c("Longitude", "Latitude", "Value")

  ggplot(data_melted, aes(y = Longitude, x = Latitude)) +
    geom_raster(aes(fill = Value)) +
    scale_fill_gradientn(
      colours = heat.colors(10), name = title,
      limits = var_range
    ) +
    theme_minimal() +
    coord_fixed()
}

# Function to compare coarse, predicted, and true fine-scale data
compare_methods <- function(coarse_data, predicted_data, true_data,
                            time_idx, title, var_range = NULL,
                            method_name = "Method") {
  # Plot coarse data
  plots <- list()

  for (t_idx in time_idx) {
    coarse_subset <- coarse_data[, , t_idx]
    coarse_melted <- melt(coarse_subset)
    colnames(coarse_melted) <- c("Longitude", "Latitude", "Value")
    
    p1 <- ggplot(coarse_melted, aes(y = Longitude, x = Latitude)) +
      geom_raster(aes(fill = Value)) +
      scale_fill_viridis_c(option = "B", name = title,
        limits = var_range, guide = "none"
      ) +
      theme_minimal() +
      coord_fixed()

    plots[[length(plots) + 1]] <- p1

    # Plot predicted data
    pred_subset <- predicted_data[, , t_idx]
    pred_melted <- melt(pred_subset)
    colnames(pred_melted) <- c("Longitude", "Latitude", "Value")

    p2 <- ggplot(pred_melted, aes(y = Longitude, x = Latitude)) +
      geom_raster(aes(fill = Value)) +
      scale_fill_viridis_c(option = "B", name = title,
        limits = var_range, guide = "none"
      ) +
      theme_minimal() +
      coord_fixed()

    plots[[length(plots) + 1]] <- p2

    # Plot true fine data
    true_subset <- true_data[, , t_idx]
    true_melted <- melt(true_subset)
    colnames(true_melted) <- c("Longitude", "Latitude", "Value")

    p3 <- ggplot(true_melted, aes(y = Longitude, x = Latitude)) +
      geom_raster(aes(fill = Value)) +
      scale_fill_viridis_c(option = "B", name = title,
        limits = var_range
      ) +
      theme_minimal() +
      coord_fixed()

    plots[[length(plots) + 1]] <- p3
  }

  # Same plot for the plots list
  gridExtra::grid.arrange(grobs = plots,
    ncol = 3,
    widths = c(0.76, 0.76, 1)
  )
}

# Function to compare multiple downscaling methods
compare_all_methods <- function(coarse_data, bcsd_preds,
                                srdrn_preds, unet_preds, true_data,
                                time_idx, title, var_range = NULL) {
  # Plot coarse data
  coarse_subset <- coarse_data[, , time_idx]
  coarse_melted <- melt(coarse_subset)
  colnames(coarse_melted) <- c("Longitude", "Latitude", "Value")

  p1 <- ggplot(coarse_melted, aes(y = Longitude, x = Latitude)) +
    geom_raster(aes(fill = Value)) +
    scale_fill_gradientn(
      colours = heat.colors(10), name = title,
      limits = var_range
    ) +
    labs(title = "Coarse Input") +
    theme_minimal() +
    coord_fixed()

  # Plot BCSD predictions
  bcsd_subset <- bcsd_preds[, , time_idx]
  bcsd_melted <- melt(bcsd_subset)
  colnames(bcsd_melted) <- c("Longitude", "Latitude", "Value")

  p2 <- ggplot(bcsd_melted, aes(y = Longitude, x = Latitude)) +
    geom_raster(aes(fill = Value)) +
    scale_fill_gradientn(
      colours = heat.colors(10), name = title,
      limits = var_range
    ) +
    labs(title = "BCSD Prediction") +
    theme_minimal() +
    coord_fixed()

  # Plot SRCRN predictions
  srdrn_subset <- srdrn_preds[, , time_idx]
  srdrn_melted <- melt(srdrn_subset)
  colnames(srdrn_melted) <- c("Longitude", "Latitude", "Value")

  p3 <- ggplot(srdrn_melted, aes(y = Longitude, x = Latitude)) +
    geom_raster(aes(fill = Value)) +
    scale_fill_gradientn(
      colours = heat.colors(10), name = title,
      limits = var_range
    ) +
    labs(title = "SRCRN Prediction") +
    theme_minimal() +
    coord_fixed()

  # Plot U-Net predictions
  unet_subset <- unet_preds[, , time_idx]
  unet_melted <- melt(unet_subset)
  colnames(unet_melted) <- c("Longitude", "Latitude", "Value")

  p4 <- ggplot(unet_melted, aes(y = Longitude, x = Latitude)) +
    geom_raster(aes(fill = Value)) +
    scale_fill_gradientn(
      colours = heat.colors(10), name = title,
      limits = var_range
    ) +
    labs(title = "U-Net Prediction") +
    theme_minimal() +
    coord_fixed()

  # Plot true fine data
  true_subset <- true_data[, , time_idx]
  true_melted <- melt(true_subset)
  colnames(true_melted) <- c("Longitude", "Latitude", "Value")

  p5 <- ggplot(true_melted, aes(y = Longitude, x = Latitude)) +
    geom_raster(aes(fill = Value)) +
    scale_fill_gradientn(
      colours = heat.colors(10), name = title,
      limits = var_range
    ) +
    labs(title = "True Fine Data") +
    theme_minimal() +
    coord_fixed()

  # Arrange the plots in a grid (with the coarse data and true fine data on the top row)
  grid.arrange(
    p1, p5,
    p2, p3, p4,
    ncol = 2, nrow = 3,
    layout_matrix = rbind(c(1, 2), c(3, 4), c(5, NA)),
    top = grid::textGrob(paste(title, "- Method Comparison"),
      gp = grid::gpar(fontsize = 16)
    )
  )
}

# Function to calculate evaluation metrics
calculate_metrics <- function(predicted, actual) {
  # Mean Absolute Error
  mae <- mean(abs(predicted - actual), na.rm = TRUE)

  # Root Mean Squared Error
  rmse <- sqrt(mean((predicted - actual)^2, na.rm = TRUE))

  # Coefficient of Determination (R-squared)
  ss_tot <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
  ss_res <- sum((actual - predicted)^2, na.rm = TRUE)
  r_squared <- 1 - (ss_res / ss_tot)

  # Structural Similarity Index (simplified version)
  # Calculate mean and variance for actual and predicted
  # Min-max scaling:
  scaled_actual <- (actual - min(actual, na.rm = TRUE)) / (max(actual, na.rm = TRUE) - min(actual, na.rm = TRUE))
  scaled_predicted <- (predicted - min(predicted, na.rm = TRUE)) / (max(predicted, na.rm = TRUE) - min(predicted, na.rm = TRUE))
  mu_actual <- mean(scaled_actual, na.rm = TRUE)
  mu_pred <- mean(scaled_predicted, na.rm = TRUE)
  var_actual <- var(as.vector(scaled_actual), na.rm = TRUE)
  var_pred <- var(as.vector(scaled_predicted), na.rm = TRUE)

  # Calculate covariance
  cov_actual_pred <- cov(as.vector(scaled_actual), as.vector(scaled_predicted),
    use = "complete.obs"
  )

  # Constants for stability
  c1 <- 0.01^2
  c2 <- 0.03^2

  # Calculate SSIM
  ssim <- ((2 * mu_actual * mu_pred + c1) * (2 * cov_actual_pred + c2)) /
    ((mu_actual^2 + mu_pred^2 + c1) * (var_actual + var_pred + c2))


  # Kling-Gupta efficiency
  valid <- complete.cases(predicted, actual)
  predicted <- predicted[valid]
  actual  <- actual[valid]
  
  # Means and standard deviations
  mu_pred <- mean(predicted)
  mu_obs  <- mean(actual)
  sd_pred <- sd(predicted)
  sd_obs  <- sd(actual)
  
  # Components
  r <- cor(predicted, actual)                # linear correlation
  alpha <- sd_pred / sd_obs          # variability ratio
  beta  <- mu_pred / mu_obs          # bias ratio
  
  # Kling–Gupta Efficiency
  kge <- 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)

  return(list(MAE = mae, RMSE = rmse, R_squared = r_squared, SSIM = ssim, KGE = kge))
}

# Function to create a table comparing metrics across methods
compare_metrics <- function(bcsd_preds, srdrn_preds, unet_preds, true_data) {
  # Calculate metrics for each variable and method
  metrics_bcsd_humidity <- calculate_metrics(
    bcsd_preds[, , 1, ],
    true_data[, , 1, ]
  )
  metrics_bcsd_temp <- calculate_metrics(
    bcsd_preds[, , 2, ],
    true_data[, , 2, ]
  )
  metrics_bcsd_precip <- calculate_metrics(
    bcsd_preds[, , 3, ],
    true_data[, , 3, ]
  )

  metrics_srdrn_humidity <- calculate_metrics(
    srdrn_preds[, , 1, ],
    true_data[, , 1, ]
  )
  metrics_srdrn_temp <- calculate_metrics(
    srdrn_preds[, , 2, ],
    true_data[, , 2, ]
  )
  metrics_srdrn_precip <- calculate_metrics(
    srdrn_preds[, , 3, ],
    true_data[, , 3, ]
  )

  metrics_unet_humidity <- calculate_metrics(
    unet_preds[, , 1, ],
    true_data[, , 1, ]
  )
  metrics_unet_temp <- calculate_metrics(
    unet_preds[, , 2, ],
    true_data[, , 2, ]
  )
  metrics_unet_precip <- calculate_metrics(
    unet_preds[, , 3, ],
    true_data[, , 3, ]
  )

  # Create data frames for each metric
  mae_df <- data.frame(
    Variable = c("Humidity", "Temperature", "Precipitation"),
    BCSD = c(
      metrics_bcsd_humidity$MAE, metrics_bcsd_temp$MAE,
      metrics_bcsd_precip$MAE
    ),
    SRCRN = c(
      metrics_srdrn_humidity$MAE, metrics_srdrn_temp$MAE,
      metrics_srdrn_precip$MAE
    ),
    UNet = c(
      metrics_unet_humidity$MAE, metrics_unet_temp$MAE,
      metrics_unet_precip$MAE
    )
  )

  rmse_df <- data.frame(
    Variable = c("Humidity", "Temperature", "Precipitation"),
    BCSD = c(
      metrics_bcsd_humidity$RMSE, metrics_bcsd_temp$RMSE,
      metrics_bcsd_precip$RMSE
    ),
    SRCRN = c(
      metrics_srdrn_humidity$RMSE, metrics_srdrn_temp$RMSE,
      metrics_srdrn_precip$RMSE
    ),
    UNet = c(
      metrics_unet_humidity$RMSE, metrics_unet_temp$RMSE,
      metrics_unet_precip$RMSE
    )
  )

  r2_df <- data.frame(
    Variable = c("Humidity", "Temperature", "Precipitation"),
    BCSD = c(
      metrics_bcsd_humidity$R_squared, metrics_bcsd_temp$R_squared,
      metrics_bcsd_precip$R_squared
    ),
    SRCRN = c(
      metrics_srdrn_humidity$R_squared, metrics_srdrn_temp$R_squared,
      metrics_srdrn_precip$R_squared
    ),
    UNet = c(
      metrics_unet_humidity$R_squared, metrics_unet_temp$R_squared,
      metrics_unet_precip$R_squared
    )
  )

  ssim_df <- data.frame(
    Variable = c("Humidity", "Temperature", "Precipitation"),
    BCSD = c(
      metrics_bcsd_humidity$SSIM, metrics_bcsd_temp$SSIM,
      metrics_bcsd_precip$SSIM
    ),
    SRCRN = c(
      metrics_srdrn_humidity$SSIM, metrics_srdrn_temp$SSIM,
      metrics_srdrn_precip$SSIM
    ),
    UNet = c(
      metrics_unet_humidity$SSIM, metrics_unet_temp$SSIM,
      metrics_unet_precip$SSIM
    )
  )

  # Return list of data frames
  return(list(MAE = mae_df, RMSE = rmse_df, R_squared = r2_df, SSIM = ssim_df))
}
# shift_resample: realign a 3D array (lon, lat, time) onto a fixed 120x120 grid
shift_resample <- function(arr, method = "bilinear") {
  # arr: 3D array [lon, lat, time]

  nt   <- dim(arr)[3]

  # Define source grid (130x130 full domain)
  lon_src <- as.numeric(colnames(arr))
  lat_src <- as.numeric(rownames(arr))

  lon_tgt <- seq(6, 19.0, by = 0.1)
  lat_tgt <- seq(47.9, 35.0, by = -0.1)

  out <- array(NA, dim = c(length(lat_tgt), length(lon_tgt), nt))
  colnames(out) <- lon_tgt
  rownames(out) <- lat_tgt

  # Loop over time slices
  for (t in 1:nt) {
    r_src <- rast(arr[, , t])  # note: terra expects [x,y]
    ext(r_src) <- ext(min(lon_src), max(lon_src), min(lat_src), max(lat_src))
    crs(r_src) <- "EPSG:4326"

    r_tgt <- rast(ncols = length(lon_tgt), nrows = length(lat_tgt),
                  xmin = min(lon_tgt), xmax = max(lon_tgt),
                  ymin = min(lat_tgt), ymax = max(lat_tgt),
                  crs = "EPSG:4326")
    
    r_aligned <- resample(r_src, r_tgt, method = method)
    out[, , t] <- as.matrix(r_aligned, wide = TRUE)
  }
  return(out)
}


################# MAIN ANALYSIS ####################

data_all <- array(NA, dim = c(120, 120, 0))
lon_min <- 6.4
lon_max <- 18.3
lat_min <- 35.5
lat_max <- 47.4
for (year in 2018:2023) {
    load(paste0("pollution_grid_data/o3_", year, ".rda"))
    if (!(year %in% c(2018, 2019, 2021))) {
        y <- shift_resample(y, method = "bilinear")
    }
    lon_names <- as.numeric(colnames(y))
    lat_names <- as.numeric(rownames(y))
    lon_inds <- which(lon_names >= lon_min & lon_names <= lon_max)
    lat_inds <- which(lat_names >= lat_min & lat_names <= lat_max)
    cropped_y <- y[lon_inds, lat_inds, ]
    data_all <- abind(data_all, cropped_y, along = 3)
}

dim(data_all)

plot_variable(data_all, time_idx = 1642, title = "Ozone (µg/m³)", var_range = c(30, 125))

# Study area plot

world <- ne_countries(scale = "large", returnclass = "sf")

bbox <- st_as_sfc(st_bbox(c(xmin = lon_min, xmax = lon_max,
                            ymin = lat_min, ymax = lat_max),
                          crs = st_crs(world)))

ggplot() +
  geom_sf(data = world, fill = "grey90") +
  geom_sf(data = bbox, fill = NA, color = "red", size = 1) +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
  theme_minimal()

block_mean_downsample <- function(arr, factor=4) {
  nx <- dim(arr)[1]; ny <- dim(arr)[2]; nt <- dim(arr)[3]
  nx_out <- nx / factor; ny_out <- ny / factor
  out <- array(NA, dim = c(nx_out, ny_out, nt))
  for (i in 1:nt) {
    r_hr <- rast(arr[, , i])           # high-res (single layer)
    r_lr <- aggregate(r_hr, fact=factor, fun=mean)
    out[, , i] <- as.array(r_lr)
  }
  out
}
dimnames_coarse1 <- names(data_all[seq(1, 120, by = 4), 1, 1])
dimnames_coarse2 <- names(data_all[1, seq(1, 120, by = 4), 1])
dimnames_coarse3 <- dimnames(data_all)[[3]]
coarse_data <- block_mean_downsample(data_all, factor=4)
dimnames(coarse_data) <- list(dimnames_coarse1, dimnames_coarse2, dimnames_coarse3)

time_points <- seq(1, dim(coarse_data)[3])

plot_variable(coarse_data, time_idx = 1642, title = "Ozone (µg/m³)", var_range = c(30, 130))

test_inds <- which(time_points > 1991)
train_data_coarse <- coarse_data[, , -test_inds]
train_data_fine <- data_all[, , -test_inds]
train_times <- time_points[-test_inds]
test_data_coarse <- coarse_data[, , test_inds]
test_data_fine <- data_all[, , test_inds]
test_times <- time_points[test_inds]

val_inds <- which(train_times > 1500)
val_data_coarse <- train_data_coarse[, , val_inds]
val_data_fine <- train_data_fine[, , val_inds]
val_times <- train_times[val_inds]
train_val_data_coarse <- train_data_coarse[, , -val_inds]
train_val_data_fine <- train_data_fine[, , -val_inds]
train_val_times <- train_times[-val_inds]

# Statistical baseline BCSD
bcsd_model <- bcsd(train_data_coarse[, , ],
  train_data_fine[, , ],
  n_quantiles = 200)

predictions_bcsd <- predict(bcsd_model, test_data_coarse[, , ])
compare_methods(test_data_coarse, predictions_bcsd, test_data_fine,
  time_idx = 1,
  title = "Ozone (µg/m³)", var_range = c(30, 120),
  method_name = "BCSD"
)

bcsd_metrics <- calculate_metrics(
  predictions_bcsd,
  test_data_fine
)

######## VALIDATION ########

seed <- 05122025
# Baseline srdrn
srdrn_model <- srdrn(train_val_data_coarse[, , ],
  train_val_data_fine[, , ],
  val_coarse_data = val_data_coarse[, , ],
  val_fine_data = val_data_fine[, , ],
  use_batch_norm = FALSE,
  epochs = 80, batch_size = 32,
  learning_rate = 0.0005,
  validation_split = 0,
  seed = seed
)

# SRDRN with time

seed <- 05122025
srdrn_model_time <- srdrn(train_val_data_coarse[, , ],
  train_val_data_fine[, , ],
  val_coarse_data = val_data_coarse[, , ],
  val_fine_data = val_data_fine[, , ],
  val_time_points = val_times,
  use_batch_norm = FALSE,
  epochs = 80, batch_size = 32,
  validation_split = 0,
  temporal_layers = c(32, 64, 128, 256),
  time_points = train_val_times,
  cyclical_period = 365.25,
  seed = seed
)

# Cos-sin transformation for seasonal component

seed <- 05122025
srdrn_model_time_cossin <- srdrn(train_val_data_coarse[, , ],
  train_val_data_fine[, , ],
  val_coarse_data = val_data_coarse[, , ],
  val_fine_data = val_data_fine[, , ],
  val_time_points = val_times,
  use_batch_norm = FALSE,
  epochs = 80, batch_size = 32,
  temporal_layers = c(32, 64, 128, 256),
  validation_split = 0,
  learning_rate = 0.0005,
  time_points = train_val_times,
  cyclical_period = 365.25,
  cos_sin_time = TRUE,
  seed = seed
)

library(ggplot2)
epoch_seq <- 2:80
loss_data <- data.frame(
  Epoch = rep(epoch_seq, 3),
  Loss = c(srdrn_model$history$metrics$val_loss[2:80], srdrn_model_time$history$metrics$val_loss[2:80], srdrn_model_time_cossin$history$metrics$val_loss[2:80]),
  Type = rep(c("SRDRN (baseline)", "SRDRN (RBF)", "SRDRN (sinusoidal)"), each = length(epoch_seq))
)
ggplot(loss_data, aes(x = Epoch, y = Loss, color = Type)) +
  labs(x = "Epoch", y = "Validation loss") +
  geom_line(linetype = "dashed") +
  # Add smoothed lines in addition to raw lines:
  geom_smooth(se = FALSE, method = "loess", span = 0.3) +
  theme_minimal()

# UNet


seed <- 30082025
unet_model <- unet(
  train_val_data_coarse[, , ],
  train_val_data_fine[, , ],
  val_coarse_data = val_data_coarse[, , ],
  val_fine_data = val_data_fine[, , ],
  initial_filters = 16,
  filters = c(32, 64, 128),
  kernel_sizes = list(c(3, 3), c(3, 3), c(3, 3)),
  batch_size = 32,
  use_batch_norm = FALSE,
  validation_split = 0,
  learning_rate = 0.0005,
  callbacks = list(),
  epochs = 60,
  seed = seed,
  verbose = 1
)

# UNet with RBF
seed <- 30082025
unet_model_time <- unet_downscale(
  train_val_data_coarse[, , ],
  train_val_data_fine[, , ],
  time_points = train_val_times,
  val_coarse_data = val_data_coarse[, , ],
  val_fine_data = val_data_fine[, , ],
  val_time_points = val_times,
  cyclical_period = 365.25,
  initial_filters = 16,
  filters = c(32, 64, 128),
  kernel_sizes = list(c(3, 3), c(3, 3), c(3, 3)),
  temporal_layers = c(32, 64, 128, 256),
  batch_size = 32,
  use_batch_norm = FALSE,
  validation_split = 0,
  learning_rate = 0.0005,
  callbacks = list(),
  seed = seed,
  epochs = 60,
  verbose = 1
)


# Cos-sin transformation for time
seed <- 30082025
unet_model_time_cossin <- unet_downscale(
  train_val_data_coarse[, , ],
  train_val_data_fine[, , ],
  time_points = train_val_times,
  val_coarse_data = val_data_coarse[, , ],
  val_fine_data = val_data_fine[, , ],
  val_time_points = val_times,
  cyclical_period = 365.25,
  initial_filters = 16,
  filters = c(32, 64, 128),
  temporal_layers = c(32, 64, 128, 256),
  kernel_sizes = list(c(3, 3), c(3, 3), c(3, 3)),
  batch_size = 32,
  use_batch_norm = FALSE,
  validation_split = 0,
  learning_rate = 0.0005,
  epochs = 60,
  verbose = 1,
  cos_sin_transform = TRUE
)

epoch_seq <- 2:60
loss_data <- data.frame(
  Epoch = rep(epoch_seq, 3),
  Loss = c(unet_model$history$metrics$val_loss[2:60], unet_model_time$history$metrics$val_loss[2:60], unet_model_time_cossin$history$metrics$val_loss[2:60]),
  Type = rep(c("UNet (baseline)", "UNet (RBF)", "UNet (sinusoidal)"), each = length(epoch_seq))
)
ggplot(loss_data, aes(x = Epoch, y = Loss, color = Type)) +
  geom_line(linetype = "dashed") +
  # Add smoothed lines in addition to raw lines:
  geom_smooth(se = FALSE, method = "loess", span = 0.3) +
  labs(x = "Epoch", y = "Validation loss") +
  theme_minimal()


######### TRAINING FINAL MODELS ##########

# Baseline SRDRN

seed <- 05122025
srcnn_model <- srdrn(train_data_coarse[, , ],
  train_data_fine[, , ],
  val_input_data = val_data_coarse[, , ],
  val_target_data = val_data_fine[, , ],
  use_batch_norm = FALSE,
  learning_rate = 0.0005,
  epochs = 80, batch_size = 32,
  validation_split = 0,
  seed = seed
)

pred_batch_size <- 32
for (i in 1:ceiling(dim(test_data_coarse)[3] / pred_batch_size)) {
  batch_inds <- ((i - 1) * pred_batch_size + 1):min(i * pred_batch_size, dim(test_data_coarse)[3])
  cat("Predicting batch", i, "with indices", batch_inds, "\n")
  if (i == 1) {
    predictions_srcnn <- drop(predict(
      srcnn_model,
      test_data_coarse[, , batch_inds],
      time_points = test_times[batch_inds]
    ))
  } else {
    batch_preds <- drop(predict(
      srcnn_model,
      test_data_coarse[, , batch_inds],
      time_points = test_times[batch_inds]
    ))
    predictions_srcnn <- abind(predictions_srcnn, batch_preds, along = 3)
  }
}

save(predictions_srcnn, file = "predictions/predictions_srcnn_large_conv.RData")
load("predictions/predictions_srcnn_large_conv.RData")

compare_methods(test_data_coarse, predictions_srcnn, test_data_fine,
  time_idx = 1,
  title = "Ozone (µg/m³)", var_range = c(30, 120),
  method_name = "SRCRN"
)

srcnn_metrics <- calculate_metrics(
  predictions_srcnn,
  test_data_fine
)

pred_erros_srcnn <- predictions_srcnn - test_data_fine

# SRCNN with time

seed <- 05122025
srcnn_model_time <- srdrn(train_data_coarse[, , ],
  train_data_fine[, , ],
  use_batch_norm = FALSE,
  epochs = 80, batch_size = 32,
  validation_split = 0,
  temporal_layers = c(32, 64, 128, 256),
  learning_rate = 0.0005,
  time_points = train_times,
  cyclical_period = 365.25,
  seed = seed
)

pred_batch_size <- 32
for (i in 1:ceiling(dim(test_data_coarse)[3] / pred_batch_size)) {
  batch_inds <- ((i - 1) * pred_batch_size + 1):min(i * pred_batch_size, dim(test_data_coarse)[3])
  cat("Predicting batch", i, "with indices", batch_inds, "\n")
  if (i == 1) {
    predictions_srcnn_time <- drop(predict(
      srcnn_model_time,
      test_data_coarse[, , batch_inds],
      time_points = test_times[batch_inds]
    ))
  } else {
    batch_preds <- drop(predict(
      srcnn_model_time,
      test_data_coarse[, , batch_inds],
      time_points = test_times[batch_inds]
    ))
    predictions_srcnn_time <- abind(predictions_srcnn_time, batch_preds, along = 3)
  }
}
save(predictions_srcnn_time, file = "predictions/predictions_srcnn_time_conv.RData")
load("predictions/predictions_srcnn_time_conv.RData")

compare_methods(test_data_coarse, predictions_srcnn_time, test_data_fine,
  time_idx = 1,
  title = "Ozone (µg/m³)", var_range = c(30, 120),
  method_name = "SRCRN"
)

srcnn_metrics_time <- calculate_metrics(
  predictions_srcnn_time,
  test_data_fine
)

pred_erros_srcnn_time <- predictions_srcnn_time - test_data_fine

# Cos-sin transformation for seasonal component

seed <- 05122025
srcnn_model_time_cossin <- srdrn(train_data_coarse[, , ],
  train_data_fine[, , ],
  use_batch_norm = FALSE,
  epochs = 80, batch_size = 32,
  temporal_layers = c(32, 64, 128, 256),
  validation_split = 0,
  time_points = train_times,
  learning_rate = 0.0005,
  cyclical_period = 365.25,
  cos_sin_time = TRUE,
  seed = seed
)

pred_batch_size <- 32
for (i in 1:ceiling(dim(test_data_coarse)[3] / pred_batch_size)) {
  batch_inds <- ((i - 1) * pred_batch_size + 1):min(i * pred_batch_size, dim(test_data_coarse)[3])
  cat("Predicting batch", i, "with indices", batch_inds, "\n")
  if (i == 1) {
    predictions_srcnn_time_cossin <- drop(predict(
      srcnn_model_time_cossin,
      test_data_coarse[, , batch_inds],
      time_points = test_times[batch_inds]
    ))
  } else {
    batch_preds <- drop(predict(
      srcnn_model_time_cossin,
      test_data_coarse[, , batch_inds],
      time_points = test_times[batch_inds]
    ))
    predictions_srcnn_time_cossin <- abind(predictions_srcnn_time_cossin, batch_preds, along = 3)
  }
}
save(predictions_srcnn_time_cossin, file = "predictions/predictions_srcnn_time_cossin_large_conv.RData")
load("predictions/predictions_srcnn_time_cossin_large_conv.RData")

compare_methods(test_data_coarse, predictions_srcnn_time_cossin, test_data_fine,
  time_idx = 1,
  title = "Ozone (µg/m³)", var_range = c(30, 120),
  method_name = "SRCRN"
)

srcnn_metrics_time_cossin <- calculate_metrics(
  predictions_srcnn_time_cossin,
  test_data_fine
)

pred_erros_srcnn_time_cossin <- predictions_srcnn_time_cossin - test_data_fine

pred_errors_srcnn_agg <- apply(abs(pred_erros_srcnn), c(1,2), mean)
pred_errors_srcnn_agg <- array(pred_errors_srcnn_agg, dim = c(dim(pred_errors_srcnn_agg), 1))
colnames(pred_errors_srcnn_agg) <- colnames(pred_erros_srcnn)
rownames(pred_errors_srcnn_agg) <- rownames(pred_erros_srcnn)
pred_errors_srcnn_time_agg <- apply(abs(pred_erros_srcnn_time), c(1,2), mean)
pred_errors_srcnn_time_agg <- array(pred_errors_srcnn_time_agg, dim = c(dim(pred_errors_srcnn_time_agg), 1))
colnames(pred_errors_srcnn_time_agg) <- colnames(pred_erros_srcnn_time)
rownames(pred_errors_srcnn_time_agg) <- rownames(pred_erros_srcnn_time)
pred_errors_srcnn_time_cossin_agg <- apply(abs(pred_erros_srcnn_time_cossin), c(1,2), mean)
pred_errors_srcnn_time_cossin_agg <- array(pred_errors_srcnn_time_cossin_agg, dim = c(dim(pred_errors_srcnn_time_cossin_agg), 1))
colnames(pred_errors_srcnn_time_cossin_agg) <- colnames(pred_erros_srcnn_time_cossin)
rownames(pred_errors_srcnn_time_cossin_agg) <- rownames(pred_erros_srcnn_time_cossin)

compare_methods(pred_errors_srcnn_agg, pred_errors_srcnn_time_agg, pred_errors_srcnn_time_cossin_agg,
  time_idx = 1,
  title = "MAE over time", var_range = c(0, 3),
  method_names = c("SRDRN (baseline)", "SRDRN (RBF)", "SRDRN (sinusoidal)")
)


# Baseline UNet


seed <- 30082025
unet_model <- unet(
  train_data_coarse[, , ],
  train_data_fine[, , ],
  initial_filters = 16,
  filters = c(32, 64, 128),
  kernel_sizes = list(c(3, 3), c(3, 3), c(3, 3)),
  batch_size = 32,
  use_batch_norm = FALSE,
  validation_split = 0,
  learning_rate = 0.0005,
  callbacks = list(),
  epochs = 40,
  seed = seed,
  verbose = 1
)

pred_batch_size <- 32
for (i in 1:ceiling(dim(test_data_coarse)[3] / pred_batch_size)) {
  batch_inds <- ((i - 1) * pred_batch_size + 1):min(i * pred_batch_size, dim(test_data_coarse)[3])
  cat("Predicting batch", i, "with indices", batch_inds, "\n")
  if (i == 1) {
    predictions_unet <- drop(predict(
      unet_model,
      test_data_coarse[, , batch_inds]
      #time_points = test_times[batch_inds]
    ))
  } else {
    batch_preds <- drop(predict(
      unet_model,
      test_data_coarse[, , batch_inds]
      #time_points = test_times[batch_inds]
    ))
    predictions_unet <- abind(predictions_unet, batch_preds, along = 3)
  }
}
save(predictions_unet, file = "predictions/predictions_unet_large_conv_40.RData")
load("predictions/preditions_unet.RData")

compare_methods(test_data_coarse, predictions_unet, test_data_fine,
  time_idx = 200,
  title = "Ozone (µg/m³)", var_range = c(0, 120),
  method_name = "UNet"
)

unet_metrics <- calculate_metrics(
  predictions_unet,
  test_data_fine
)
pred_erros_unet <- predictions_unet - test_data_fine

# UNet with RBF
seed <- 30082025
unet_model_time <- unet(
  train_data_coarse[, , ],
  train_data_fine[, , ],
  time_points = train_times,
  cyclical_period = 365.25,
  initial_filters = 16,
  filters = c(32, 64, 128),
  kernel_sizes = list(c(3, 3), c(3, 3), c(3, 3)),
  batch_size = 32,
  temporal_layers = c(32, 64, 128, 256),
  use_batch_norm = FALSE,
  validation_split = 0,
  learning_rate = 0.0005,
  callbacks = list(),
  seed = seed,
  epochs = 40,
  verbose = 1
)

pred_batch_size <- 32
for (i in 1:ceiling(dim(test_data_coarse)[3] / pred_batch_size)) {
  batch_inds <- ((i - 1) * pred_batch_size + 1):min(i * pred_batch_size, dim(test_data_coarse)[3])
  cat("Predicting batch", i, "with indices", batch_inds, "\n")
  if (i == 1) {
    predictions_unet_time <- drop(predict(
      unet_model_time,
      test_data_coarse[, , batch_inds],
      time_points = test_times[batch_inds]
    ))
  } else {
    batch_preds <- drop(predict(
      unet_model_time,
      test_data_coarse[, , batch_inds],
      time_points = test_times[batch_inds]
    ))
    predictions_unet_time <- abind(predictions_unet_time, batch_preds, along = 3)
  }
}

save(predictions_unet_time, file = "predictions/predictions_unet_time_large_conv_40.RData")
load("predictions/predictions_unet_time_large_conv_40.RData")

compare_methods(test_data_coarse, predictions_unet_time, test_data_fine,
  time_idx = 1,
  title = "Ozone (µg/m³)", var_range = c(30, 120),
  method_name = "UNet"
)

unet_metrics_time <- calculate_metrics(
  predictions_unet_time,
  test_data_fine
)
pred_errors_unet_time <- predictions_unet_time - test_data_fine

# Cos-sin transformation for time
seed <- 30082025
unet_model_time_cossin <- unet(
  train_data_coarse[, , ],
  train_data_fine[, , ],
  time_points = train_times,
  cyclical_period = 365.25,
  initial_filters = 16,
  filters = c(32, 64, 128),
  temporal_layers = c(32, 64, 128, 256),
  kernel_sizes = list(c(3, 3), c(3, 3), c(3, 3)),
  batch_size = 32,
  use_batch_norm = FALSE,
  validation_split = 0,
  learning_rate = 0.0005,
  epochs = 40,
  verbose = 1,
  seed = seed,
  cos_sin_transform = TRUE
)

pred_batch_size <- 32
for (i in 1:ceiling(dim(test_data_coarse)[3] / pred_batch_size)) {
  batch_inds <- ((i - 1) * pred_batch_size + 1):min(i * pred_batch_size, dim(test_data_coarse)[3])
  cat("Predicting batch", i, "with indices", batch_inds, "\n")
  if (i == 1) {
    predictions_unet_time_cossin <- drop(predict(
      unet_model_time_cossin,
      test_data_coarse[, , batch_inds],
      time_points = test_times[batch_inds]
    ))
  } else {
    batch_preds <- drop(predict(
      unet_model_time_cossin,
      test_data_coarse[, , batch_inds],
      time_points = test_times[batch_inds]
    ))
    predictions_unet_time_cossin <- abind(predictions_unet_time_cossin, batch_preds, along = 3)
  }
}
save(predictions_unet_time_cossin, file = "predictions/predictions_unet_time_cossin_large_conv_40.RData")
load("predictions/predictions_unet_time_cossin_large_conv_40.RData")

compare_methods(test_data_coarse, predictions_unet_time_cossin, test_data_fine,
  time_idx = 1,
  title = "Ozone (µg/m³)", var_range = c(30, 120),
  method_name = "UNet"
)

unet_metrics_time_cossin <- calculate_metrics(
  predictions_unet_time_cossin,
  test_data_fine
)   
pred_errors_unet_time_cossin <- predictions_unet_time_cossin - test_data_fine

pred_errors_unet_agg <- apply(abs(pred_erros_unet), c(1,2), mean)
pred_errors_unet_agg <- array(pred_errors_unet_agg, dim = c(dim(pred_errors_unet_agg), 1))
colnames(pred_errors_unet_agg) <- colnames(pred_erros_unet)
rownames(pred_errors_unet_agg) <- rownames(pred_erros_unet)
pred_errors_unet_time_agg <- apply(abs(pred_errors_unet_time), c(1,2), mean)
pred_errors_unet_time_agg <- array(pred_errors_unet_time_agg, dim = c(dim(pred_errors_unet_time_agg), 1))
colnames(pred_errors_unet_time_agg) <- colnames(pred_errors_unet_time)
rownames(pred_errors_unet_time_agg) <- rownames(pred_errors_unet_time)
pred_errors_unet_time_cossin_agg <- apply(abs(pred_errors_unet_time_cossin), c(1,2), mean)
pred_errors_unet_time_cossin_agg <- array(pred_errors_unet_time_cossin_agg, dim = c(dim(pred_errors_unet_time_cossin_agg), 1))
colnames(pred_errors_unet_time_cossin_agg) <- colnames(pred_errors_unet_time_cossin)
rownames(pred_errors_unet_time_cossin_agg) <- rownames(pred_errors_unet_time_cossin)

compare_methods(pred_errors_unet_agg, pred_errors_unet_time_agg, pred_errors_unet_time_cossin_agg,
  time_idx = 1,
  title = "MAE over time", var_range = c(0, 3),
  method_names = c("UNet (baseline)", "UNet (RBF)", "UNet (sinusoidal)")
)
