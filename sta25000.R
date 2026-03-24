library(clock)
library(stringr)
library(stringi)
library(dplyr)
library(lubridate)
library(dlookr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(patchwork)
library(lwgeom)
library(fdasrvf)
library(future.apply)
library(future)
library(purrr)
library(tidyr)
library(splines)
library(refund)
library(janitor)
library(tidyverse)
sessionInfo()
pretty_names <- c(
  lon_start_z      = "Starting longitude",
  lat_start_z      = "Starting latitude",
  t_hours_total_z  = "Storm lifetime",
  seasonSpring     = "Spring",
  seasonSummer     = "Summer",
  seasonAutumn     = "Autumn",
  year_start_z     = "Year of formation"
)

source("C:\\Users\\Administrator\\Downloads\\trevor-harris-elasticdepth-8938a92\\R\\util.R")
source("C:\\Users\\Administrator\\Downloads\\trevor-harris-elasticdepth-8938a92\\R\\depths.R")
source("C:\\Users\\Administrator\\Downloads\\trevor-harris-elasticdepth-8938a92\\R\\outliers.R")


path  <- "C:\\Users\\Administrator\\Downloads\\1945-2023_english.csv"
df0 <- read_csv(path, show_col_types = FALSE) %>% clean_names()

parse_time_to_numeric <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "^\ufeff", "")
  x <- str_replace_all(x, "[\u200B-\u200D\uFEFF]", "")
  x <- str_replace_all(x, "[[:cntrl:]]", "")
  x <- str_trim(x)
  x <- str_replace(x, "(Z$|[+-]\\d\\d:\\d\\d$)", "")
  x <- str_replace(x, "\\.\\d+$", "")
  x <- str_replace(x, "T", " ")
  date_str <- str_sub(x, 1, 10)
  time_str <- str_sub(x, 12, 19)
  time_str[is.na(time_str) | time_str == "" | nchar(time_str) < 8] <- "00:00:00"
  d <- as.Date(date_str, format = "%Y-%m-%d")
  hh <- as.integer(str_sub(time_str, 1, 2))
  mm <- as.integer(str_sub(time_str, 4, 5))
  ss <- as.integer(str_sub(time_str, 7, 8))
  sec_of_day <- hh * 3600 + mm * 60 + ss
  num_sec <- as.numeric(d - as.Date("1970-01-01")) * 86400 + sec_of_day
  num_sec
}


df1 <- df0 %>%
  mutate(start_num = parse_time_to_numeric(typhoon_start_time),
         end_num = parse_time_to_numeric(typhoon_end_time),
         obs_num = parse_time_to_numeric(observation_time)) %>%
  select(-movement_direction, -wind_scale_level, -intensity)


df2 <- df1 %>%
  group_by(typhoon_id) %>%
  arrange(obs_num, .by_group = TRUE) %>%
  mutate(
    t0_num = min(obs_num, na.rm = TRUE),
    t_hours = (obs_num - t0_num) / 3600, 
    dur_hours = max(t_hours, na.rm = TRUE),
    t_scaled = if_else(dur_hours > 0, t_hours / dur_hours, 0),
    dt_hours = (obs_num - lag(obs_num)) / 3600) %>%
  ungroup() %>%
  select(-t0_num, -dur_hours) %>%
  relocate(t_hours, t_scaled, dt_hours, .after = observation_time)

df2 <- df2 %>%
  mutate(
    time_num = parse_time_to_numeric(observation_time),
    obs_time = as.POSIXct(time_num, origin = "1970-01-01", tz = "UTC"),
    year = as.integer(format(obs_time, "%Y")),
    month = as.integer(format(obs_time, "%m"))) %>%
  relocate(year, month, .after = observation_time)

df_pull <- df2 %>%
  group_by(typhoon_id) %>%
  summarise(n_obs = n(), span_h = max(t_hours, na.rm = TRUE),
            max_gap_h = max(dt_hours, na.rm = TRUE), .groups = "drop") %>%
  filter(n_obs >= 24, span_h >= 120, max_gap_h <= 6) %>% pull(typhoon_id)


df2 <- df2 %>% filter(typhoon_id %in% df_pull) %>% 
  mutate(observation_time = obs_time) %>% 
  select(-obs_time, -start_num, -obs_num, -time_num, -end_num, -dt_hours)


## df2 <- df2 %>% arrange(observation_time) %>% mutate(is_na = is.na(central_pressure))
## df2 <- df2 %>% mutate(block_id = cumsum(is_na != lag(is_na, default = FALSE)))
## block_info <- df2 %>% filter(is_na) %>% group_by(typhoon_id) %>%
##                       summarise(block_len = n(), start_time = min(observation_time), 
##                                 end_time = max(observation_time), .groups = "drop")

df2 %>% diagnose_numeric()

df2 <- df2 %>% filter(year >= 1949, central_pressure >= 800, !is.na(wind_speed), !is.na(central_pressure))

disp <- df2 %>% group_by(typhoon_id) %>%
  mutate(lon0 = longitude[which.min(observation_time)], lat0 = latitude [which.min(observation_time)],
         dlon = longitude - lon0, dlat = latitude  - lat0) %>%
  summarise(max_abs_lon = max(abs(dlon), na.rm = TRUE), max_abs_lat = max(abs(dlat), na.rm = TRUE),
            max_disp = pmax(max_abs_lon, max_abs_lat), .groups = "drop") %>%
  filter(max_disp >= 10)

keep_ids <- disp %>% filter(max_disp >= 10) %>% pull(typhoon_id)
df2 <- df2 %>% filter(typhoon_id %in% keep_ids)

wgs84_params <- function() {
  a <- 6378137.0 
  f <- 1 / 298.257223563
  e2 <- f * (2 - f)
  list(a = a, f = f, e2 = e2)
}
deg2rad <- function(x) x * pi / 180

llh_to_ecef_wgs84 <- function(lon_deg, lat_deg, h_m = 0, params = wgs84_params()) {
  stopifnot(length(lon_deg) == length(lat_deg))
  lon <- deg2rad(lon_deg)
  lat <- deg2rad(lat_deg)
  a <- params$a
  e2 <- params$e2
  sin_lat <- sin(lat)
  cos_lat <- cos(lat)
  cos_lon <- cos(lon)
  sin_lon <- sin(lon)
  N <- a / sqrt(1 - e2 * sin_lat^2)
  x <- (N + h_m) * cos_lat * cos_lon
  y <- (N + h_m) * cos_lat * sin_lon
  z <- (N * (1 - e2) + h_m) * sin_lat
  cbind(x_km = x / 1000, y_km = y / 1000, z_km = z / 1000)
}

ecef_diff_to_enu_wgs84 <- function(dx_km, dy_km, dz_km, lon0_deg, lat0_deg) {
  lon0 <- deg2rad(lon0_deg)
  lat0 <- deg2rad(lat0_deg)
  sin_lon0 <- sin(lon0); cos_lon0 <- cos(lon0)
  sin_lat0 <- sin(lat0); cos_lat0 <- cos(lat0)
  e <- -sin_lon0 * dx_km +  cos_lon0 * dy_km
  n <- -sin_lat0 * cos_lon0 * dx_km - sin_lat0 * sin_lon0 * dy_km + cos_lat0 * dz_km
  u <- cos_lat0 * cos_lon0 * dx_km + cos_lat0 * sin_lon0 * dy_km + sin_lat0 * dz_km
  cbind(e_km = e, n_km = n, u_km = u)
}


add_enu_from_first_wgs84 <- function(d, lon_col = "longitude", lat_col = "latitude", h_col = NULL) {
  lon <- as.numeric(d[[lon_col]])
  lat <- as.numeric(d[[lat_col]])
  ok <- is.finite(lon) & is.finite(lat)
  d <- d[ok, , drop = FALSE]
  lon0 <- lon[ok][1]
  lat0 <- lat[ok][1]
  h_m <- if (!is.null(h_col) && h_col %in% names(d)) as.numeric(d[[h_col]]) else 0
  p0 <- llh_to_ecef_wgs84(lon0, lat0, h_m = if (length(h_m) >= 1) h_m[1] else 0)
  p <- llh_to_ecef_wgs84(lon_deg = lon[ok], lat_deg = lat[ok], h_m = h_m)
  dp <- sweep(p, 2, p0, "-")
  enu <- ecef_diff_to_enu_wgs84(dp[,1], dp[,2], dp[,3], lon0, lat0)
  enu <- sweep(enu, 2, enu[1,], "-")
  d$x_km <- enu[, "e_km"]
  d$y_km <- enu[, "n_km"]
  d
}

df <- df2 %>%
  mutate(longitude = as.numeric(longitude), latitude  = as.numeric(latitude)) %>%
  arrange(typhoon_id, observation_time) %>%
  group_by(typhoon_id) %>%
  group_modify(\(d, key) add_enu_from_first_wgs84(d, lon_col="longitude", lat_col="latitude")) %>%
  relocate(x_km, y_km, .after = latitude) %>%
  ungroup()

df$typhoon_id <- as.factor(df$typhoon_id)


starts <- df %>% group_by(typhoon_id) %>% slice_min(observation_time, n = 1, with_ties = FALSE) %>% ungroup()
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
  st_make_valid() |> st_union() |> st_as_sf()



land_proj <- st_transform(land, 3857)
land_buffer <- st_buffer(land_proj, dist = 50000)
df_sf <- df |> st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) |> st_transform(3857)
df$is_land <- lengths(st_intersects(df_sf, land_buffer)) > 0
ids_primary <- starts %>% filter(longitude >= 120, longitude <= 165, latitude >= 5, latitude <=30) %>% pull(typhoon_id)
data <- df %>% filter(typhoon_id %in% ids_primary)




ends <- data %>% group_by(typhoon_id) %>% slice_max(observation_time, n = 1, with_ties = FALSE) %>% ungroup()
ends2 <- ends %>% mutate(nw_score = as.numeric(scale(latitude)) - as.numeric(scale(longitude)))
bad_ids <- ends2 %>% filter(nw_score >= sort(nw_score, decreasing = TRUE)[4] | longitude > 183 |
                              (latitude < 21 & longitude > 140) | (latitude < 25 & longitude < 95)) %>% pull(typhoon_id)

bad_tracks <- data %>% filter(typhoon_id %in% bad_ids)
good_tracks <- data %>% filter(!typhoon_id %in% bad_ids)


land_ll <- st_transform(land, 4326) %>% st_make_valid()
sf_use_s2(FALSE)
land_cut <- st_wrap_dateline(land_ll, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE)
land_cut <- st_collection_extract(land_cut, "POLYGON")
land_360 <- st_shift_longitude(land_cut)
sf_use_s2(TRUE)
land_360 <- st_make_valid(land_360)


p1 <- ggplot() + 
  geom_sf(data = land_360, linewidth = .2, fill = "grey90") +
  geom_point(data = ends, aes(longitude, latitude), size = 1.6, alpha = .7) +
  geom_point(data = ends %>% filter(typhoon_id %in% bad_ids), 
             aes(longitude, latitude), color = "red", size = 2) +
  coord_sf(xlim = c(80, 200), ylim = c(0, 65), expand = FALSE) +
  labs(title = "Typhoon Track End Locations", 
       x = "Longitude", y = "Latitude") +
  theme_minimal()

p2 <- ggplot() + 
  geom_sf(data = land_360, linewidth = 0.2, fill = "grey90") +
  geom_path(data = good_tracks, aes(longitude, latitude, group = typhoon_id), 
            color = "grey60", alpha = 0.25, linewidth = 0.3) +
  geom_path(data = bad_tracks, aes(longitude, latitude, group = typhoon_id, color = factor(typhoon_id)), linewidth = 0.5) +
  scale_color_viridis_d(guide = "none") +
  geom_point(data = ends %>% filter(typhoon_id %in% bad_ids), 
             aes(longitude, latitude), color = "red", size = 2) +
  coord_sf(xlim = c(80, 200), ylim = c(0, 65), expand = FALSE) +
  labs(title = "Outlying Typhoon Tracks", x = "Longitude",y = "Latitude") +
  theme_minimal()
p2
p1 | p2
data <- good_tracks



make_tmat_r2 <- function(dat, m = 101, id_col = "typhoon_id",
                         time_col = "t_scaled", x_col = "x_km", y_col = "y_km") {
  stopifnot(all(c(id_col, time_col, x_col, y_col) %in% names(dat)))
  grid <- seq(0, 1, length.out = m)
  ids <- sort(unique(dat[[id_col]]))
  N <- length(ids)
  tmat <- array(NA_real_, dim = c(2, m, N))
  for (i in seq_along(ids)) {
    d <- dat[dat[[id_col]] == ids[i], ]
    d <- d[order(d[[time_col]]), ]
    d <- d[!duplicated(d[[time_col]]), ] 
    tt <- as.numeric(d[[time_col]])
    xx <- as.numeric(d[[x_col]])
    yy <- as.numeric(d[[y_col]])
    tmat[1, , i] <- approx(tt, xx, xout = grid, rule = 2)$y
    tmat[2, , i] <- approx(tt, yy, xout = grid, rule = 2)$y
  }
  
  list(tmat = tmat, ids = ids, grid = grid)
}

obj <- make_tmat_r2(data, m = 101)
tmat <- obj$tmat
ids <- obj$ids



future::plan(future::multisession, workers = 8)  

dep <- depth.R2(tmat)  
outs <- elastic_outliers(dep, ka = 2, thresh = 0.8)


res <- dplyr::tibble(typhoon_id = ids, amp_depth = dep$amplitude, phs_depth = dep$phase,
                     amp_outlier = as.logical(outs$amp), phs_outlier = as.logical(outs$phs),
                     any_outlier = amp_outlier | phs_outlier)

res
p3 <- ggplot(res, aes(amp_depth, phs_depth)) +
  geom_point(color = "salmon", alpha = 0.4) +
  geom_point(data = dplyr::filter(res, any_outlier),
             color = "blue", size = 2) +
  geom_text(data = dplyr::filter(res, any_outlier),
            aes(label = typhoon_id),
            hjust = -0.05, size = 2.5) +
  labs(title = "Outlying Typhoon Tracks by Elastic Depth", x = "Amplitude", y = "Phase") +
  theme_minimal()

out_ids <- res %>% filter(any_outlier) %>% pull(typhoon_id)
bad_tracks  <- data %>% filter(typhoon_id %in% out_ids)
good_tracks <- data %>% filter(!typhoon_id %in% out_ids)
ends_all <- data %>% group_by(typhoon_id) %>% slice_max(observation_time, n = 1, with_ties = FALSE) %>% ungroup()

p4 <- ggplot() + 
  geom_sf(data = land_360, linewidth = 0.2, fill = "grey90") +
  geom_path(data = good_tracks, aes(longitude, latitude, group = typhoon_id),
            color = "grey70", alpha = 0.15, linewidth = 0.25) +
  geom_path(data = bad_tracks, aes(longitude, latitude, group = typhoon_id, color = typhoon_id),
            linewidth = 0.6) +
  scale_color_viridis_d(guide = "none") +
  geom_point(data = ends_all %>% filter(typhoon_id %in% out_ids), aes(longitude, latitude),
             color = "red", size = 2) +
  coord_sf(xlim = c(80, 180), ylim = c(0, 65), expand = FALSE) +
  labs(title = "Mapping of Outliers", x = "Longitude", y = "Latitude") +
  theme_minimal()

combined_plot <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Typhoon Track Outlier Analysis",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )
combined_plot
(p1 | p2) / (p3 | p4)
dt <- good_tracks
## cor(df$wind_speed, df$central_pressure, use = "complete.obs", method = "spearman")
df_pca <- dt %>% select(wind_speed, central_pressure)
pc <- prcomp(df_pca, center = TRUE, scale. = TRUE)
summary(pc)
dt$intensity <- pc$x[, 1]
pc$rotation
## cor(df$movement_speed, df$intensity, use = "complete.obs", method = "spearman")
dt <- dt %>% dplyr::select(-movement_speed)





m <- 101
grid <- seq(0, 1, length.out = m)

smooth_to_grid <- function(t, y, grid, smooth_df =48) {
  ok <- is.finite(t) & is.finite(y)
  t <- t[ok]
  y <- y[ok]
  if (length(t) < 6 || length(unique(t)) < 4) return(rep(NA_real_, length(grid)))
  fit <- lm(y ~ bs(t, df = smooth_df))
  as.numeric(predict(fit, newdata = data.frame(t = grid)))
}

dt_int_smooth <- dt %>%
  filter(is.finite(intensity), is.finite(t_scaled)) %>%
  group_by(typhoon_id) %>%
  summarise(
    intensity_smooth = list(smooth_to_grid(t_scaled, intensity, grid = grid, smooth_df = 8)),
    .groups = "drop"
  )

ok_dt <- dt_int_smooth %>%
  transmute(typhoon_id, ok_intensity = map_lgl(intensity_smooth, ~ any(is.finite(.x))))

table(ok_dt$ok_intensity)

good_ids_dt <- ok_dt %>% filter(ok_intensity) %>% pull(typhoon_id)

dt_int_smooth2 <- dt_int_smooth %>%
  filter(typhoon_id %in% good_ids_dt)

Y_wide_dt <- dt_int_smooth2 %>%
  transmute(typhoon_id, y = intensity_smooth) %>%
  unnest_longer(y, indices_to = "k") %>%
  mutate(t_reg = grid[k]) %>%
  select(-k) %>%
  pivot_wider(names_from = t_reg, values_from = y) %>%
  arrange(typhoon_id)

ids_dt <- Y_wide_dt$typhoon_id
Y_dt <- as.matrix(Y_wide_dt[, -1, drop = FALSE])




Ly <- lapply(seq_len(nrow(Y_dt)), function(i) as.numeric(Y_dt[i, ]))
Lt <- lapply(seq_len(nrow(Y_dt)), function(i) as.numeric(grid))
fp <- refund::fpca.sc(Y_dt, argvals = grid, npc = 8)


`%||%` <- function(a, b) if (!is.null(a)) a else b
mu <- fp$mu %||% fp$meanfd %||% fp$mu_hat
  
phi <- fp$efunctions %||% fp$phi %||% fp$eigenfunctions
lam <- fp$evalues %||% fp$lambda %||% fp$eigenvalues
scores <- fp$scores %||% fp$xiEst %||% fp$scores_hat
evalues <- fp$evalues
fve <- evalues / sum(evalues)

phi[,1] <- -phi[,1]
scores[,1]  <- -scores[,1]

phi[,2] <- -phi[,2]
scores[,2]  <- -scores[,2]
reconstruct_fpca_safe <- function(mu, phi, scores, K) {
  K_avail <- min(ncol(scores), ncol(phi))
  K_eff <- min(K, K_avail)
  mu_mat <- matrix(mu, nrow = nrow(scores), ncol = length(mu), byrow = TRUE)
  mu_mat + scores[, 1:K_eff, drop = FALSE] %*% t(phi[, 1:K_eff, drop = FALSE])
}

plot_fpca_like_demo_safe <- function(Y, grid, fp, K_use = 3, seed = 1, n_recon = 3) {
  stopifnot(is.matrix(Y), length(grid) == ncol(Y))
  parts <- get_fpca_parts(fp)
  mu_hat <- parts$mu
  phi_hat <- parts$phi
  lambda_hat <- parts$lam
  scores_hat <- parts$sc
  Yc <- sweep(Y, 2, mu_hat, "-")
  G_hat <- crossprod(Yc) / (nrow(Y) - 1)
  K_avail <- min(ncol(scores_hat), ncol(phi_hat))
  K_eff <- min(K_use, K_avail)
  X_hat <- reconstruct_fpca_safe(mu_hat, phi_hat, scores_hat, K = K_eff)
  R_hat <- Y - X_hat
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  
  matplot(grid, t(Y), type = "l", lty = 1, col = gray(0.7, 0.25),
          xlab = "t", ylab = "Y(t)", main = "Observed curves + mean")
  lines(grid, mu_hat, col = "red", lwd = 2)
  
  image(grid, grid, G_hat, xlab = "s", ylab = "t",
        main = "Estimated covariance surface G_hat(s,t)")
  contour(grid, grid, G_hat, add = TRUE)
  
 lam <- as.numeric(lambda_hat)
lam <- lam[is.finite(lam)]    
lam_pos <- pmax(lam, 0)
K_scree <- min(8, length(lam_pos))
k <- seq_len(K_scree)

pve  <- lam_pos / sum(lam_pos)
cpve <- cumsum(pve)

op <- par(no.readonly = TRUE)
on.exit(par(op), add = TRUE)
par(mar = c(5, 4, 4, 4) + 0.1)

plot(k, lam_pos[k], type = "b", pch = 1,
     xlab = "k", ylab = expression(hat(lambda)[k]),
     main = sprintf("Scree + CPVE (first %d)", K_scree),
     xaxt = "n")
axis(1, at = k, labels = k)

par(new = TRUE)
plot(k, cpve[k], type = "b", pch = 2,
     axes = FALSE, xlab = "", ylab = "",
     xlim = range(k), ylim = c(0, 1))
axis(4, at = seq(0, 1, 0.2))
mtext("Cumulative PVE", side = 4, line = 2.5)

abline(h = 0.9, lty = 3)
k90 <- which(cpve >= 0.9)[1]
if (!is.na(k90)) abline(v = k90, lty = 2)

legend("bottomright", bty = "n",
       legend = c(expression(hat(lambda)[k]), "CPVE", "90%",
                  if (!is.na(k90)) paste0("K90=", k90) else "K90=NA"),
       lty = c(1, 1, 3, 2),
       pch = c(1, 2, NA, NA))

  set.seed(12356)
  idx <- sample(seq_len(nrow(Y)), min(n_recon, nrow(Y)))
  matplot(grid, t(Y[idx, , drop = FALSE]), type = "l", lty = 1,
          col = gray(0.6, 0.6),
          xlab = "t", ylab = "value",
          main = paste0("Reconstruction (K=", K_eff, ")"))
  matlines(grid, t(X_hat[idx, , drop = FALSE]), col = "red", lwd = 2,
           lty = seq_along(idx))
  legend("topleft", bty = "n", lty = c(1, 1),  col = c(gray(0.6, 0.6), "red"),
         legend = c("Observed Y", "Reconstructed X_hat"))
}


png("fpca.png",  width = 16, height = 9, units = "in", res = 100, pointsize = 15)
plot_fpca_like_demo_safe(Y = Y_dt, grid = grid, fp = fp, K_use = 3)
dev.off()


get_fpca_parts <- function(fp) {
  mu   <- fp$mu %||% fp$meanfd %||% fp$mu_hat
  phi  <- fp$efunctions %||% fp$phi %||% fp$eigenfunctions
  lam  <- fp$evalues %||% fp$lambda %||% fp$eigenvalues
  sc   <- fp$scores %||% fp$xiEst %||% fp$scores_hat
  list(mu = as.numeric(mu), phi = as.matrix(phi), lam = as.numeric(lam), sc = as.matrix(sc))
}
plot_fpca_scores_pairs <- function(fp, pcs = c(1, 2, 3), col = "steelblue", pch = 19, cex = 0.7) {
  parts <- get_fpca_parts(fp)
  sc <- parts$sc
  stopifnot(max(pcs) <= ncol(sc))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  plot(sc[, pcs[1]], sc[, pcs[2]], xlab = paste0("PC", pcs[1]),
       ylab = paste0("PC", pcs[2]), main = paste0("PC", pcs[1], " vs PC", pcs[2]),
       pch = pch, col = col, cex = cex)
  abline(h = 0, v = 0, lty = 3, col = "gray50")
  
  plot(sc[, pcs[1]], sc[, pcs[3]], xlab = paste0("PC", pcs[1]),
       ylab = paste0("PC", pcs[3]), main = paste0("PC", pcs[1], " vs PC", pcs[3]),
       pch = pch, col = col, cex = cex)
  abline(h = 0, v = 0, lty = 3, col = "gray50")
  
  plot(sc[, pcs[2]], sc[, pcs[3]], xlab = paste0("PC", pcs[2]),
       ylab = paste0("PC", pcs[3]), main = paste0("PC", pcs[2], " vs PC", pcs[3]),
       pch = pch, col = col, cex = cex)
  abline(h = 0, v = 0, lty = 3, col = "gray50")
  invisible(sc)
}

plot_fpca_scores_pairs(fp)


plot_fpca_typhoon_components <- function(grid, phi, scores, fve, K = 3, ids = NULL, 
                                         max_lines = Inf, add_extreme_labels = TRUE) {
  stopifnot(nrow(phi) == length(grid))
  K <- min(K, ncol(phi), ncol(scores))
  n <- nrow(scores)
  if (is.null(ids)) ids <- seq_len(n)
  keep <- seq_len(n)
  if (is.finite(max_lines) && n > max_lines) {
    set.seed(1)
    keep <- sort(sample(keep, max_lines))
  }
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  par(mfrow = c(2,K), mar = c(4, 4, 3, 1))
  
  for (k in seq_len(K)) {
    y1 <- phi[, k]
    plot(grid, y1, type = "l", lwd = 2,
         xlab = "t (scaled)", ylab = paste0("Eigenfunction ", k),
         main = sprintf("Eigenfunction %d (FVE: %.2f%%)", k, 100 * fve[k]))
    abline(h = 0, col = "grey80", lwd = 2)
    contrib <- scores[, k] %o% phi[, k]
    ylim <- range(contrib[keep, , drop = FALSE], na.rm = TRUE)
    plot(NA, xlim = range(grid), ylim = ylim,
         xlab = "t (scaled)", ylab = paste0("Score * EF ", k),
         main = paste0("Typhoon trajectories on PC", k))
    abline(h = 0, col = "grey85", lwd = 2)
    apply(contrib[keep, , drop = FALSE], 1, function(v) {
      lines(grid, v, col = "grey40", lwd = 1)
    })
    if (add_extreme_labels) {
      kk <- keep
      i_pos <- kk[which.max(scores[kk, k])]
      i_neg <- kk[which.min(scores[kk, k])]
      lines(grid, contrib[i_pos, ], col = "red", lwd = 2)
      lines(grid, contrib[i_neg, ], col = "blue", lwd = 2)
    }
  }
}

plot_fpca_typhoon_components(grid = grid, phi = phi, scores = scores, fve = fve, K = 3, ids = ids_dt)


dt0 <- dt %>%
  filter(typhoon_id %in% ids_dt) %>%
  arrange(typhoon_id, t_hours)

month_to_season <- function(m) {
  dplyr::case_when(
    m %in% c(3, 4, 5)  ~ "Spring",
    m %in% c(6, 7, 8)  ~ "Summer",
    m %in% c(9,10,11)  ~ "Autumn",
    m %in% c(12,1,2)  ~ "Winter",
    TRUE ~ NA_character_
  )
}

X_dt <- dt0 %>%
  group_by(typhoon_id) %>%
  summarise(lon_start = first(longitude),
    lat_start = first(latitude),
    t_hours_total = max(t_hours, na.rm = TRUE) - 
      min(t_hours, na.rm = TRUE),
    year_start  = first(year),
    month_start = first(month),
    season = month_to_season(month_start),
    .groups = "drop") %>%
  mutate(season = factor(season, levels = c("Spring", "Summer", "Autumn", "Winter")))

scores <- fp$scores[, 1:3]
phi[,1] <- -phi[,1]
scores[,1]  <- -scores[,1]

phi[,2] <- -phi[,2]
scores[,2]  <- -scores[,2]
dat_reg <- bind_cols(X_dt %>% select(-typhoon_id, -month_start), as.data.frame(scores))%>%
  mutate(season = relevel(season, ref = "Winter"))
names(dat_reg)[(ncol(dat_reg)-2):ncol(dat_reg)] <- paste0("PC", 1:3)





K <-3
Ly <- lapply(seq_len(nrow(Y_dt)), function(i) as.numeric(Y_dt[i, ]))
Lt <- lapply(seq_len(nrow(Y_dt)), function(i) as.numeric(grid))


workGrid <- as.numeric(fp$workGrid)
dat_reg <- bind_cols(X_dt %>% select(-typhoon_id, -month_start), as.data.frame(scores))%>%
  mutate(season = relevel(season, ref = "Winter"))
names(dat_reg)[(ncol(dat_reg)-2):ncol(dat_reg)] <- paste0("PC", 1:3)

dat_reg_std <- dat_reg %>%
  mutate(
    lon_start_z = scale(lon_start)[,1],
    lat_start_z = scale(lat_start)[,1],
    t_hours_total_z = scale(t_hours_total)[,1],
    year_start_z = scale(year_start)[,1]
  )


fits <- lapply(1:K, function(k) {
  lm(as.formula(paste0("PC", k, " ~ season + lon_start_z + lat_start_z + ",
             "t_hours_total_z + year_start_z")),data = dat_reg_std)
})

X <- model.matrix(
  ~ lon_start_z + lat_start_z +
    t_hours_total_z + season + year_start_z,
  data = dat_reg_std
)[,-1]

cor(X, method = "spearman")

for (k in seq_along(fits)) {
  cat("\n============================\n")
  cat("Summary for PC", k, "\n")
  cat("============================\n")
  print(summary(fits[[k]]))
}


diagnose_lm_2plots <- function(fit, name = "") {
    plot(fitted(fit), resid(fit), main = paste(name, "Residuals vs Fitted"),
         xlab = "", ylab = "", pch = 19, col = rgb(0,0,0,0.4))
    abline(h = 0, lty = 2)
    qqnorm(resid(fit), main = paste(name, "Q-Q Plot"),
           xlab = "", ylab = "", pch = 19, col = rgb(0,0,0,0.4))
    qqline(resid(fit), col = "red")
  }

op <- par(no.readonly = TRUE)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

for (k in seq_along(fits)) {
  diagnose_lm_2plots(fits[[k]], paste0("PC", k))
}

par(op)


phi
phi[1, ] <- -phi[1, ] 
phi[2, ] <- -phi[2, ] 
xi <- scores
dat_sc <- bind_cols(dat_reg_std, as.data.frame(xi))

K_use <- 3
phi3  <- phi[, 1:K_use, drop = FALSE]

X <- model.matrix(~ lon_start_z + lat_start_z + t_hours_total_z + season + year_start_z,
                  data = dat_sc)
cnX <- colnames(X)

coefMat <- matrix(NA_real_, nrow = K_use, ncol = ncol(X),
                  dimnames = list(paste0("PC", 1:K_use), cnX))
for (k in 1:K_use) {
  ck <- coef(fits[[k]])
  coefMat[k, names(ck)] <- ck
}

alpha_hat <- coefMat[, "(Intercept)"]
Gamma_hat <- coefMat[, setdiff(cnX, "(Intercept)"), drop = FALSE]
mu <- as.numeric(mu)
mu_star <- as.numeric(mu + phi3 %*% alpha_hat)

beta_funcs <- lapply(colnames(Gamma_hat), function(nm) {
  as.numeric(phi3 %*% Gamma_hat[, nm])
})
names(beta_funcs) <- colnames(Gamma_hat)
coefMat
beta_funcs
K_use <- length(fits)
pnames <- setdiff(colnames(X), "(Intercept)")
vcov_list <- lapply(fits, vcov)
beta_ci <- list()
for (nm in pnames) {
  gamma_hat <- sapply(fits, function(fit) coef(fit)[nm])
  var_gamma <- diag(sapply(1:K_use, function(k) {
    vcov_list[[k]][nm, nm]
  }))
  var_beta_t <- apply(phi[,1:K_use, drop=FALSE], 1, function(phi_t) {
    t(phi_t) %*% var_gamma %*% phi_t
  })
  se_beta_t <- sqrt(var_beta_t)
  beta_hat_t <- as.numeric(phi[,1:K_use, drop=FALSE] %*% gamma_hat)
  beta_ci[[nm]] <- list(
    beta = beta_hat_t,
    lower = beta_hat_t - 1.96 * se_beta_t,
    upper = beta_hat_t + 1.96 * se_beta_t
  )
}

workGrid <- grid 

par(mfrow = c(3,3), mar = c(4,4,3,1))

for (nm in names(beta_ci)) {
  obj <- beta_ci[[nm]]
  title_use <- pretty_names[nm]
  if (is.na(title_use)) title_use <- nm 
  plot(workGrid, obj$beta, type="l", lwd=2, ylim=range(c(obj$lower, obj$upper)),
       xlab="t", ylab=expression(beta(t)), main = title_use)
  lines(workGrid, obj$lower, lty=2)
  lines(workGrid, obj$upper, lty=2)
  
  abline(h=0, lty=3, col="gray50")
}





m <- 101
grid <- seq(0, 1, length.out = m)

smooth_to_grid <- function(t, y, grid, smooth_df = 8) {
  ok <- is.finite(t) & is.finite(y)
  t <- t[ok]; y <- y[ok]
  if (length(t) < 6 || length(unique(t)) < 4) return(rep(NA_real_, length(grid)))
  fit <- lm(y ~ bs(t, df = smooth_df))
  as.numeric(predict(fit, newdata = data.frame(t = grid)))
}
int_smooth <- smooth_to_grid(
  t = dt$t_scaled,
  y = dt$intensity,
  grid = grid,
  smooth_df = 6
)


dt_lon_smooth <- dt %>%
  filter(typhoon_id %in% ids_dt, is.finite(x_km), is.finite(t_scaled)) %>%
  group_by(typhoon_id) %>%
  summarise(
    lon_smooth = list(smooth_to_grid(t_scaled, x_km, grid = grid, smooth_df = 8)),
    .groups = "drop"
  )

dt_lat_smooth <- dt %>%
  filter(typhoon_id %in% ids_dt, is.finite(y_km), is.finite(t_scaled)) %>%
  group_by(typhoon_id) %>%
  summarise(
    lat_smooth = list(smooth_to_grid(t_scaled, y_km, grid = grid, smooth_df = 8)),
    .groups = "drop"
  )

Xlon_wide <- dt_lon_smooth %>%
  transmute(typhoon_id, x = lon_smooth) %>%
  unnest_longer(x, indices_to = "k") %>%
  mutate(t_reg = grid[k]) %>%
  select(-k) %>%
  pivot_wider(names_from = t_reg, values_from = x)

Xlat_wide <- dt_lat_smooth %>%
  transmute(typhoon_id, x = lat_smooth) %>%
  unnest_longer(x, indices_to = "k") %>%
  mutate(t_reg = grid[k]) %>%
  select(-k) %>%
  pivot_wider(names_from = t_reg, values_from = x)

Xlon_wide2 <- Xlon_wide %>% slice(match(ids_dt, typhoon_id))
Xlat_wide2 <- Xlat_wide %>% slice(match(ids_dt, typhoon_id))
lon_dt <- as.matrix(Xlon_wide2[, -1, drop = FALSE])
lat_dt <- as.matrix(Xlat_wide2[, -1, drop = FALSE])

drop_start <- 3
drop_end <- 0
m <- 101
grid <- seq(0, 1, length.out = m)
keep <- seq(drop_start + 1, m - drop_end)
grid  <- grid[keep]
lon_dt <- lon_dt[, keep, drop = FALSE]
lat_dt <- lat_dt[, keep, drop = FALSE]
Y_dt <- Y_dt[, keep, drop = FALSE]

## --- (B) standardize predictors (overwrite; column-wise over storms) ---
lon_dt <- scale(lon_dt, center = TRUE, scale = TRUE)
lat_dt <- scale(lat_dt, center = TRUE, scale = TRUE)
lon_dt <- as.matrix(lon_dt)
lat_dt <- as.matrix(lat_dt)


dat <- list(Xlon = lon_dt, Xlat = lat_dt, Y = Y_dt)
res_pt <- ptFCReg(tGrid = grid, dat = dat)

bw <- 2.5 / (length(grid) - 1)
res_sm <- smPtFCRegCoef(res_pt, bw = bw, kernel_type = "epan")

beta0_hat <- res_sm$beta0
beta_lon  <- res_sm$beta[1, ]
beta_lat  <- res_sm$beta[2, ]

fit_once <- function(idx){
  dat_boot <- list(
    Xlon = lon_dt[idx, ],
    Xlat = lat_dt[idx, ],
    Y = Y_dt[idx, ]
  )
  res_pt <- ptFCReg(tGrid = grid, dat = dat_boot)
  bw <- 2.5 / (length(grid) - 1)
  res_sm <- smPtFCRegCoef(res_pt, bw = bw, kernel_type = "epan")
  list(
    beta_lon = res_sm$beta[1, ],
    beta_lat = res_sm$beta[2, ]
  )
}

set.seed(123)
B <- 300
n <- nrow(Y_dt)

boot_lon <- matrix(NA, B, length(grid))
boot_lat <- matrix(NA, B, length(grid))

for(b in 1:B){
  idx <- sample(1:n, n, replace = TRUE)
  fit_b <- fit_once(idx)
  boot_lon[b, ] <- fit_b$beta_lon
  boot_lat[b, ] <- fit_b$beta_lat
}

lon_se <- apply(boot_lon, 2, sd)
lat_se <- apply(boot_lat, 2, sd)

z <- qnorm(0.975)

lon_upper <- beta_lon + z * lon_se
lon_lower <- beta_lon - z * lon_se

lat_upper <- beta_lat + z * lat_se
lat_lower <- beta_lat - z * lat_se

den <- beta_lon^2 + beta_lat^2
p_lon <- ifelse(den > 0, beta_lon^2 / den, NA_real_)
p_lat <- ifelse(den > 0, beta_lat^2 / den, NA_real_)


par(mfrow=c(2,2), mar=c(4.5,4.5,2,1))

plot(grid, beta_lon, type="l", lwd=2,
     ylim=range(lon_lower, lon_upper),
     xlab="t", ylab=expression(beta[lon](t)),
     main="Longitude effect")
lines(grid, lon_upper, lty=2)
lines(grid, lon_lower, lty=2)
abline(h=0, lty=3)

plot(grid, beta_lat, type="l", lwd=2, ylim=range(lat_lower, lat_upper),
     xlab="t", ylab=expression(beta[lat](t)), main="Latitude effect")
lines(grid, lat_upper, lty=2)
lines(grid, lat_lower, lty=2)
abline(h=0, lty=3)

cols <- colorRampPalette(c("blue","red"))(length(grid))
plot(beta_lon, beta_lat, col=cols, pch=16,
     xlab=expression(beta[x](t)),
     ylab=expression(beta[y](t)),
     main="Effect path")
points(beta_lon[1], beta_lat[1], pch=19, col="blue", cex=1.5)
points(beta_lon[length(keep)], beta_lat[length(keep)], pch=19, col="red", cex=1.5)
legend("topright", legend=c("Start (t=0)", "End (t=1)"), col=c("blue","red"), pch=19, bty="n")

plot(grid, p_lon, type="l", lwd=2,ylim=c(0,1),
     xlab="t",
     ylab="Relative contribution",
     main="Longitude vs Latitude dominance")
polygon(c(grid, rev(grid)), c(p_lon, rep(0, length(p_lon))), col=rgb(0,0,1,0.3), border=NA) 
polygon(c(grid, rev(grid)), c(rep(1, length(p_lon)), rev(p_lon)), col=rgb(1,0,0,0.3), border=NA) 
lines(grid, p_lon, lwd=2) 
abline(h=0.5, lty=2) 
legend("bottomright", legend=c("Longitude", "Latitude"), fill=c(rgb(0,0,1,0.3), rgb(1,0,0,0.3)), bty="n", cex=0.9)
abline(h=0.5, lty=2)

par(mfrow=c(1,1))





Y_hat_c <- beta0_hat +
  lon_dt * matrix(beta_lon, nrow=n, ncol=length(grid), byrow=TRUE) +
  lat_dt * matrix(beta_lat, nrow=n, ncol=length(grid), byrow=TRUE)

R_c <- Y  - Y_hat_c

par(mfrow=c(1,3), mar=c(3,4,2,1))
matplot(grid, t(R_c), type="l", lty=1, col=rgb(0,0,0,0.15),
        xlab="t", ylab="r_i(t)", main="Residual curves (Y centered)")
plot(grid, colMeans(R_c), type="l", lwd=2,
     xlab="t", ylab=expression(bar(r)(t)), main="Mean residual"); abline(h=0,lty=2)
plot(grid, apply(R_c,2,var), type="l", lwd=2,
     xlab="t", ylab="Var(r(t))", main="Residual variance")
par(mfrow=c(1,1))











ids5 <- dt %>%
  filter(is.finite(t_scaled), is.finite(x_km), is.finite(y_km), is.finite(intensity)) %>%
  group_by(typhoon_id) %>%
  summarise(n  = n(), ut = n_distinct(t_scaled), peak = max(intensity, na.rm = TRUE), 
    .groups = "drop") %>% filter(n >= 6, ut >= 4) %>%
  mutate(bin = dplyr::ntile(peak, 5)) %>%
  group_by(bin) %>% slice_sample(n = 1) %>% ungroup() %>%
  pull(typhoon_id)

plot_overlay_raw_smooth <- function(data, ycol, title){
  
  raw <- data %>%
    filter(typhoon_id %in% ids5,
           is.finite(t_scaled), is.finite(.data[[ycol]])) %>%
    select(typhoon_id, t_scaled, y = all_of(ycol)) %>%
    arrange(typhoon_id, t_scaled)
  
  sm <- raw %>%
    group_by(typhoon_id) %>%
    summarise(y_smooth = list(smooth_to_grid(t_scaled, y, grid=grid, smooth_df=8)), .groups="drop") %>%
    unnest_longer(y_smooth) %>%
    group_by(typhoon_id) %>%
    mutate(t_grid = grid[row_number()]) %>%
    ungroup()
  
  ggplot() +
    geom_line(data = raw, aes(t_scaled, y, group = typhoon_id), color = "grey40", alpha = 0.25, linewidth = 0.4) +
    geom_point(data = raw, aes(t_scaled, y, group = typhoon_id),
      color = "grey10", alpha = 0.35, size = 0.8) +
    geom_line(data = sm, aes(t_grid, y_smooth, group = typhoon_id, color = typhoon_id),
      linewidth = 1.2, show.legend = FALSE) +
    labs(x = "Scaled lifetime (t)", y = title) +
    theme_bw() +
    theme(legend.position = "none", plot.title = element_blank())
}
 
p_lon <- plot_overlay_raw_smooth(dt, "x_km", "Longitude")
p_lat <- plot_overlay_raw_smooth(dt, "y_km", "Latitude")
p_int <- plot_overlay_raw_smooth(dt, "intensity", "Intensity")
(p_lon | p_lat | p_int)

(p_lon | p_lat | p_int) +
  plot_annotation(
    title = "Raw vs Smoothed Functional Trajectories (5 Sample Typhoons)",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  )




check_diff_one <- function(dt, id, ycol, grid, smooth_df=8){
  tmp <- dt %>% filter(typhoon_id==id, is.finite(t_scaled), is.finite(.data[[ycol]])) %>%
    arrange(t_scaled)
  s <- smooth_to_grid(tmp$t_scaled, tmp[[ycol]], grid=grid, smooth_df=smooth_df)
  raw_on_grid <- approx(tmp$t_scaled, tmp[[ycol]], xout=grid, rule=2)$y
  c(max_abs = max(abs(raw_on_grid - s), na.rm=TRUE),
    mean_abs = mean(abs(raw_on_grid - s), na.rm=TRUE))
}

sapply(ids_dt[1:5], \(id) check_diff_one(dt, id, "x_km", grid))