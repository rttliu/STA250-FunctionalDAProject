# Functional Representation and Regression of Typhoon Evolution in Western North Pacific Area

## Overview
This repository contains the data processing pipelines, functional data analysis (FDA) models, and regression frameworks used to investigate typhoon intensity trajectories in the Western North Pacific area. 

By analyzing 1,101 historical storms, this project represents each storm as a smooth intensity curve over its normalized lifetime and applies advanced statistical modeling to uncover stable empirical patterns. The project contrasts with computationally expensive physics-based models by leveraging statistical methods in an $L^2$ space to preserve and analyze the full life-cycle shape of storm evolution.

## Questions Addressed
1. What are the dominant modes of variation in the evolution of typhoon intensity?
2. How do climatic and environmental factors (like origin location, season, etc.) shift these evolution patterns?
3. How does the spatial-intensity relationship (influence of eastward and northward displacement) evolve during the life cycle of a typhoon?

## Data Source
* **Dataset:** International Best Track Archive for Climate Stewardship (IBTrACS), Version 4.
* **Provider:** NOAA National Centers for Environmental Information (NCEI).
* **Features Used:** Time-stamped longitude, latitude, wind speed, and central pressure. Observations were filtered to include records after 1949 and standardized to UTC.

## Methodology
The project employs a comprehensive Functional Data Analysis (FDA) framework:

1. **Data Preprocessing & Smoothing:** * A 1D intensity index was constructed using Principal Component Analysis (PCA) on wind speed and central pressure (explaining 95.9% of the variance).
   * The longitudinal, latitudinal, and intensity components were smoothed using penalized B-splines and evaluated on a common uniform grid.
   * An elastic depth approach (square-root velocity framework) was utilized to identify and remove functional outliers.

2. **Functional Principal Component Analysis (FPCA):**
   * FPCA decomposed the trajectory variation into orthogonal modes around the population mean.
   * Over 90% of the variance is explained by just three dominant modes: **overall intensity magnitude**, **timing shift of the peak**, and **high-frequency oscillation**.

3. **Function-on-Scalar Regression**
This method examines how single-value storm characteristics (like starting location, season, and lifetime) shape the *entire* continuous intensity curve, rather than just a single peak value.
* Overall Magnitude: A storm's starting location and season heavily dictate its lifetime strength (e.g., storms forming further south/east are consistently stronger).
* Peak Timing: Seasonality is the strongest driver of *when* a storm reaches its maximum intensity.
* Fluctuations: Starting longitude and formation year influence rapid, short-term intensity changes.

4. **Concurrent Functional Regression**
This approach models the instantaneous, moment-by-moment relationship between a storm's physical movement and its intensity at that exact same time. 
* Early Stage: Spatial displacement has little to no measurable effect on intensity.
* Mid-Life (Development): Moving northward strongly correlates with intensification, while moving eastward is linked to weakening.
* Late Stage (Decay): As the storm recurves and dies down, further displacement in *either* direction (north or east) is associated with continued weakening.
