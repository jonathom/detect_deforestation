Introduction
============

Methods
=======

The time series Landsat 8 satellite images used for this research were
already provided in a state where cloud cover were already removed. The
data covers the period 01-01-2013 to 31-12-2019. The Landsat 8 data
covered row 001, paths 066 and 067.

The time series data was aggregated to a monthly and quarterly NDVI
(Normalized Difference Vegetation Index) data where the median pixel
values per row and column were selected. This resulted to 12 NDVI images
per year in the former approach and 4 NDVI images per year in the latter
approach

Prediction of missing values in satellite data are carried out using
gapfill package in R . The gapfill approach was designed to carry out
predictions on satellite data that were recorded at equally spaced
points of time. Based on Gerber et.al 2016 , they applied the algorithm
to MODIS NDVI data with cloud cover scenarios of 50% missing data. The
method was further compared to Gapfill-Python and TIMESAT and it
provided the most accurate prediction in terms of RMSE.

Gapfill was appealing to this research because it’s capable of handling
large spatio-temporal data, it’s user friendly and capability tailored
to specific features of different satellite images. The predictions of
the missing values are based on a subset-predict procedure, i.e., each
missing value is predicted separately by (1) selecting subsets of the
data that are in a neighborhood around the missing point and (2)
predicting the missing value based on the subset(Gerber, 2016).

In this research we also explored to tailor gapfill by customizing the
iMax parameter which is the maximum number of iterations until NA is
returned as predicted value (Gerber, 2016). The research compares the
iMax parameter at value 5 (five) for five iterations versus undefined
which results to infinite iterations by default.

BFAST(Breaks For Additive Season and Trend)

Input Data and Preparations
---------------------------

-   Package description
-   PRODES data sorted by years can be found here: [PRODES yearly
    deforestation](http://terrabrasilis.dpi.inpe.br/download/dataset/legal-amz-prodes/vector/yearly_deforestation.zip)

As seen in “final.Rmd”. The subdirectory `L8cube_subregion` contains a
NDVI time series as single `.tif` files, a file per acquisition, as
input data.

``` r
library(stars)
library(gapfill)
library(bfast)
library(zoo)
library(raster)
subdir = "landsat_monthly"
f = paste0(subdir, "/", list.files(subdir))
st = merge(read_stars(f)) # make stars object
plot(st)
subdir = "landsat_quarterly"
f = paste0(subdir, "/", list.files(subdir))
st_q = merge(read_stars(f)) # make stars object
plot(st_q)
```

<img src="main_files/figure-markdown_github/load-data-1.png" width="50%" /><img src="main_files/figure-markdown_github/load-data-2.png" width="50%" />

``` r
# load PRODES data
prod <- read_sf("./yearly_deforestation/yearly_deforestation.shp")
prod_3857 <- st_make_valid(st_transform(prod, crs = st_crs(st)))
prod_crop <- st_crop(prod_3857, st) # clip
write_sf(prod_crop, "./yearly_deforestation/PRODES_cropped.shp", overwrite = TRUE)
```

``` r
prod <- read_sf("./yearly_deforestation/PRODES_cropped.shp")
# prod <- prod[prod$YEAR < 2019,]
```

``` r
# whats in it?
# plot(st)
plot(prod["YEAR"], axes = TRUE, main = "PRODES Deforestation Data Colored by Year")
```

![](main_files/figure-markdown_github/plot-AOI-1.png)

Prepare for Gapfill
-------------------

`Gapfill` documentation tells us that as input, a 4-dimensional numeric
array is needed, with dimensions x, y, seasonal index (doy), year.

``` r
prep_gapfill <- function(st, doy, ts) {
  # st is stars object, doy is day of year vector, ts is number of timesteps per year
  
  # get pixels of whole dataset
  imgdata <- c(st[,,,][[1]])

  # make labels
  xlab <- seq(from = attr(st, "dimensions")[[1]]$offset, by = attr(st, "dimensions")[[1]]$delta, length.out = attr(st, "dimensions")[[1]]$to)
  ylab <- seq(from = attr(st, "dimensions")[[2]]$offset, by = attr(st, "dimensions")[[2]]$delta, length.out = attr(st, "dimensions")[[2]]$to)
  years <- seq(2013,2019,1)

  # make array, transpose
  h <- array(imgdata, dim = c(140, 140, ts, 7), dimnames = list(xlab, ylab, doy, years))
  # x, y is switched between stars and these arrays
  h <- aperm(h, c(2,1,3,4))
  return(h)
}

doy_12 <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
doy_4 <- c(1, 91, 182, 274)

ma_monthly <- prep_gapfill(st, doy_12, 12)
ma_quarter <- prep_gapfill(st_q, doy_4, 4)
```

Gapfill
-------

``` r
d <- Gapfill(ma_monthly, iMax = 5)
saveRDS(d, "./monthly_iMax5_140_gapfilled.rds")
e <- Gapfill(ma_quarter, iMax = 5)
saveRDS(e, "./quarterly_iMax5_140_gapfilled.rds")
f <- Gapfill(ma_quarter) # iMax defaults to infinite
saveRDS(f, "./quarterly_iMaxInf_140_gapfilled.rds")
```

### Gapfill Results

``` r
gf_monthly <- readRDS("monthly_iMax5_140_gapfilled.rds")
Image(gf_monthly$fill, zlim = c(0.2, 1)) + ggtitle("Gapfilled Monthly Data")
gf_quarterly <- readRDS("quarterly_iMax5_140_gapfilled.rds")
Image(gf_quarterly$fill, zlim = c(0.2, 1)) + ggtitle("Gapfilled Quarterly Data")
gf_quarterly_inf <- readRDS("quarterly_iMaxInf_140_gapfilled.rds")
Image(gf_quarterly_inf$fill, zlim = c(0.2, 1)) + ggtitle("Gapfilled Quarterly Data with iMax = infinite")
```

<img src="main_files/figure-markdown_github/load-gapfill-1.png" width="33%" /><img src="main_files/figure-markdown_github/load-gapfill-2.png" width="33%" /><img src="main_files/figure-markdown_github/load-gapfill-3.png" width="33%" />

### Gapfill Results - Closeup

Here, October to December of 2013 are plotted for comparison. First, The
input data is plotted. Below that, the gapfilled dataset are plotted.

``` r
Image(ma_monthly[,,10:12,1], zlim = c(0.2, 1)) + ggtitle("Monthly Input Data, Oct - Dec 2013")
Image(ma_quarter[,,4,1], zlim = c(0.2, 1)) + ggtitle("Quarterly Input Data, Last Quarter 2013")
```

<img src="main_files/figure-markdown_github/zoom-gapfill-input-1.png" width="50%" /><img src="main_files/figure-markdown_github/zoom-gapfill-input-2.png" width="50%" />

``` r
Image(gf_monthly$fill[,,10:12,1], zlim = c(0.2, 1)) + ggtitle("Monthly Gapfilled Data, Oct - Dec 2013, iMax = 5")
Image(gf_quarterly$fill[,,4,1], zlim = c(0.2, 1)) + ggtitle("Quarterly Gapfilled Data, Last Quarter 2013, iMax = 5")
Image(gf_quarterly_inf$fill[,,4,1], zlim = c(0.2, 1)) + ggtitle("Quarterly Gapfilled Data, Last Quarter 2013, iMax = infinite")
```

<img src="main_files/figure-markdown_github/zoom_gapfill-1.png" width="33%" /><img src="main_files/figure-markdown_github/zoom_gapfill-2.png" width="33%" /><img src="main_files/figure-markdown_github/zoom_gapfill-3.png" width="33%" />

Calculate BFAST Over AOI
------------------------

``` r
bfast_on_tile <- function(gapfill_matrix, by, ts, order) {
  dims <- dim(gapfill_matrix)
  result <- matrix(rep(FALSE, dims[1]*dims[2]), ncol = dims[1])
  for (i in 1:dims[1]) { # looping through x
    for (j in 1:dims[2]) { # loops through y
      raw_px_ts <- as.vector(gapfill_matrix[i,j,,]) # create pixel timeseries vector
      px_ts_obj <- as.ts(zoo(raw_px_ts, seq(2013, by = by, length.out = ts))) # make into ts object
      bfm_obj <- bfastmonitor(px_ts_obj, start = 2019, order = order) # bfastmonitor of pixel timeseries
      brkpoint <- bfm_obj$breakpoint
      if(!is.na(brkpoint)) { # if breakpoint is available..
        result[i,j] <- TRUE # .. write TRUE to solution raster
      } else {
        # FALSE
      }
    }
  }
  return(result)
}

bfast_monthly3 <- bfast_on_tile(gf_monthly$fill, by = .08333333, ts = 84, order = 3)
bfast_monthly2 <- bfast_on_tile(gf_monthly$fill, by = .08333333, ts = 84, order = 2)
bfast_quarter2 <- bfast_on_tile(gf_quarterly$fill, by = 0.25, ts = 28, order = 2)
# order = 2 was chosen because order 3 doesn't work on our quarterly aggreggated data
saveRDS(bfast_monthly2, "bfast_monthly2.rds")
saveRDS(bfast_quarter2, "bfast_quarter2.rds")
# warning: too few observations in history period
```

``` r
bfast_monthly <- readRDS("bfast_monthly2.rds")
bfast_quarter <- readRDS("bfast_quarter2.rds")
```

To eliminate errors that may appear in already deforested areas, these
arreas are simply excluded, according to PRODES reference data.

``` r
# rasterize PRODES data
# 2019 = TRUE, !2019 = FALSE
ras <- rasterize(prod, as(st[,,,5], "Raster"), "YEAR")
prodes <- aperm(matrix(ras[], ncol = 140), c(2,1))
prodes[prodes < 2019] <- FALSE
prodes[prodes == 2019] <- TRUE
prodes[is.na(prodes)] <- FALSE

# to mask out previous deforestation
# <2019 = TRUE, !<2019 = FALSE
prodes_prev <- aperm(matrix(ras[], ncol = 140), c(2,1))
prodes_prev[prodes_prev < 2019] <- TRUE
prodes_prev[prodes_prev == 2019] <- FALSE
prodes_prev[is.na(prodes_prev)] <- FALSE

bfast_monthly[prodes_prev == 1] <- FALSE

bfast_quarter[prodes_prev == 1] <- FALSE
```

``` r
table1 <- addmargins(table(bfast_monthly, prodes))
table2 <- addmargins(table(bfast_quarter, prodes))

accuracies <- function(table1) {
  # overall accuracy
  P0 <- (table1[1] + table1[5]) / table1[9]
  # producer's accuracy, Probability of classifying a pixel correctly
  pa_f <- table1[1] / table1[3] # FALSE
  pa_t <- table1[5] / table1[6] # TRUE
  # user's accuracy, Probability of a pixel being the classified type
  ua_f <- table1[1] / table1[7] # FALSE
  ua_t <- table1[5] / table1[8] # TRUE
  # kappa
  # chance that both TRUE / FALSE randomly
  tr <- (table1[8] / table1[9]) * (table1[6] / table1[9])
  fr <- (table1[7] / table1[9]) * (table1[3] / table1[9])
  Pe <- tr + fr
  kappa <- (P0 - Pe) / (1 - Pe)
  
  return(list("Overall Accuracy" = P0*100, "Prod. Acc. FALSE" = pa_f*100, "Prod. Acc. TRUE" = pa_t*100, "User's Acc. FALSE" = ua_f*100, "User's Acc. TRUE" = ua_t*100, "Kappa" = kappa))
}
```

Results
=======

``` r
Image(bfast_monthly) + ggtitle("Monthly Data") + theme(plot.title = element_text(size=22))
Image(bfast_quarter) + ggtitle("Quarterly Data") + theme(plot.title = element_text(size=22))
Image(prodes) + ggtitle("PRODES Data") + theme(plot.title = element_text(size=22))

table1
```

    ##              prodes
    ## bfast_monthly     0     1   Sum
    ##         FALSE 15905   170 16075
    ##         TRUE   1428  2097  3525
    ##         Sum   17333  2267 19600

``` r
table2
```

    ##              prodes
    ## bfast_quarter     0     1   Sum
    ##         FALSE 15720   374 16094
    ##         TRUE   1613  1893  3506
    ##         Sum   17333  2267 19600

``` r
array(c(accuracies(table1), accuracies(table2)), dim = c(6,2), dimnames = list(c("Overall Accuracy", "Prod. Acc. FALSE", "Prod. Acc. TRUE", "User's Acc. FALSE", "User's Acc. TRUE", "Kappa"), c("monthly", "quarterly")))
```

    ##                   monthly   quarterly
    ## Overall Accuracy  91.84694  89.86224 
    ## Prod. Acc. FALSE  91.76138  90.69405 
    ## Prod. Acc. TRUE   92.5011   83.50243 
    ## User's Acc. FALSE 98.94246  97.67615 
    ## User's Acc. TRUE  59.48936  53.99315 
    ## Kappa             0.6788956 0.5995541

<img src="main_files/figure-markdown_github/results-1.png" width="33%" /><img src="main_files/figure-markdown_github/results-2.png" width="33%" /><img src="main_files/figure-markdown_github/results-3.png" width="33%" />

Discussion
==========

Conclusion
==========

References
==========

Gerber F, Furrer R, Schaepman-Strub G, de Jong R, Schaepman ME (2016)
Predicting missing values in spatio-temporal satellite data.

Verbesselt, J., Hyndman, R., Newnham, G., & Culvenor, D. (2010).
Detecting trend and seasonal changes in satellite image time series.

Verbesselt, J., Zeileis, A., & Herold, M. (2013). Near real-time
disturbance detection using satellite image time series.
