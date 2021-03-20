-   Observe structure: Intro, Methods, Results, Discussion, Conclusion

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
subdir = "landsat_quarterly"
f = paste0(subdir, "/", list.files(subdir))
st = merge(read_stars(f)) # make stars object
plot(st)
```

![](main_files/figure-markdown_github/load-data-1.png)

``` r
# load PRODES data
# prod <- read_sf("./yearly_deforestation/yearly_deforestation.shp")
# prod_3857 <- st_make_valid(st_transform(prod, crs = st_crs(st)))
# prod_crop <- st_crop(prod_3857, st) # clip
# write_sf(prod_crop, "./yearly_deforestation/PRODES_cropped.shp", overwrite = TRUE)
prod <- read_sf("./yearly_deforestation/PRODES_cropped.shp")
# prod <- prod[prod$YEAR < 2019,]
```

``` r
# whats in it?
# plot(st)
plot(as(st[,,,6], "Raster"))
plot(prod["YEAR"], add = TRUE)
```

![](main_files/figure-markdown_github/plot-AOI-1.png)

``` r
plot(prod["YEAR"])
```

![](main_files/figure-markdown_github/plot-AOI-2.png)

Prepare for Gapfill
-------------------

`Gapfill` documentation tells us that as input, a 4-dimensional numeric
array is needed, with dimensions x, y, seasonal index (doy), year.

``` r
# get pixels of whole dataset
imgdata <- c(st[,,,][[1]])

# make labels
xlab <- seq(from = attr(st, "dimensions")[[1]]$offset, by = attr(st, "dimensions")[[1]]$delta, length.out = attr(st, "dimensions")[[1]]$to)
ylab <- seq(from = attr(st, "dimensions")[[2]]$offset, by = attr(st, "dimensions")[[2]]$delta, length.out = attr(st, "dimensions")[[2]]$to)
doi <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
doi <- c(1, 91, 182, 274)
years <- seq(2013,2019,1)

# make array, transpose
h <- array(imgdata, dim = c(140, 140, 4, 7), dimnames = list(xlab, ylab, doi, years))
h <- aperm(h, c(2,1,3,4))

# all
Image(h[,,2,3:4])
```

![](main_files/figure-markdown_github/create-gapfill-input-1.png)

``` r
# original stars
plot(st[,,,c(5, 17)])
```

![](main_files/figure-markdown_github/create-gapfill-input-2.png)

Gapfill
-------

``` r
# d <- Gapfill(h, iMax = 5)
# saveRDS(d, "./iMax5_140_gapfilled_quarterly.rds")
gf_monthly <- readRDS("monthly_iMax5_140_gapfilled.rds")
Image(gf_monthly$fill)
```

![](main_files/figure-markdown_github/do-gapfill-1.png)

``` r
gf_quarterly <- readRDS("quarterly_iMax5_140_gapfilled.rds")
Image(gf_quarterly$fill)
```

![](main_files/figure-markdown_github/do-gapfill-2.png)

BFAST Test
----------

``` r
# ts of a pixel
x <- as.vector(gf_monthly$fill[15,20,,])
# time must be given scaled to 1. -> monthly -> 1/12 = .08333333
# zoo must also be used in between, exactly as in ?bfastmonitor example
y <- as.ts(zoo(x, seq(2013, by = .08333333, length.out = 84))) 
bf <- bfastmonitor(y, start = 2019, order = 3, verbose = TRUE) # 2019,6 works but monitoring period is then off chart
plot(bf)

x1 <- as.vector(gf_quarterly$fill[15,20,,])
y1 <- as.ts(zoo(x1, seq(2013, by = .24999999, length.out = 28))) 
bf1 <- bfastmonitor(y1, start = 2019, order = 2, verbose = TRUE) # 2019,6 works but monitoring period is then off chart
plot(bf1)
```

Calculate BFAST Over AOI
------------------------

``` r
bfast_on_tile <- function(gapfill_matrix, by, ts, order) {
  dims <- dim(gapfill_matrix)
  result <- matrix(rep(FALSE, dims[1]*dims[2]), ncol = dims[1])
  for (i in 1:dims[1]) { # looping through x
    for (j in 1:dims[2]) { # loops through y
      raw_px_ts <- as.vector(gapfill_matrix[i,j,,])
      px_ts_obj <- as.ts(zoo(raw_px_ts, seq(2013, by = by, length.out = ts)))
      bfm_obj <- bfastmonitor(px_ts_obj, start = 2019, order = order)
      brkpoint <- bfm_obj$breakpoint
      if(!is.na(brkpoint)) {
        result[i,j] <- TRUE
      } else {
        # FALSE
      }
    }
  }
  return(result)
}

bfast_monthly <- bfast_on_tile(gf_monthly$fill, by = .08333333, ts = 84, order = 3)
bfast_quarter <- bfast_on_tile(gf_quarterly$fill, by = 0.25, ts = 28, order = 2)
# warning: too few observations in history period
```

``` r
# plot(prod["YEAR"])
ras <- rasterize(prod, as(st[,,,5], "Raster"), "YEAR")
ras.m <- aperm(matrix(ras[], ncol = 140), c(2,1))
ras.m[ras.m < 2019] <- FALSE
ras.m[ras.m == 2019] <- TRUE
ras.m[is.na(ras.m)] <- FALSE

Image(bfast_monthly)
```

![](main_files/figure-markdown_github/plot-bfast-aoi-1.png)

``` r
Image(bfast_quarter)
```

![](main_files/figure-markdown_github/plot-bfast-aoi-2.png)

``` r
Image(ras.m)
```

![](main_files/figure-markdown_github/plot-bfast-aoi-3.png)

``` r
table(bfast_monthly, ras.m)
```

    ##              ras.m
    ## bfast_monthly     0     1
    ##         FALSE 15646   162
    ##         TRUE   1687  2105

``` r
table(bfast_quarter, ras.m)
```

    ##              ras.m
    ## bfast_quarter     0     1
    ##         FALSE 15662   374
    ##         TRUE   1671  1893

old code
--------

``` r
# single pixel
defo_test <- st[st_geometry(prod[4,])]
defo_test_ts <- as.vector(defo_test[,3,3,][[1]])
test_ts <- as.ts(zoo(defo_test_ts, seq(2013, by = .08333333, length.out = 84))) 
test <- bfastmonitor(test_ts, start = c(2018, 1))
plot(test)

# mean of this polygon
mean_test <- st_apply(st[st_geometry(prod[2,])], "attributes", mean,  na.rm = TRUE) # aggregate
defo_test_ts <- mean_test[][[1]]
test_ts <- as.ts(zoo(defo_test_ts, seq(2013, by = .08333333, length.out = 84))) 
test <- bfastmonitor(test_ts, start = c(2018, 1))
plot(test)

# calculate breaks for every single pixel
vec <- c(0)
count <- 1
for (i in 1:5) {
  for (j in 1:12) {
    ts <- defo_test[,i,j,][[1]]
    test_ts <- as.ts(zoo(ts, seq(2013, by = .08333333, length.out = 84)))
    
    if (length(unique(test_ts)) > 1) {
      test_bf <- bfastmonitor(test_ts, start = c(2018, 1))
      vec[count] <- test_bf$breakpoint
    } else {
        # vec[count] <- 99
    }
    count <- count + 1
  }
}
vec
```

Results
=======

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
