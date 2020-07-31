---
title: 'Lineages Through Space and Time in R: a worked example using the global distribution of dragon lizards (family Agamidae)'
author: "Alexander Skeels"

output: html_document
---






## Required Data

Generating LTST plots requires information on ancestral geographic states obtained from ancestral state estimation analysis using BioGeoBEARs. Once users have performed a BioGeoBEARS analysis, have selected a model of best-fit, and generated biogeographic stochastic maps from this model and its parameters, users should have the four objects necessary to continue.

First load the required packages



```r
library(BioGeoBEARS)
library(phytools)
library(reshape2)
library(ggplot2)
library(stringr)
library(devtools)
library(RColorBrewer)
```


Now install ltstR from the github repository




```r
devtools::install_github(repo = "alexskeels/ltstR")

library(ltstR)


# users can query the different functions to get a description of what they
# do. e.g.,

`?`(getEventTiming)

`?`(getRangeState)
```


```r
# users can alternatively load the R functions from the appendix

path <- "C:/path/to/appendix/directory"
setwd(path)

functions <- list("getLTSTDataTable.R","getLTSTFromStochasticMaps.R","getRangeState.R","LTST_vignette.html","slidingWindowForStochasticMaps.R","timeRangeAcrossStochasticMaps.R","getEventTiming.R")

lapply(functions, source)
```




For further information on each function or to get the code to modify these functions users can find the commented code in the R folder on the gitub page: 
www.github.com/alexskeels/ltstR

The ltstR package comes with the data necessary to create LTST plots for the agamid lizards.


```r
# first load the data
data(Agamidae)
```

```r
# or alternatively from the appenidx
setwd(path)
load("Agamidae.RData")
```

```r
# view the objects

ls()
```

```
## [1] "AET"       "CET"       "tipranges" "tr"
```

```r
# can have a look at the different objects

# tr is the phylogeny
# tipranges is the BioGeoBEARs geography object
# CET is the cladogenetic events table
# AET is the anagenetic events table

# the phylogeny
plot(tr, cex=0.1,  edge.width = 0.6)
```

<img src="ltstR_vignette_files/figure-html/unnamed-chunk-7-1.png" width="672" />

```r
# the geographic data
# W, I, P, and A are all single letter variable names associated with four biogeographic realms: 
# W = Wallacea/Australasia
# I = Indomalaya
# P = Palearctic
# A = Afrotropics
head(tipranges@df)
```

```
##                         W A P I
## Hydrosaurus_weberi      1 0 0 0
## Hydrosaurus_pustulatus  0 0 0 1
## Hydrosaurus_amboinensis 1 0 0 1
## Hypsilurus_papuensis    1 0 0 0
## Hypsilurus_bruijnii     1 0 0 0
## Hypsilurus_nigrigularis 1 0 0 0
```

```r
# the cladogenetic events table list(CET)
class(CET)
```

```
## [1] "list"
```

```r
class(CET[[1]])
```

```
## [1] "data.frame"
```

```r
head(CET[[1]][1:10])
```

```
##   node ord_ndname node_lvl node.type parent_br edge.length ancestor
## 1    1          1        6       tip         1    3.120884      426
## 2    2          2        6       tip         2    3.120884      426
## 3    3          3        5       tip       272    7.397286      425
## 4    4          4       11       tip       437   16.649613      433
## 5    5          5       12       tip       273   12.804597      434
## 6    6          6       13       tip         3   11.809053      435
##   daughter_nds  node_ht     time_bp
## 1              113.3192 3.90000e-05
## 2              113.3192 3.90000e-05
## 3              113.3192 3.90000e-05
## 4              113.3192 5.04270e-05
## 5              113.3192 5.04290e-05
## 6              113.3192 4.90214e-05
```

```r
# the anagenetic events table list(CET)
class(AET)
```

```
## [1] "list"
```

```r
class(AET[[1]])
```

```
## [1] "data.frame"
```

```r
head(AET[[1]][1:10])
```

```
##     node ord_ndname node_lvl node.type parent_br edge.length ancestor
## 3      3          3        5       tip       272    7.397286      425
## 93    93         93       12       tip       721   30.650273      522
## 96    96         96       16       tip        64    7.829888      526
## 99    99         99       17       tip       307   13.279947      529
## 107  107        107       17       tip       461   10.650706      538
## 124  124        124       18       tip       316    3.111069      551
##     daughter_nds  node_ht     time_bp
## 3                113.3192 3.90000e-05
## 93               113.3192 6.77341e-05
## 96               113.3192 5.75121e-05
## 99               113.3192 5.80816e-05
## 107              113.3192 5.75131e-05
## 124              113.3192 5.68981e-05
```

# LTST plot from a single biogeographic stochastic map (BSM)

Firstly we will load data from a single BSM and transform it into an LTST data table that has time steps as rows, and regions as columns, with each region/time-step being populated by a value of species diversity. In other words, species diversity through time for each region.







