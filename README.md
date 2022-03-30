
# RNAseqStat2

<!-- badges: start -->
<!-- badges: end -->

The goal of RNAseqStat2 is to do DEG analysis

## Installation

You can install the development version of RNAseqStat2 from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xiayh17/RNAseqStat2")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(RNAseqStat2)
library(airway)

data("airway")

row_counts <- as.data.frame(assay(airway))

group_list <- as.character(colData(airway)$dex)

data_i <- Create_DEGContainer(species = "Human",
                              dataType = "Counts",
                              idType = "ENSEMBL",
                              expMatrix = row_counts,
                              groupInfo = group_list,
                              caseGroup = "trt")

# run all in one line
data_o <- runALL(object = data_i,dir = "output_test")
```

更多参考文档
https://www.yuque.com/xiayonghe/rnaseqstat2
