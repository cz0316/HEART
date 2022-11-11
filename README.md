---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# HEART

<!-- badges: start -->
<!-- badges: end -->

HEART, a statistical combination test, 
	can easily capture the various source of the difference beyond mean 
	expression changesand facilitate the accuracy, robustness, 
	and efficiency of single-cell DE analysis.

## Installation


First, install the packages: data.table, fBasics, matrixTests.
Then, install the package HEART.

``` r
install.packages("data.table")
install.packages("fBasics")
install.packages("matrixTests")
install.packages("~path/HEART_1.1.0.tar.gz")


```




## Example

This is a basic example which shows you how to use HEART.
```{r example}
library(data.table)
library(fBasics)
library(matrixTests)
library(HEART)
data("example")

expression[1:5,1:5]

head(cell.data)
```

"example" data has two datasets: 
  expression: a counts matrix.
  cell.data: a data frame containing information on individual cells.






```{r }

result=heart(counts = expression,cell.data = cell.data)

head(result)

```

