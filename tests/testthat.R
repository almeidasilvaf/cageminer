Sys.setenv("R_TESTS" = "")
library(testthat)
library(cageminer)

test_check("cageminer")
