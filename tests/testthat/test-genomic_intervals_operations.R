
#----Load data----
data(gwas)
data(maize_gr)

genes_ranges <- maize_gr[maize_gr$type == "gene", ]
marker_ranges <- split(gwas, gwas$trait)

#----Start tests----
test_that("simulate_windows() plots number of genes per sliding window", {
    p <- simulate_windows(genes_ranges, marker_ranges)
    expect_equal(class(p), c("gg", "ggplot"))
})
