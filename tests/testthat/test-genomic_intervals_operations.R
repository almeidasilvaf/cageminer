
#----Load data----
data(snp_pos)
data(gene_ranges)

#----Start tests----
test_that("simulate_windows() plots number of genes per sliding window", {
    p <- simulate_windows(gene_ranges, snp_pos)
    expect_equal(class(p), c("gg", "ggplot"))
})


test_that("get_all_candidates() returns a character vector of gene IDs", {
    genes <- get_all_candidates(gene_ranges, snp_pos, window = 2)
    expect_true(is(genes, "GRanges"))
})
