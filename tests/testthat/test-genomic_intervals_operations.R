
#----Load data------------------------------------------------------------------
data(snp_pos)
data(gene_ranges)

#----Start tests----------------------------------------------------------------
test_that("simulate_windows() plots number of genes per sliding window", {
    snp_pos_list <- GenomicRanges::GRangesList(
        snp_pos[1:20],
        snp_pos[21:40],
        snp_pos[41:60]
    )
    p <- simulate_windows(gene_ranges, snp_pos)
    p2 <- simulate_windows(gene_ranges, snp_pos_list)

    expect_equal(class(p), c("gg", "ggplot"))
    expect_equal(class(p2), c("gg", "ggplot"))
})
