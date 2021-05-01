
#----Load data----
data(snp_pos)
data(gene_ranges)
data(chr_length)

#----Start tests----
test_that("plot_snp_circos() plots SNP positions in the genome", {
    p <- plot_snp_circos(chr_length, gene_ranges, snp_pos)
    expect_equal(class(p), c("gg", "ggplot"))
})

test_that("plot_snp_distribution() plots SNP distribution", {
    p <- plot_snp_distribution(snp_pos)
    expect_equal(class(p), c("gg", "ggplot"))
})

