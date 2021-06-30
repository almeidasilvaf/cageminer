
#----Load data----
data(snp_pos)
data(gene_ranges)
data(chr_length)

#----Start tests----
test_that("plot_snp_circos() plots SNP positions in the genome", {
    snp_pos_list <- GenomicRanges::GRangesList(
        snp_pos[1:20],
        snp_pos[21:40],
        snp_pos[41:60],
        snp_pos[61:80],
        snp_pos[81:100],
        snp_pos[101:116]
    )
    names(snp_pos_list) <- c("Trait1", "Trait2", "Trait3",
                              "Trait4", "Trait5", "Trait6")
    p1 <- plot_snp_circos(chr_length, gene_ranges, snp_pos)
    p2 <- plot_snp_circos(chr_length, gene_ranges, snp_pos_list)
    expect_true(methods::is(p1, "GGbio"))
    expect_true(methods::is(p2, "GGbio"))
})

test_that("plot_snp_distribution() plots SNP distribution", {
    p <- plot_snp_distribution(snp_pos)
    expect_equal(class(p), c("gg", "ggplot"))
})

