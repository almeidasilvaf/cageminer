
#----Load data------------------------------------------------------------------
data(snp_pos)
data(gene_ranges)
data(chr_length)

## Simulate GRangesList of multiple traits
snp_pos_list <- GenomicRanges::GRangesList(
    snp_pos[1:20],
    snp_pos[21:40],
    snp_pos[41:60],
    snp_pos[61:80],
    snp_pos[81:100],
    snp_pos[101:116]
)
names(snp_pos_list) <- c("Trait1", "Trait2", "Trait3", "Trait4",
                         "Trait5", "Trait6")


#----Start tests----------------------------------------------------------------
test_that("custom_pal() returns a color palette in a character vector", {
    pal1 <- custom_pal(1)
    pal2 <- custom_pal(2)

    expect_equal(class(pal1), "character")
    expect_equal(class(pal2), "character")
    expect_error(custom_pal(3))
})

test_that("check_input_circos() can properly check input data classes", {

    expect_error(
        check_input_circos(as.data.frame(chr_length))
    )

    expect_error(
        check_input_circos(chr_length, gene_ranges = as.data.frame(gene_ranges))
    )

    expect_error(
        check_input_circos(
            chr_length, gene_ranges, as.data.frame(snp_pos)
        )
    )
})

test_that("plot_snp_circos() plots SNP positions in the genome", {
    p1 <- plot_snp_circos(chr_length, gene_ranges, snp_pos)
    p2 <- plot_snp_circos(chr_length, gene_ranges, snp_pos_list)
    expect_true(methods::is(p1, "GGbio"))
    expect_true(methods::is(p2, "GGbio"))
})

test_that("plot_snp_distribution() plots SNP distribution", {
    p <- plot_snp_distribution(snp_pos)
    p2 <- plot_snp_distribution(snp_pos_list)

    expect_equal(class(p), c("gg", "ggplot"))
    expect_equal(class(p2), c("gg", "ggplot"))
    expect_error(plot_snp_distribution(as.data.frame(snp_pos)))
})

