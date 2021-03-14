
#----Load data----
data(gwas)
data(maize_gr)
data(chr_length)

#----Start tests----
test_that("plot_circos() plots SNP positions in the genome", {
    gwas_list <- split(gwas, gwas$trait)
    genes_ranges <- maize_gr[maize_gr$type == "gene", ]
    p <- plot_snp_circos(chr_length, genes_ranges, gwas_list)
    expect_equal(class(p), c("gg", "ggplot"))
})



