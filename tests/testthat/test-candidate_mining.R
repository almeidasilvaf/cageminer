
#----Load data----
data(pepper_se)
data(snp_pos)
data(gene_ranges)
data(guides)
data(tfs)

# Get candidates
genes <- mine_step1(gene_ranges, snp_pos, window = 2)$ID

# Create expression with only candidates
set.seed(1)
# power = 12
gcn <- BioNERO::exp2gcn(pepper_se, net_type = "signed", cor_method = "pearson",
                        module_merging_threshold = 0.8, SFTpower = 12)

hubs <- BioNERO::get_hubs_gcn(pepper_se, gcn)

#----Start tests----
test_that("mine_step1() returns a character vector of gene IDs", {
    genes <- mine_step1(gene_ranges, snp_pos, window = 2)
    expect_true(is(genes, "GRanges"))
})

test_that("mine_step2() returns a list with 2 elements", {
    mine2 <- mine_step2(pepper_se, gcn = gcn, guides = guides$Gene,
                        candidates = rownames(pepper_se))
    expect_equal(class(mine2), "list")
    expect_equal(names(mine2), c("candidates", "enrichment"))
})

# test_that("mine_step3() returns a data frame of mined candidates", {
#     mine2 <- mine_step2(pepper_se, gcn = gcn, guides = guides$Gene,
#                         candidates = rownames(pepper_se))
#     mine3 <- mine_step3(pepper_se, candidates = mine2$candidates,
#                         sample_group = "PRR_stress")
#     expect_equal(class(mine3), "data.frame")
# })


# test_that("mine_candidates() integrates mine_step* functions", {
#     hc_genes <- mine_candidates(gene_ranges, snp_pos, exp = pepper_se,
#                                 gcn = gcn, guides = guides,
#                                 sample_group = "PRR_stress")
#     expect_equal(class(hc_genes), "data.frame")
# })


test_that("score_genes() returns a data frame", {
    hc_genes <- mine_candidates(gene_ranges, snp_pos, exp = pepper_se,
                                gcn = gcn, guides = guides,
                                sample_group = "PRR_stress")
    scored <- score_genes(hc_genes, hubs$Gene, tfs$Gene_ID)
    expect_equal(class(scored), "data.frame")
})
