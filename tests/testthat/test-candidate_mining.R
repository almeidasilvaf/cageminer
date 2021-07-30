
#----Load data----
data(pepper_se)
data(snp_pos)
data(gene_ranges)
data(mined_candidates)
data(guides)
data(mine2)
data(tfs)
data(gcn)
data(hubs)

# Get candidates
genes <- mine_step1(gene_ranges, snp_pos, window = 2)$ID

# Create expression with only candidates
set.seed(1)

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
#     mine3 <- mine_step3(pepper_se, candidates = mine2$candidates,
#                         sample_group = "PRR_stress")
#     expect_equal(class(mine3), "data.frame")
# })

test_that("score_genes() returns a data frame", {
    hc_genes <- mined_candidates
    scored <- score_genes(hc_genes, hubs$Gene, tfs$Gene_ID)
    expect_equal(class(scored), "data.frame")
})
