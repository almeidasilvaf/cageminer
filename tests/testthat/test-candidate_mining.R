
#----Load data------------------------------------------------------------------
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

#----Start tests----------------------------------------------------------------
test_that("handle_intervals() returns a GRanges object", {
    sim_markers <- snp_pos + 50
    h1 <- handle_intervals(snp_pos)
    h2 <- handle_intervals(sim_markers, expand_intervals = FALSE)

    expect_true("GRanges" %in% class(h1))
    expect_true("GRanges" %in% class(h2))
})

test_that("mine_step1() returns a character vector of gene IDs", {
    genes <- mine_step1(gene_ranges, snp_pos, window = 2)
    snp_list <- GenomicRanges::GRangesList(
        Trait1 = sample(snp_pos, 20),
        Trait2 = sample(snp_pos, 20)
    )
    genes2 <- mine_step1(gene_ranges, snp_list)

    expect_true(is(genes, "GRanges"))
    expect_true(
        "GRangesList" %in% class(genes2) |
            "CompressedGRangesList" %in% class(genes2)
    )
    expect_error(mine_step1(gene_ranges, as.data.frame(snp_pos)))

})

test_that("mine_step2() returns a list with 2 elements", {
    mine2 <- mine_step2(
        pepper_se, gcn = gcn, guides = guides$Gene,
        candidates = rownames(pepper_se)
    )

    expect_equal(class(mine2), "list")
    expect_equal(names(mine2), c("candidates", "enrichment"))

    expect_error(
        mine_step2(
            pepper_se, gcn = gcn, guides = guides$Gene[1:5],
            candidates = rownames(pepper_se)
        )
    )
})

test_that("mine_step3() returns a data frame of mined candidates", {
    mine3 <- mine_step3(
        pepper_se, candidates = mine2$candidates,
        sample_group = "PRR_stress"
    )

    expect_equal(class(mine3), "data.frame")
})

test_that("mine_candidates() returns a data frame", {
    mc1 <- mine_candidates(
        gene_ranges, snp_pos, exp = pepper_se, gcn = gcn,
        guides = guides$Gene, sample_group = "PRR_stress"
    )

    expect_equal(class(mc1), "data.frame")
    expect_equal(ncol(mc1), 5)
})

test_that("score_genes() returns a data frame", {
    hc_genes <- mined_candidates
    scored <- score_genes(hc_genes, hubs$Gene, tfs$Gene_ID)

    # Artitficially add scored genes as hubs and TFs to test weighing
    scored2 <- score_genes(
        hc_genes,
        c(hubs$Gene, "CA12g18400", "CA03g03310"),
        c(tfs$Gene_ID, "CA03g03310", "CA10g02780")
    )

    scored3 <- score_genes(
        hc_genes,
        c(hubs$Gene, "CA12g18400", "CA03g03310"),
        c(tfs$Gene_ID, "CA03g03310", "CA10g02780"),
        pick_top = 2
    )

    expect_equal(class(scored), "data.frame")
    expect_equal(class(scored2), "data.frame")
    expect_equal(nrow(scored3), 2)
    expect_error(
        score_genes(hc_genes)
    )
})
