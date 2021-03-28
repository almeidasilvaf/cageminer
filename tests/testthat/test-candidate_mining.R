
#----Load data----
data(maize_exp)
data(gwas)
data(maize_gr)
data(guides)

genes_ranges <- maize_gr[maize_gr$type == "gene", ]
marker_ranges <- gwas[gwas$trait == "Sucrose"]

# Get candidates
genes <- get_all_candidates(genes_ranges, marker_ranges, window = 2)$ID

# Create expression with only candidates
exp <- maize_exp
exp <- BioNERO::exp_preprocess(exp, cor_method = "pearson", min_exp=5)
gcn <- BioNERO::exp2gcn(exp, net_type = "signed", cor_method = "pearson",
                        module_merging_threshold = 0.8, SFTpower = 12)

#----Start tests----
test_that("mine_candidates() returns high-confidence genes", {
    hc_genes <- mine_candidates(exp, gcn = gcn, guides = guides,
                                sample_group = "leaf")

})

