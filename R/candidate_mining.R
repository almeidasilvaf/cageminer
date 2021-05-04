


#' Mine high-confidence candidate genes
#'
#' @param exp Expression data frame with genes in row names and samples in
#' column names or a SummarizedExperiment object.
#' @param metadata Sample metadata with samples in row names and sample
#' information in the first column. Ignored if `exp` is a SummarizedExperiment
#' object, as the colData will be extracted from the object.
#' @param gcn Gene coexpression network returned by \code{BioNERO::exp2gcn()}.
#' @param guides Guide genes as a character vector or as a data frame with
#' genes in the first column and gene annotation class in the second column.
#' @param candidates Character vector of all candidates genes to be inspected.
#' @param min_cor Minimum correlation value for
#' \code{BioNERO::gene_significance()}. Default: 0.2
#' @param sample_group Level of sample metadata to be used for filtering
#' in gene-trait correlation.
#' @param alpha Numeric indicating significance level. Default: 0.05
#' @param continuous Logical indicating whether metadata is continuous or not.
#' Default: FALSE
#' @param verbose Logical indicating the verbosity of the function.
#' Default: TRUE
#' @return A data frame with high-confidence genes, their correlations to
#' `sample_group` and correlation p-values.
#' @details DETAILS
#' @examples
#' data(pepper_se)
#' data(snp_pos)
#' data(gene_ranges)
#' data(guides)
#' set.seed(1)
#' # Previously selected power = 12
#' sft <- BioNERO::SFT_fit(pepper_se, net_type = "signed",
#'                         cor_method = "pearson")
#' gcn <- BioNERO::exp2gcn(pepper_se, net_type = "signed", cor_method = "pearson",
#'                         module_merging_threshold = 0.8, SFTpower = 12)
#' hc_genes <- mine_candidates(pepper_se, gcn = gcn, guides = guides$Gene,
#'                             candidates = rownames(pepper_se),
#'                             sample_group = "PRR_stress")
#' @seealso
#'  \code{BioNERO::module_enrichment},
#'  \code{BioNERO::gene_significance}
#' @rdname mine_candidates
#' @export
#' @importFrom BioNERO module_enrichment gene_significance
#' @importFrom methods is
mine_candidates <- function(exp, metadata, gcn, guides, candidates,
                            sample_group,
                            min_cor = 0.2, alpha = 0.05,
                            continuous = FALSE,
                            verbose = TRUE) {
    if(is.character(guides)) {
        guides <- data.frame(Gene = guides, Class="guide")
    }
    bkgenes <- rownames(exp)
    annotation <- merge(guides, as.data.frame(bkgenes), by=1, all.y=TRUE)
    annotation[,2][is.na(annotation[,2])] <- "None"
    enrichment <- BioNERO::module_enrichment(
        gcn, background_genes = bkgenes, annotation = annotation
    )
    enrichment <- enrichment[enrichment$TermID != "None", ]
    if(verbose) {
        m <- lapply(seq_len(nrow(enrichment)), function(x) {
            message("Module ", enrichment[x, "Module"], ": ",
                    enrichment[x, "TermID"])
        })
    }
    key_modules <- unique(enrichment$Module)
    if(is.null(key_modules)) { stop("No enrichment. Exiting...")}
    genes_f1 <- gcn$genes_and_modules$Genes[gcn$genes_and_modules$Modules %in%
                                                key_modules]
    genes_f1 <- genes_f1[genes_f1 %in% candidates]
    exp_f1 <- exp[genes_f1, , drop = FALSE]
    genes_f2 <- BioNERO::gene_significance(
        exp_f1, genes = genes_f1, alpha = alpha, min_cor = min_cor,
        continuous_trait = continuous, metadata = metadata
    )
    genes_f2 <- genes_f2$filtered_corandp[
        genes_f2$filtered_corandp$trait %in% sample_group,
        ]
    genes_final <- genes_f2[order(-genes_f2$cor), ]
    return(genes_final)
}



#' Score candidate genes and select the top
#'
#' @param mined_candidates Data frame resulting from \code{mine_candidates}.
#' @param hubs Character vector of hub genes.
#' @param tfs Character vector of transcription factors.
#' @param pick_top Number of top genes to select. Default: 10.
#' @return Data frame with top n candidates.
#' @export
#' @rdname score_genes
#' @examples
#' data(pepper_se)
#' data(snp_pos)
#' data(gene_ranges)
#' data(guides)
#' data(tfs)
#' set.seed(1)
#' # Previously selected power = 12
#' sft <- BioNERO::SFT_fit(pepper_se, net_type = "signed",
#'                         cor_method = "pearson")
#' gcn <- BioNERO::exp2gcn(pepper_se, net_type = "signed", cor_method = "pearson",
#'                         module_merging_threshold = 0.8, SFTpower = 12)
#' hc_genes <- mine_candidates(pepper_se, gcn = gcn, guides = guides$Gene,
#'                             candidates = rownames(pepper_se),
#'                             sample_group = "PRR_stress")
#' hubs <- BioNERO::get_hubs_gcn(pepper_se, gcn)
#' scored <- score_genes(hc_genes, hubs$Gene, tfs$Gene_ID)
score_genes <- function(mined_candidates, hubs=NULL, tfs=NULL,
                        pick_top=10) {
    if(is.null(hubs) & is.null(tfs)) {stop("Neither hubs nor TFs were provided.")}
    candidates <- mined_candidates[!duplicated(mined_candidates$cor), ]
    cand_hubs <- candidates[candidates$gene %in% hubs, "gene"]
    cand_tfs <- candidates[candidates$gene %in% tfs, "gene"]
    cand_both <- intersect(cand_hubs, cand_tfs)
    scored <- candidates
    scored$score <- vapply(seq_len(nrow(scored)), function(x) {
        if(scored$gene[x] %in% cand_hubs) {
            y <- scored$cor[x] * 2
        } else if(scored$gene[x] %in% cand_tfs) {
            y <- scored$cor[x] * 2
        } else if(scored$gene[x] %in% cand_both) {
            y <- scored$cor[x] * 3
        } else {
            y <- scored$cor[x]
        }
        return(y)
    }, numeric(1))
    scored$abscore <- abs(scored$score)
    scored <- scored[order(-scored$abscore), ]
    scored <- scored[!duplicated(scored$gene), ]
    scored$abscore <- NULL
    if(nrow(scored) < pick_top) {
        message("Number of genes < 'pick_top'. Picking all genes.")
    } else {
        scored <- scored[1:pick_top, ]
    }
    return(scored)
}









