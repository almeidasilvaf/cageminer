


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
#' @return A data frame with high-confidence genes, their correlations to
#' `sample_group` and correlation p-values.
#' @details DETAILS
#' @examples
#' data(maize_exp)
#' data(gwas)
#' data(maize_gr)
#' data(guides)
#' genes_ranges <- maize_gr[maize_gr$type == "gene", ]
#' marker_ranges <- gwas[gwas$trait == "Sucrose"]
#' genes <- get_all_candidates(genes_ranges, marker_ranges, window = 2)$ID
#' exp <- maize_exp
#' exp <- BioNERO::exp_preprocess(exp, cor_method = "pearson", min_exp=5)
#' # SFT power previously calculated
#' gcn <- BioNERO::exp2gcn(exp, net_type = "signed", cor_method = "pearson",
#'                         module_merging_threshold = 0.8, SFTpower = 12)
#' @seealso
#'  \code{BioNERO::module_enrichment},
#'  \code{BioNERO::gene_significance}
#' @rdname mine_candidates
#' @export
#' @importFrom BioNERO module_enrichment gene_significance
mine_candidates <- function(exp, metadata, gcn, guides, candidates,
                            min_cor = 0.2, sample_group,
                            alpha = 0.05, continuous = FALSE) {
    if(is.character(guides)) {
        guides <- data.frame(Gene = guides, Class="guide")
    }
    bkgenes <- rownames(exp)
    annotation <- merge(guides, as.data.frame(bkgenes), by=1, all.y=TRUE)
    annotation[,2][is.na(annotation[,2])] <- "None"
    enrichment <- BioNERO::module_enrichment(
        gcn, background_genes = bkgenes, annotation = annotation
    )
    key_modules <- unique(enrichment$Module)
    genes_f1 <- gcn$genes_and_modules$Genes[gcn$genes_and_modules$Modules %in%
                                                key_modules]
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
