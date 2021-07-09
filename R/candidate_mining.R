
#' Handle markers represented by intervals instead of SNPs
#'
#' @param marker_ranges Genomic positions of SNPs. For a single trait,
#' a GRanges object. For multiple traits, a GRangesList or CompressedGRangesList
#' object, with each element of the list representing SNP positions for a
#' particular trait.
#' @param window Sliding window (in Mb) upstream and downstream relative
#' to each SNP. Default: 2.
#' @param expand_intervals Logical indicating whether or not to expand markers
#' that are represented by intervals. This is particularly useful
#' if users want to use a custom interval defined by linkage disequilibrium,
#' for example. Default: TRUE.
#'
#' @return A GRanges object with genomic intervals relative to each SNP.
#' @importFrom GenomicRanges width
#' @noRd
handle_intervals <- function(marker_ranges, window = 2,
                             expand_intervals = TRUE) {
    window <- window * 10^6
    widths <- GenomicRanges::width(marker_ranges)
    if(any(widths > 1) & expand_intervals == FALSE) {
        indices <- which(widths == 1)
        marker_window <- marker_ranges
        marker_window[indices] <- marker_window[indices] + window
    } else {
        marker_window <- marker_ranges + window
    }
    return(marker_window)
}

#' Step 1: Get all putative candidate genes for a given sliding window
#'
#' For a user-defined sliding window relative to each SNP, this function will
#' subset all genes whose genomic positions overlap with the sliding window.
#'
#' @param gene_ranges A GRanges object with genomic coordinates
#' of all genes in the genome.
#' @param marker_ranges Genomic positions of SNPs. For a single trait,
#' a GRanges object. For multiple traits, a GRangesList or CompressedGRangesList
#' object, with each element of the list representing SNP positions for a
#' particular trait.
#' @param window Sliding window (in Mb) upstream and downstream relative
#' to each SNP. Default: 2.
#' @param expand_intervals Logical indicating whether or not to expand markers
#' that are represented by intervals. This is particularly useful
#' if users want to use a custom interval defined by linkage disequilibrium,
#' for example. Default: TRUE.
#'
#' @return A GRanges or GRangesList object with the genomic positions of
#' all putative candidate genes.
#' @seealso
#'  \code{\link[IRanges]{findOverlaps-methods}}
#' @rdname mine_step1
#' @export
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges GRangesList
#' @examples
#' data(snp_pos)
#' data(gene_ranges)
#' genes <- mine_step1(gene_ranges, snp_pos, window = 2)
mine_step1 <- function(gene_ranges, marker_ranges, window = 2,
                       expand_intervals = TRUE) {
    allowed_classes <- c("GRanges", "GRangesList", "CompressedGRangesList")
    if(!class(marker_ranges) %in% allowed_classes) {
        stop("Argument 'marker_ranges' must be a GRanges, GRangesList, or
             CompressedGRangesList object.")
    }
    if(is(marker_ranges, "GRanges")) {
        snp_granges <- handle_intervals(marker_ranges, window = window,
                                        expand_intervals = expand_intervals)
        genes <- unique(IRanges::subsetByOverlaps(gene_ranges, snp_granges))
    } else {
        snp_granges <- lapply(marker_ranges, function(x) {
          return(handle_intervals(x, window = window,
                                  expand_intervals = expand_intervals))
        })
        genes <- lapply(snp_granges, function(x) {
            return(unique(IRanges::subsetByOverlaps(gene_ranges, x)))
        })
        genes <- GenomicRanges::GRangesList(genes)
    }
    return(genes)
}

#' Step 2: Get candidates in modules enriched in guide genes
#'
#' @param exp Expression data frame with genes in row names and samples in
#' column names or a SummarizedExperiment object.
#' @param gcn Gene coexpression network returned by \code{BioNERO::exp2gcn()}.
#' @param guides Guide genes as a character vector or as a data frame with
#' genes in the first column and gene annotation class in the second column.
#' @param candidates Character vector of all candidates genes to be inspected.
#'
#' @return A list of 2 elements:
#' \describe{
#'   \item{candidates}{Character vector of candidates after step 2}
#'   \item{enrichment}{Data frame of results for enrichment analysis}
#' }
#' @rdname mine_step2
#' @export
#' @importFrom BioNERO module_enrichment
#' @examples
#' data(pepper_se)
#' data(snp_pos)
#' data(gene_ranges)
#' data(guides)
#' data(gcn)
#' set.seed(1)
#' mine2 <- mine_step2(pepper_se, gcn = gcn, guides = guides$Gene,
#'                     candidates = rownames(pepper_se))
mine_step2 <- function(exp, gcn, guides, candidates) {
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
    key_modules <- unique(enrichment$Module)
    if(is.null(key_modules)) { stop("No modules enriched in guide genes.") }
    genes_f1 <- gcn$genes_and_modules$Genes[gcn$genes_and_modules$Modules %in%
                                                key_modules]
    genes_f1 <- genes_f1[genes_f1 %in% candidates]
    result_list <- list(candidates = genes_f1,
                        enrichment = enrichment)
    return(result_list)
}

#' Step 3: Select candidates based on gene significance
#'
#' @param exp Expression data frame with genes in row names and samples in
#' column names or a SummarizedExperiment object.
#' @param metadata Sample metadata with samples in row names and sample
#' information in the first column. Ignored if `exp` is a SummarizedExperiment
#' object, as the colData will be extracted from the object.
#' @param candidates Character vector of candidate genes to be inspected.
#' @param sample_group Level of sample metadata to be used for filtering
#' in gene-trait correlation.
#' @param min_cor Minimum correlation value for
#' \code{BioNERO::gene_significance()}. Default: 0.2
#' @param alpha Numeric indicating significance level. Default: 0.05
#' @param continuous Logical indicating whether metadata is continuous or not.
#' Default: FALSE
#' @return A data frame with mined candidate genes and their correlation to
#' the condition of interest.
#'
#' @importFrom BioNERO gene_significance
#' @rdname mine_step3
#' @export
#' @examples
#' data(pepper_se)
#' data(snp_pos)
#' data(gene_ranges)
#' data(guides)
#' data(gcn)
#' set.seed(1)
#' mine2 <- mine_step2(pepper_se, gcn = gcn, guides = guides$Gene,
#'                     candidates = rownames(pepper_se))
#' mine3 <- mine_step3(pepper_se, candidates = mine2$candidates,
#'                     sample_group = "PRR_stress")
mine_step3 <- function(exp, metadata, candidates, sample_group,
                       min_cor = 0.2, alpha = 0.05,
                       continuous = FALSE) {
    exp <- exp[candidates, , drop = FALSE]
    genes_f2 <- BioNERO::gene_significance(
        exp, genes = candidates, alpha = alpha, min_cor = min_cor,
        continuous_trait = continuous, metadata = metadata
    )
    genes_f2 <- genes_f2$filtered_corandp[
        genes_f2$filtered_corandp$trait %in% sample_group,
    ]
    genes_final <- genes_f2[order(-genes_f2$cor), ]
    return(genes_final)
}

#' Mine high-confidence candidate genes in a single step
#'
#' @param gene_ranges A GRanges object with genomic coordinates
#' of all genes in the genome.
#' @param marker_ranges Genomic positions of SNPs. For a single trait,
#' a GRanges object. For multiple traits, a GRangesList or CompressedGRangesList
#' object, with each element of the list representing SNP positions for a
#' particular trait.
#' @param window Sliding window (in Mb) upstream and downstream relative
#' to each SNP. Default: 2.
#' @param expand_intervals Logical indicating whether or not to expand markers
#' that are represented by intervals. This is particularly useful
#' if users want to use a custom interval defined by linkage disequilibrium,
#' for example. Default: TRUE.
#' @param gene_col Column of the GRanges object containing gene ID.
#' Default: "ID", the default for gff/gff3 files imported with
#' rtracklayer::import.
#' @param exp Expression data frame with genes in row names and samples in
#' column names or a SummarizedExperiment object.
#' @param gcn Gene coexpression network returned by \code{BioNERO::exp2gcn()}.
#' @param guides Guide genes as a character vector or as a data frame with
#' genes in the first column and gene annotation class in the second column.
#' @param metadata Sample metadata with samples in row names and sample
#' information in the first column. Ignored if `exp` is a SummarizedExperiment
#' object, as the colData will be extracted from the object.
#' @param sample_group Level of sample metadata to be used for filtering
#' in gene-trait correlation.
#' @param min_cor Minimum correlation value for
#' \code{BioNERO::gene_significance()}. Default: 0.2
#' @param alpha Numeric indicating significance level. Default: 0.05
#' @param continuous Logical indicating whether metadata is continuous or not.
#' Default: FALSE
#' @return A data frame with mined candidate genes and their correlation to
#' the condition of interest.
#' @importFrom GenomicRanges mcols
#' @export
#' @rdname mine_candidates
#' @examples
#' data(pepper_se)
#' data(snp_pos)
#' data(gene_ranges)
#' data(guides)
#' data(gcn)
#' set.seed(1)
#' candidates <- mine_candidates(gene_ranges, snp_pos, exp = pepper_se,
#'                               gcn = gcn, guides = guides,
#'                               sample_group = "PRR_stress")
mine_candidates <- function(gene_ranges=NULL, marker_ranges=NULL, window = 2,
                            expand_intervals = TRUE,
                            gene_col = "ID",
                            exp=NULL, gcn=NULL, guides=NULL,
                            metadata, sample_group,
                            min_cor = 0.2, alpha = 0.05,
                            continuous = FALSE) {

    # Step 1
    candidates1 <- mine_step1(gene_ranges, marker_ranges, window = window,
                              expand_intervals = expand_intervals)
    candidates1 <- GenomicRanges::mcols(candidates1)[[gene_col]]

    # Step 2
    candidates2 <- mine_step2(exp, gcn, guides, candidates1)

    # Step 3
    candidates3 <- mine_step3(exp, metadata, candidates2$candidates,
                              sample_group, min_cor = min_cor, alpha = alpha,
                              continuous = continuous)
    return(candidates3)
}


#' Score candidate genes and select the top n genes
#'
#' @param mined_candidates Data frame resulting from \code{mine_candidates()}
#' or \code{mine_step()}.
#' @param hubs Character vector of hub genes.
#' @param tfs Character vector of transcription factors.
#' @param pick_top Number of top genes to select. Default: 10.
#' @return Data frame with top n candidates and their scores.
#' @export
#' @rdname score_genes
#' @examples
#' \donttest{
#' data(pepper_se)
#' data(snp_pos)
#' data(gene_ranges)
#' data(guides)
#' data(tfs)
#' set.seed(1)
#' # sft <- BioNERO::SFT_fit(pepper_se, net_type = "signed",
#' #                         cor_method = "pearson")
#' # Previously selected power = 12
#' gcn <- BioNERO::exp2gcn(pepper_se, net_type = "signed", cor_method = "pearson",
#'                         module_merging_threshold = 0.8, SFTpower = 12)
#' candidates <- mine_candidates(gene_ranges, snp_pos, exp = pepper_se,
#'                               gcn = gcn, guides = guides,
#'                               sample_group = "PRR_stress")
#' hubs <- BioNERO::get_hubs_gcn(pepper_se, gcn)
#' scored <- score_genes(candidates, hubs$Gene, tfs$Gene_ID)
#' }
score_genes <- function(mined_candidates, hubs=NULL, tfs=NULL,
                        pick_top=10) {
    if(is.null(hubs) & is.null(tfs)) {stop("Neither hubs nor TFs were provided.")}
    candidates <- mined_candidates
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









