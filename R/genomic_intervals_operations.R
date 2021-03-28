

#' Simulate number of genes for each sliding window
#'
#' This function counts genes that are contained in sliding windows related to
#' each SNP.
#'
#' @param genes_ranges A GRanges object with genomic coordinates
#' of all genes in the genome.
#' @param marker_ranges A GRanges or GRangesList object with
#' positions of molecular markers.
#' @param windows Sliding windows (in Mb) upstream and downstream relative
#' to each SNP. Default: seq(0.1, 2, by = 0.1).
#' @return A ggplot object summarizing the results of the simulations.
#' @details
#' By default, the function creates 20 sliding windows by expanding upstream
#' and downstream boundaries for each SNP from 0.1 Mb (100 kb) to 2 Mb.
#' @examples
#' data(gwas)
#' data(maize_gr)
#' genes_ranges <- maize_gr[maize_gr$type == "gene", ]
#' marker_ranges <- split(gwas, gwas$trait)
#' simulate_windows(genes_ranges, marker_ranges)
#' @seealso
#'  \code{\link[IRanges]{findOverlaps-methods}}
#' @rdname simulate_windows
#' @export
#' @importFrom IRanges subsetByOverlaps
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_manual labs theme element_text
#' @importFrom methods is
simulate_windows <- function(genes_ranges, marker_ranges,
                             windows = seq(0.1, 2, by=0.1)) {

    windows <- windows * 10^6
    if(is(marker_ranges, "GRanges")) {
        sim <- lapply(windows, function(y) {
            return(marker_ranges + y)
        })
        gene_count <- data.frame(count=unlist(lapply(sim, function(x) {
            return(length(unique(IRanges::subsetByOverlaps(genes_ranges, x))))
        })))
        legend_pos <- "none"
    } else {
        sim <- lapply(marker_ranges, function(x) {
            runs <- lapply(windows, function(y) {
                return(x + y)
            })
            return(runs)
        })
        gene_count <- as.data.frame(lapply(sim, function(x) {
            count <- unlist(lapply(x, function(y) {
                return(length(unique(IRanges::subsetByOverlaps(genes_ranges, y))))
            }))
            return(count)
        }))
        legend_pos <- "bottom"
    }
    gene_count$windows <- as.factor(seq(0.1, 2, by = 0.1))
    gene_count_melt <- reshape2::melt(gene_count)

    p <- ggplot2::ggplot(gene_count_melt,
                         ggplot2::aes(x=windows, y=value,
                                      group=variable, color=variable)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::scale_color_manual(values = custom_pal(1)) +
        ggplot2::labs(x="Sliding windows (Mb)", y="Gene count",
                      title="Genes per sliding window",
                      color="Trait") +
        ggplot_theme() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1),
                       legend.position = legend_pos)
    return(p)
}


#' Get candidate genes for a given sliding window
#'
#' For an user-defined sliding window relative to each SNP, this function will
#' subset all genes whose genomic positions overlap with the sliding window.
#'
#' @param genes_ranges A GRanges object with genomic coordinates
#' of all genes in the genome.
#' @param marker_ranges A GRanges or GRangesList object with
#' positions of molecular markers.
#' @param window Sliding window (in Mb) upstream and downstream relative
#' to each SNP. Default: 2.
#' @return A GRanges or GRangesList object with the genomic positions of
#' candidate genes.
#' @examples
#' data(gwas)
#' data(maize_gr)
#' genes_ranges <- maize_gr[maize_gr$type == "gene", ]
#' marker_ranges <- split(gwas, gwas$trait)
#' genes <- get_all_candidates(genes_ranges, marker_ranges, window = 2)
#' @seealso
#'  \code{\link[IRanges]{findOverlaps-methods}}
#' @rdname get_all_candidates
#' @export
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges GRangesList
get_all_candidates <- function(genes_ranges, marker_ranges, window = 2) {
    window <- window * 10^6
    if(is(marker_ranges, "GRanges")) {
        snp_granges <- marker_ranges + window
        genes <- unique(IRanges::subsetByOverlaps(genes_ranges, snp_granges))
    } else {
        snp_granges <- lapply(marker_ranges, function(x) return(x + window))
        genes <- lapply(snp_granges, function(x) {
            return(unique(IRanges::subsetByOverlaps(genes_ranges, x)))
        })
        genes <- GenomicRanges::GRangesList(genes)
    }
    return(genes)
}



