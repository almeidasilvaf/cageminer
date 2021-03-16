

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
#' @return OUTPUT_DESCRIPTION
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
