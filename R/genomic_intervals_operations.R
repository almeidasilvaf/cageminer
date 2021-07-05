

#' Simulate number of genes for each sliding window
#'
#' This function counts genes that are contained in sliding windows related to
#' each SNP.
#'
#' @param gene_ranges A GRanges object with genomic coordinates
#' of all genes in the genome.
#' @param marker_ranges Genomic positions of SNPs. For a single trait,
#' a GRanges object. For multiple traits, a GRangesList or CompressedGRangesList
#' object, with each element of the list representing SNP positions for a
#' particular trait.
#' @param windows Sliding windows (in Mb) upstream and downstream relative
#' to each SNP. Default: seq(0.1, 2, by = 0.1).
#' @param expand_intervals Logical indicating whether or not to expand markers
#' that are represented by intervals. This is particularly useful
#' if users want to use a custom interval defined by linkage disequilibrium,
#' for example. Default: TRUE.
#'
#' @return A ggplot object summarizing the results of the simulations.
#' @details
#' By default, the function creates 20 sliding windows by expanding upstream
#' and downstream boundaries for each SNP from 0.1 Mb (100 kb) to 2 Mb.
#' @seealso
#'  \code{\link[IRanges]{findOverlaps-methods}}
#' @rdname simulate_windows
#' @export
#' @importFrom IRanges subsetByOverlaps
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes_ geom_line geom_point scale_color_manual labs theme element_text
#' @importFrom methods is
#' @examples
#' data(snp_pos)
#' data(gene_ranges)
#' simulate_windows(gene_ranges, snp_pos)
simulate_windows <- function(gene_ranges, marker_ranges,
                             windows = seq(0.1, 2, by=0.1),
                             expand_intervals = TRUE) {
    if(is(marker_ranges, "GRanges")) {
        sim <- lapply(windows, function(y) {
            return(handle_intervals(marker_ranges, y,
                                    expand_intervals = TRUE))
        })
        gene_count <- data.frame(count=unlist(lapply(sim, function(x) {
            return(length(unique(IRanges::subsetByOverlaps(gene_ranges, x))))
        })))
        legend_pos <- "none"
    } else {
        sim <- lapply(marker_ranges, function(x) {
            runs <- lapply(windows, function(y) {
                return(handle_intervals(x, y,
                                        expand_intervals = TRUE))
            })
            return(runs)
        })
        gene_count <- as.data.frame(lapply(sim, function(x) {
            count <- unlist(lapply(x, function(y) {
                return(length(unique(IRanges::subsetByOverlaps(gene_ranges, y))))
            }))
            return(count)
        }))
        legend_pos <- "bottom"
    }
    gene_count$windows <- as.factor(seq(0.1, 2, by = 0.1))
    gene_count_melt <- reshape2::melt(gene_count, id.vars = "windows")

    p <- ggplot2::ggplot(gene_count_melt,
                         ggplot2::aes_(x=~windows, y=~value,
                                      group=~variable, color=~variable)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::scale_color_manual(values = custom_pal(1)) +
        ggplot2::labs(x="Sliding windows (Mb)", y="Gene count",
                      title="Genes per sliding window",
                      color="Trait") +
        ggplot_theme() +
        ggplot2::theme(legend.position = legend_pos)
    return(p)
}
