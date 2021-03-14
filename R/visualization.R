
#' Wrapper for the circos plot
#'
#' @param genome_ranges A GRanges object with chromosome lengths.
#' @param genes_ranges A GRanges object with genomic coordinates
#' of all genes in the genome.
#' @param marker_ranges A GRanges object with genomic positions of
#' molecular markers.
#' @return A base plot to be converted to ggplot in \code{plot_snp_circos()}.
#'
#' @noRd
#' @seealso
#'  \code{\link[circlize]{circos.par}},
#'  \code{\link[circlize]{circos.genomicInitialize}},
#'  \code{\link[circlize]{circos.genomicDensity}},
#'  \code{\link[circlize]{circos.genomicTrack}},
#'  \code{\link[circlize]{circos.genomicRect}}
#' @importFrom circlize circos.par circos.genomicInitialize circos.genomicDensity circos.genomicTrack circos.genomicRect
circos_plot <- function(genome_ranges, genes_ranges, marker_ranges) {
    pal <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF",
             "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF")
    if(is.list(marker_ranges)) {
        if(length(marker_ranges) > 6) {
            stop("The input GRangesList contains more than 6 elements.
                 Reduce the number of elements in the list.")
        }
        marker_rangesdf <- lapply(marker_ranges, function(x) {
            return(as.data.frame(x[!duplicated(x)]))
        })
        height <- 0.2 - 0.035 * length(marker_ranges)
        circlize::circos.par("track.height" = height)
        circlize::circos.genomicInitialize(genome_ranges)
        circlize::circos.genomicDensity(as.data.frame(genes_ranges), col="grey50")
        x <- lapply(seq_along(marker_rangesdf), function(x) {
            circlize::circos.genomicTrack(marker_rangesdf[[x]], ylim=c(0,1),
                                          panel.fun=function(region, value, ...) {
                                              circlize::circos.genomicRect(
                                                  region, value,
                                                  ytop = 1, ybottom = 0,
                                                  border=pal[x], col=pal[x]
                                              )
                                          })
        })
    } else {
        marker_rangesdf <- as.data.frame(marker_ranges)
        circlize::circos.par("track.height" = 0.2)
        circlize::circos.genomicInitialize(genome_ranges)
        circlize::circos.genomicDensity(as.data.frame(genes_ranges), col="grey50")
        circlize::circos.genomicTrack(marker_rangesdf, ylim=c(0,1),
                                      panel.fun=function(region, value, ...) {
                                          circlize::circos.genomicRect(
                                              region, value,
                                              ytop = 1, ybottom = 0,
                                              border=pal[1], col=pal[1]
                                          )
                                      })
    }
    return(NULL)
}



#' Circos plot of molecular marker distribution across chromosomes
#'
#' @param genome_ranges A GRanges object with chromosome lengths.
#' @param genes_ranges A GRanges object with genomic coordinates
#' of all genes in the genome.
#' @param marker_ranges A GRanges object with genomic positions of
#' molecular markers.
#'
#' @return A ggplot object with a circos plot of molecular marker distribution
#' across chromosomes.
#' @rdname plot_snp_circos
#' @seealso
#'  \code{\link[ggplotify]{as.ggplot}}
#'  \code{\link[ggtext]{element_textbox}}
#' @export
#' @importFrom ggplotify as.ggplot
#' @importFrom ggplot2 labs theme margin
#' @importFrom ggtext element_textbox_simple
#' @examples
#' data(gwas)
#' data(maize_gr)
#' data(chr_length)
#' gwas_list <- split(gwas, gwas$trait)
#' genes_ranges <- maize_gr[maize_gr$type == "gene", ]
#' p <- plot_snp_circos(chr_length, genes_ranges, gwas_list)
plot_snp_circos <- function(genome_ranges, genes_ranges, marker_ranges) {

    pal <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF",
             "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF")

    p <- ggplotify::as.ggplot(function() circos_plot(
        genome_ranges, genes_ranges, marker_ranges)
        )

    traits <- unlist(lapply(seq_along(marker_ranges), function(i) {
        return(paste0(
            "<span style='color:", pal[i], "'>", names(marker_ranges)[i],"</span>"
        ))
    }))
    title <- paste0(
        "<b>SNP distribution across chromosomes</b><br>",
        "<span style = 'font-size:10pt'>",
        "<span style='color:#404040'>Gene density </span>and SNPs associated with ",
        paste(traits, collapse=", "),
        ".</span>")

    p2 <- p +
        ggplot2::labs(title = title) +
        ggplot2::theme(plot.title.position = "plot",
                       plot.title = ggtext::element_textbox_simple(
                           size = 13,
                           lineheight = 1,
                           padding = ggplot2::margin(5.5, 5.5, 5.5, 5.5),
                           margin = ggplot2::margin(0, 0, 5.5, 0),
                       ))
    return(p2)
}

