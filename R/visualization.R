
#' Set custom ggplot2 theme
#'
#' @return Theme for ggplot graphics.
#' @noRd
#' @importFrom ggplot2 theme_bw theme element_text element_blank
ggplot_theme <- function() {
    theme <- ggplot2::theme_bw() +
        ggplot2::theme(
        plot.title = ggplot2::element_text(face="bold", size=12),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
    )
    return(theme)
}


#' Custom color palette
#'
#' @param pal Color palette to use. One of 1 or 2. Default: 1.
#' @return A color palette as a character vector.
#' @noRd
custom_pal <- function(pal=1) {
    if(pal == 1) {
        col <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
                 "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF",
                 "#7E6148FF", "#B09C85FF")
    } else if(pal == 2) {
        col <- c("#393B79FF", "#637939FF", "#8C6D31FF", "#843C39FF",
                 "#7B4173FF", "#5254A3FF", "#8CA252FF", "#BD9E39FF",
                 "#AD494AFF", "#A55194FF", "#6B6ECFFF", "#B5CF6BFF",
                 "#E7BA52FF", "#D6616BFF", "#CE6DBDFF", "#9C9EDEFF",
                 "#CEDB9CFF", "#E7CB94FF", "#E7969CFF", "#DE9ED6FF")
    } else {
        stop("pal must be 1 or 2.")
    }
    return(col)
}



#' Wrapper function to check input for plot_snp_circos()
#'
#' @param genome_ranges A GRanges object with chromosome lengths.
#' @param gene_ranges A GRanges object with genomic coordinates
#' of all genes in the genome.
#' @param marker_ranges Genomic positions of SNPs. For a single trait,
#' a GRanges object. For multiple traits, a GRangesList or CompressedGRangesList
#' object, with each element of the list representing SNP positions for a
#' particular trait.
#'
#' @return If input objects are not as expected, it will throw an error.
#' Otherwise, nothing happens.
#' @noRd
check_input_circos <- function(genome_ranges, gene_ranges, marker_ranges) {
    if(!methods::is(genome_ranges, "GRanges")) {
        stop("Argument 'genome_ranges' must be a GRanges object.")
    }
    if(!methods::is(gene_ranges, "GRanges")) {
        stop("Argument 'gene_ranges' must be a GRanges object.")
    }
    if(!class(marker_ranges) %in%
       c("GRanges", "GRangesList", "CompressedGRangesList")) {
        stop(
            "Argument 'marker_ranges' must be a GRanges, GRangesList or
            CompressedGRangesList object."
            )
    }
}


#' Wrapper function to add seqlengths to ranges for plot_snp_circos()
#'
#' @param genome_ranges A GRanges object with chromosome lengths.
#' @param gene_ranges A GRanges object with genomic coordinates
#' of all genes in the genome.
#' @param marker_ranges Genomic positions of SNPs. For a single trait,
#' a GRanges object. For multiple traits, a GRangesList or CompressedGRangesList
#' object, with each element of the list representing SNP positions for a
#' particular trait.
#'
#' @return The genome_ranges, gene_ranges, or marker_ranges object with
#' seqlengths included.
#' @importFrom GenomeInfoDb seqlevels seqlengths keepSeqlevels
#' @importFrom GenomicRanges end
#' @noRd
add_seqlen <- function(genome_ranges, gene_ranges, marker_ranges) {
    GenomeInfoDb::seqlengths(gene_ranges) <- GenomicRanges::end(genome_ranges)
    GenomeInfoDb::seqlengths(genome_ranges) <- GenomicRanges::end(genome_ranges)

    levels <- GenomeInfoDb::seqlevels(marker_ranges)
    filt_genome <- GenomeInfoDb::keepSeqlevels(genome_ranges, levels,
                                               pruning.mode = "tidy")
    GenomeInfoDb::seqlengths(marker_ranges) <- GenomicRanges::end(filt_genome)
    result <- list(genome = genome_ranges,
                   genes = gene_ranges,
                   markers = marker_ranges)
    return(result)
}

#' Circos plot of SNP distribution across chromosomes
#'
#' @param genome_ranges A GRanges object with chromosome lengths.
#' @param gene_ranges A GRanges object with genomic coordinates
#' of all genes in the genome.
#' @param marker_ranges Genomic positions of SNPs. For a single trait,
#' a GRanges object. For multiple traits, a GRangesList or CompressedGRangesList
#' object, with each element of the list representing SNP positions for a
#' particular trait.
#'
#' @return A ggplot object with a circos plot of molecular marker distribution
#' across chromosomes.
#' @rdname plot_snp_circos
#' @export
#' @importFrom rlang .data
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom ggplot2 labs theme margin aes scale_fill_gradient
#' @importFrom ggtext element_textbox_simple
#' @importFrom GenomicRanges end tileGenome countOverlaps
#' @importFrom ggbio circle ggbio
#' @examples
#' data(snp_pos)
#' data(gene_ranges)
#' data(chr_length)
#' p <- plot_snp_circos(chr_length, gene_ranges, snp_pos)
plot_snp_circos <- function(genome_ranges, gene_ranges, marker_ranges) {
    check_input_circos(genome_ranges, gene_ranges, marker_ranges)
    pal <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF",
             "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF")
    gene_ranges <- add_seqlen(genome_ranges, gene_ranges, marker_ranges)$genes
    genome_ranges <- add_seqlen(genome_ranges,gene_ranges,marker_ranges)$genome
    marker_ranges <- add_seqlen(genome_ranges,gene_ranges,marker_ranges)$markers
    windows <- GenomicRanges::tileGenome(GenomeInfoDb::seqinfo(gene_ranges),
                                         tilewidth = 10^6,
                                         cut.last.tile.in.chrom = TRUE)
    windows$ngenes <- GenomicRanges::countOverlaps(windows, gene_ranges)
    p <- ggbio::ggbio()
    if(!is(marker_ranges, "GRanges")) {
        traits <- paste0(unlist(lapply(seq_along(marker_ranges), function(i) {
            return(paste0(
                "<span style='color:", pal[i], "'>",
                names(marker_ranges)[i],"</span>"))
        })), collapse=", ")
        for(n in seq_along(marker_ranges)) {
            p <- p + ggbio::circle(marker_ranges[[n]], geom = "rect", color = pal[n])
        }
    } else {

        traits <- "trait"
        p <- p + ggbio::circle(marker_ranges, geom = "rect", color = pal[1])
    }
    title <- paste0(
        "<b>SNP distribution across chromosomes</b><br>",
        "<span style = 'font-size:10pt'>Gene density and SNPs associated with ",
        traits, ".</span>")
    p <- p +
        ggbio::circle(
            windows, geom = "bar", stat="identity",
            aes(fill = .data$ngenes, y = .data$ngenes),
            color = NA
        ) +
        ggplot2::scale_fill_gradient(low="#74c476", high="darkgreen") +
        ggbio::circle(genome_ranges, geom = "scale", size = 2) +
        ggbio::circle(
            genome_ranges, geom = "text",
            aes(label = .data$seqnames), vjust = -1, size = 3
        ) +
        labs(title = title, fill = "Gene frequency") +
        theme(
            plot.title.position = "plot",
            plot.title = ggtext::element_textbox_simple(
                size = 13,
                lineheight = 1,
                padding = ggplot2::margin(5.5, 5.5, 5.5, 5.5),
                margin = ggplot2::margin(0, 0, 5.5, 0),
                )
        )

    return(p)
}

#' Plot a barplot of SNP distribution across chromosomes
#'
#' @param marker_ranges Genomic positions of SNPs. For a single trait,
#' a GRanges object. For multiple traits, a GRangesList or CompressedGRangesList
#' object, with each element of the list representing SNP positions for a
#' particular trait. List elements must have names for proper labelling.
#' @return A ggplot object.
#' @rdname plot_snp_distribution
#' @export
#' @importFrom ggplot2 ggplot aes geom_bar facet_wrap coord_flip theme_bw
#' ggtitle theme element_text element_blank scale_fill_manual
#' @importFrom methods is
#' @examples
#' data(snp_pos)
#' p <- plot_snp_distribution(snp_pos)
plot_snp_distribution <- function(marker_ranges) {
    if(!class(marker_ranges) %in% c("GRanges","GRangesList",
                                    "CompressedGRangesList")) {
        stop("Argument 'marker_ranges' must be a GRanges,
             GRangesList, or CompressedGRangesList object.")
    }
    if(methods::is(marker_ranges, "GRanges")) {
        marker_df <- as.data.frame(marker_ranges)
        marker_df <- as.data.frame(table(marker_df$seqnames))
        colnames(marker_df) <- c("Chromosome", "Frequency")
        wrap <- NULL
        bar <- ggplot2::geom_bar(stat = "identity")
        cols <- NULL
    } else {
        marker_df <- lapply(marker_ranges, as.data.frame)
        marker_df <- Reduce(rbind, lapply(seq_along(marker_df), function(x) {
            return(cbind(marker_df[[x]], names(marker_df)[x]))
        }))
        marker_df <- as.data.frame(table(marker_df[[ncol(marker_df)]],
                                         marker_df[["seqnames"]]))
        colnames(marker_df) <- c("Trait", "Chromosome", "Frequency")
        ntr <- nlevels(marker_df$Trait)
        wrap <- facet_wrap(~Trait, ncol = ntr)
        bar <- geom_bar(aes(fill = .data$Trait), stat = "identity",
                                 show.legend = FALSE)
        cols <- ggplot2::scale_fill_manual(values = custom_pal(2)[seq_len(ntr)])
    }
    p <- ggplot(marker_df, aes(x = .data$Chromosome, y = .data$Frequency)) +
        bar + cols + wrap +
        coord_flip() +
        ggplot2::ggtitle("SNP distribution across chromosomes") +
        ggplot_theme()
    return(p)
}











