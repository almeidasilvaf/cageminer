
#' Capsicum annuum SNPs associated with resistance to Phytophthora root rot.
#'
#' The SNPs in this data set were retrieved from Siddique et al., 2019, and
#' they are associated to resistance to Phytophthora root rot.
#'
#' @name snp_pos
#' @format A GRanges object.
#' @references
#' Siddique, M.I., Lee, HY., Ro, NY. et al. Identifying candidate genes for
#' Phytophthora capsici resistance in pepper (Capsicum annuum) via
#' genotyping-by-sequencing-based QTL mapping and genome-wide association
#' study. Sci Rep 9, 9962 (2019). https://doi.org/10.1038/s41598-019-46342-1
#' @examples
#' data(snp_pos)
#' @usage data(snp_pos)
"snp_pos"


#' Genomic coordinates of pepper genes
#'
#' GRanges object with genomic coordinates of pepper genes downloaded from
#' http://peppergenome.snu.ac.kr/download.php.
#'
#' @name gene_ranges
#' @format A GRanges object
#' @examples
#' data(gene_ranges)
#' @usage data(gene_ranges)
"gene_ranges"


#' Pepper chromosome lengths
#'
#' Lengths of pepper chromosomes 1-12 in a GRanges object.
#' The genome for which lengths were calculated (v1.55) was downloaded from
#' http://peppergenome.snu.ac.kr/download.php
#'
#' @name chr_length
#' @format A GRanges object
#' @examples
#' data(chr_length)
#' @usage data(chr_length)
"chr_length"


#' Gene expression data from Kim et al., 2018.
#'
#' The data were filtered to keep only the top 4000 genes with highest RPKM
#' values in PRR stress-related samples.
#'
#' @name pepper_se
#' @format A SummarizedExperiment object.
#' @references
#' Kim, MS., Kim, S., Jeon, J. et al. Global gene expression profiling for
#' fruit organs and pathogen infections in the pepper, Capsicum annuum L..
#' Sci Data 5, 180103 (2018). https://doi.org/10.1038/sdata.2018.103
#' @examples
#' data(pepper_se)
"pepper_se"


#' Guide genes associated with defense and resistance to oomycetes
#'
#' The GO annotation was retrieved from PLAZA 4.0 Dicots.
#'
#' @name guides
#' @format A data frame with genes in the first column and GO description
#' in the second column.
#' @references
#' Van Bel, M., Diels, T., Vancaester, E., Kreft, L., Botzki, A.,
#' Van de Peer, Y., ... & Vandepoele, K. (2018). PLAZA 4.0: an integrative
#' resource for functional, evolutionary and comparative plant
#' genomics. Nucleic acids research, 46(D1), D1190-D1196.
#' @examples
#' data(guides)
"guides"

#' Pepper transcription factors
#'
#' Pepper transcription factors and their families retrieved from PlantTFDB 4.0.
#'
#' @name tfs
#' @format A data frame with gene IDs in the first column and TF families in
#' the second column.
#' @references
#' Jin, J., Tian, F., Yang, D. C., Meng, Y. Q., Kong, L., Luo, J., &
#' Gao, G. (2016). PlantTFDB 4.0: toward a central hub for transcription
#' factors and regulatory interactions in plants. Nucleic acids research, gkw982.
#' @examples
#' data(tfs)
"tfs"
