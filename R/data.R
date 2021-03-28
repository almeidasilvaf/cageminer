
#' GWAS-derived SNPs for maize traits
#'
#' The SNPs in this data set were retrieved from Wallace et al., 2014, and
#' they are associated to chlorophyll A, chlorophyll B, nitrate and sucrose.
#' The SNP positions and associated traits are stored in a GRanges object.
#'
#' @name gwas
#' @format A GRanges object.
#' @references
#' Wallace JG, Bradbury PJ, Zhang N, Gibon Y, Stitt M, Buckler ES (2014)
#' Association Mapping across Numerous Traits Reveals Patterns of Functional
#' Variation in Maize. PLoS Genet 10(12): e1004845.
#' doi:10.1371/journal.pgen.1004845
#' @examples
#' data(gwas)
#' @usage data(gwas)
"gwas"


#' Genomic coordinates for maize features
#'
#' The genomic coordinates for maize features (e.g., genes, transcripts)
#' are stored in a GRanges object. Some columns of the .gff file were considered
#' unnecessary and, hence, they were removed to reduce package size. The .gff
#' file was obtained from PLAZA 4.0 Monocots.
#'
#' @name maize_gr
#' @format A GRanges object
#' @references
#' Van Bel, M., Diels, T., Vancaester, E., Kreft, L., Botzki, A.,
#' Van de Peer, Y., ... & Vandepoele, K. (2018). PLAZA 4.0: an integrative
#' resource for functional, evolutionary and comparative plant
#' genomics. Nucleic acids research, 46(D1), D1190-D1196.
#' @examples
#' data(maize_gr)
#' @usage data(maize_gr)
"maize_gr"


#' Maize chromosome lengths
#'
#' Lengths of maize chromosomes 1-10 calculated as a GRanges object.
#' The genome for which lengths were calculated (B73.v4) was downloaded from
#' PLAZA 4.0 Monocots.
#'
#' @name chr_length
#' @format A GRanges object
#' @references
#' Van Bel, M., Diels, T., Vancaester, E., Kreft, L., Botzki, A.,
#' Van de Peer, Y., ... & Vandepoele, K. (2018). PLAZA 4.0: an integrative
#' resource for functional, evolutionary and comparative plant
#' genomics. Nucleic acids research, 46(D1), D1190-D1196.
#' @examples
#' data(chr_length)
#' @usage data(chr_length)
"chr_length"


#' Maize expression data from Stelpflug et al., 2016
#'
#' The data were filtered to remove genes with median expression levels <2
#' for package size issues.
#'
#' @name maize_exp
#' @format A SummarizedExperiment object.
#' @references
#' Stelpflug, S.C., Sekhon, R.S., Vaillancourt, B., Hirsch, C.N., Buell,
#' C.R., de Leon, N. and Kaeppler, S.M. (2016), An Expanded Maize Gene
#' Expression Atlas based on RNA Sequencing and its Use to Explore Root
#' Development. The Plant Genome, 9: plantgenome2015.04.0025.
#' https://doi.org/10.3835/plantgenome2015.04.0025
#' @examples
#' data(maize_exp)
"maize_exp"


#' Guide genes associated with sucrose metabolism
#'
#' The MapMan annotation was downloaded from PLAZA Monocots.
#'
#' @name guides
#' @format A data frame with genes in the first column and MapMan in the
#' second column.
#' @references
#' Van Bel, M., Diels, T., Vancaester, E., Kreft, L., Botzki, A.,
#' Van de Peer, Y., ... & Vandepoele, K. (2018). PLAZA 4.0: an integrative
#' resource for functional, evolutionary and comparative plant
#' genomics. Nucleic acids research, 46(D1), D1190-D1196.
#' @examples
#' data(guides)
"guides"

