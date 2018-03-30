#' Array of all Metagenes
#'
#' list of metagenes
#'
#' @format list of 11 elements
#' @source \url{http://www.cell.com/cell-reports/abstract/S2211-1247(14)00904-8}
"Biton.list"

#' Array of all Metagenes
#'
#' list of 22 immune cell type profiles
#'
#' @format list of 22 \code{data.frames}
#' @source \url{http://www.cell.com/cell-reports/abstract/S2211-1247(14)00904-8}
"LM22.list"

#' Cell type signatures
#'
#' Imported in correct format with \code{\link[ACSNMineR]{format_from_gmt}} and parsed
#'
#' @format a dataframe with the
#' \describe{
#'   \item{module}{first column}
#'   \item{module length}{second column}
#'   \item{gene names}{third column}
#'   }
#' @source \url{http://www.immgen.org}
"ImmgenHUGO"

#'Example of a cancer dataset
#'
#'A sample of 10000 genes and 59 randomly selected samples from transcriptome  of inflammatory
#'breast cance (IBC)
#'
#'Bekhouche I, Finetti P, Adelaïde J, Ferrari A et al. High-resolution
#'comparative genomic hybridization of inflammatory breast cancer and
#'identification of candidate genes. PLoS One 2011 Feb 9;6(2):e16950. PMID:
#'21339811
#'
#'@format a dataframe with the \describe{ \item{rows}{10000} \item{columns}{60}
#'  \item{first column}{is related to GENE names} }
#'@source
#'  \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23720}
#'
"Example_ds"
