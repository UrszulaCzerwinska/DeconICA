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
#'A a sample 60 randomly selected samples from transcriptome of inflammatory
#'breast cance (IBC). Data were centred and in tranformed in log2 before sampling
#'
#'Bekhouche I, Finetti P, Adelaïde J, Ferrari A et al. High-resolution
#'comparative genomic hybridization of inflammatory breast cancer and
#'identification of candidate genes. PLoS One 2011 Feb 9;6(2):e16950. PMID:
#'21339811
#'
#'@format a dataframe with the \describe{ \item{rows}{21320} \item{columns}{61}
#'  \item{first column}{is related to GENE names} }
#'@source
#'  \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23720}
#'
"Example_ds"

#'TIMER signatures
#'
#'Signatures of 9 cell types published as part of TIMER tool
#'
#'Li, Bo, et al. "Comprehensive analyses of tumor immunity: implications for cancer immunotherapy." Genome biology 17.1 (2016): 174.
#'
#'@format a dataframe with the
#' \describe{
#'   \item{module}{first column}
#'   \item{module length}{second column}
#'   \item{gene names}{third column}
#'   }
#'@source
#'  \url{http://cistrome.org/TIMER/}
#'
"TIMER_cellTypes"

#'Decomposition of a transcriptome data
#'
#'A dataset overdecomposed (into 100 components). Data were downloaded from GEO,
#'then \code{\link{run_fastica}} using MATLAB algorithm with stabilisation as
#'applied.
#'
#'Source:
#'Bekhouche I, Finetti P, Adelaïde J, Ferrari A et al. High-resolution
#'comparative genomic hybridization of inflammatory breast cancer and
#'identification of candidate genes. PLoS One 2011 Feb 9;6(2):e16950. PMID:
#'21339811
#'
#'@format a list containg \describe{
#'  \item{A}{A ICA matrix (sample scores)}
#'  \item{S}{S ICA matrix (gene scores)}
#'  \item{names}{gene names}
#'  \item{samples}{sample names}
#'  \item{counts}{raw counts (non centered)}
#'  \item{log.counts}{log2 counts (non centered)}
#'  }
#'@source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23720}
#'
"BEK_ica_overdecompose"
