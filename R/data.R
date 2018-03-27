#' Array of all Metagenes
#'
#' list of metagenes
#'
#' @format array of 11 elements
#' @source \url{zezez}
"data.list"

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

#' Example of a cancer dataset
#'
#' A sample of 5000 genes and 59 samples from transcriptome  of lymph-node-negative
#' breast cancer published in Wang et al. (2005)
#'
#' Wang, Y., Klijn, J.G., Zhang, Y., Sieuwerts, A.M., Look, M.P., Yang, F., Talantov, D.,
#' Timmermans, M., Meijer-Van Gelder, M.E., Yu, J., Jatkoe, T., Berns, E.M., Atkins, D.,
#' Foekens, J.A.: Gene-expression profiles to predict distant metastasis of
#' lymph-node-negative primary breast cancer. Lancet 365(9460) (2005) 671–679
#'
#' @format a dataframe with the
#' \describe{
#'   \item{rows}{5000}
#'   \item{columns}{60}
#'   \item{first column}{is related to GENE names}
#'  }
#' @source \url{http://www.thelancet.com/journals/lancet/article/PIIS0140-6736(05)17947-1}
"Example_ds"
