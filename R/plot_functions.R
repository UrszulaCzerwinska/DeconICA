#' Radar plot of correlations
#'
#' Wrapper using ggplot2 to plot correlations between ICs and given metagenes or
#' pure profiles
#'
#' @param df output of function \code{\link{correlate_metagenes}} - correlation
#'   matrix with correlation and p-values
#' @param ax.size define general size of letters, 15 by default
#' @param size.el.txt  define size of axis labels, adapts automatically by default
#'
#' @return It shows radar plots for each input IC and matagene/profile to
#'   corelate with, it also returns a list containing the data.frame \code{df} used to generate the plot - long
#'   format - and the plot itself \code{p}.
#'
#' @seealso \code{\link[ggplot2]{ggplot}} components,
#'   \code{\link[ggplot2]{geom_point}}, \code{\link[ggplot2]{coord_polar}},
#'   \code{\link[ggplot2]{theme_bw}}, \code{\link[ggplot2]{facet_wrap}},
#'   \code{\link[ggplot2]{scale_color_distiller}}, \code{\link[ggplot2]{theme}},
#'   \code{\link[ggplot2]{element_text}}
#'
#' @export
#'
#' @examples
#'
#' #radar_plot_markers(df, alpha =0.2)
#'
#'
#'
#'
#'
radar_plot_markers <-
  function(df,
           ax.size = NULL,
           size.el.txt = 15) {
    rows <- ceiling(sqrt(ncol(df$r)))
    if (is.null(ax.size)) {
      if (nrow(df$r) < 30L) {
        ax.size <- 8
      } else {
        ax.size <- (2 / nrow(df$r)) * 100
      }
    }
    # Setting the variables to NULL first for R CMD
    IC <- M <- time <- NULL
    chart <-
      data.frame(IC = row.names(df$r[, 1, drop = FALSE]), df$r)
    chart.long  <-
      do.call("rbind", lapply(2:(ncol(chart) - 1), function(i)
        data.frame(chart[, 1], chart[, i], colnames(chart)[i])))
    colnames(chart.long) <- c("IC", "M", "time")
    ICs <- chart.long[order(chart.long[, "IC"]), "IC"]

    (
      p <- ggplot2::ggplot(chart.long, ggplot2::aes(
        x = IC,
        y = M,
        group = time,
        color = M
      )) +
        ggplot2::geom_point(size = 5) +
        ggplot2::coord_polar() +
        ggplot2::theme_bw() +
        ggplot2::facet_wrap( ~ time, nrow = rows, ncol = rows) +
        ggplot2::scale_color_distiller(palette = "Spectral") +
        ggplot2::theme(
          text = ggplot2::element_text(size = size.el.txt),
          axis.text.x = ggplot2::element_text(size = ax.size)
        )
    )
  return(list(df = ICs, p = p))
  }
