#' Radar plot of correlations
#'
#' Wrapper using ggplot2 to plot correlations between components and given metagenes or
#' pure profiles
#'
#' @param df output of function \code{\link{correlate_metagenes}} - correlation
#'   matrix with correlation and p-values
#' @param ax.size define size of axis labels, adapts automatically by default
#' @param size.el.txt define general size of letters, 15 by default
#' @param point.size \code{size} parameter in \code{\link[ggplot2]{geom_point}}
#'
#' @return Radar plots for correlations of each input component with matagene/profile,
#' Returns a list containing the \code{data.frame} \code{df} used to generate the plot
#' - long format - and the plot itself \code{p}.
#'
#' @seealso \code{\link[ggplot2]{ggplot}},
#'   \code{\link[ggplot2]{geom_point}}, \code{\link[ggplot2]{coord_polar}},
#'   \code{\link[ggplot2]{theme_bw}}, \code{\link[ggplot2]{facet_wrap}},
#'   \code{\link[ggplot2]{scale_color_distiller}}, \code{\link[ggplot2]{theme}},
#'   \code{\link[ggplot2]{element_text}}
#'
#' @export
#'
#' @examples
#'res_run_ica <- run_fastica (
#'  Example_ds,
#'  overdecompose = FALSE,
#'  n.comp = 20,
#'  with.names = TRUE
#')
#'corr <- correlate_metagenes(
#'    S = res_run_ica$S,
#'    gene.names = res_run_ica$names)
#'
#' radar_plot_corr(corr)
#' data <- radar_plot_corr(corr)$df
#'
#' #change plot
#'radar_plot_corr(corr)$p +
#'ggplot2::labs(title="11 Biton et al. metagenes vs my ICA components",
#'                            subtitle="Pearson correlation coefficients")
#'
#'radar_plot_corr(corr, point.size = 1)
#'
#'radar_plot_corr(corr, point.size = 0)$p +
#'ggplot2::geom_point(size = 8, alpha = 0.4)
radar_plot_corr <-
  function(df,
           ax.size = NULL,
           size.el.txt = 15,
           point.size = 5) {
    package <- "ggplot2"
    new.packages <-
      package[!(package %in% utils::installed.packages()[, "Package"])]
    if (length(new.packages))
      utils::install.packages(new.packages)

    rows <- ceiling(sqrt(ncol(df$r)))
    if (is.null(ax.size)) {
      if (nrow(df$r) < 30L) {
        ax.size <- 8
      } else {
        ax.size <- (2 / nrow(df$r)) * 100
      }
    }
    # Setting the variables to NULL first for R CMD
    component <- correlation <- metagene <- NULL
    chart <-
      data.frame(IC = row.names(df$r[, 1, drop = FALSE]), df$r)

    chart.long  <-
      do.call("rbind", lapply(2:ncol(chart), function(i)
        data.frame(chart[, 1], chart[, i], colnames(chart)[i])))
    colnames(chart.long) <-
      c("component", "correlation", "metagene")
    ICs <- chart.long[order(chart.long[, "component"]), "component"]
    (
      p <- ggplot2::ggplot(
        chart.long,
        ggplot2::aes(
          x = component,
          y = correlation,
          group = metagene,
          color = correlation
        )
      ) +
        ggplot2::geom_point(size = point.size) +
        ggplot2::coord_polar() +
        ggplot2::theme_bw() +
        ggplot2::facet_wrap( ~ metagene, nrow = rows, ncol = rows) +
        ggplot2::scale_color_distiller(palette = "Spectral") +
        ggplot2::theme(
          text = ggplot2::element_text(size = size.el.txt),
          axis.text.x = ggplot2::element_text(size = ax.size)
        )
    )
    return(list(df = chart.long, p = p))
  }
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#'Plot results of density test
#'
#'Wrapper over \code{\link[ggplot2:ggplot]{ggplot}} plotting either rank or
#'density versus selected value in \code{\link{dist_test_samples}} (p.value or
#'test statistics)
#'
#'@param df \code{data.frame} in long format
#'@param plot.type can be either \code{"line"} or \code{"density"}
#'
#'@return
#'returns a line or density plot of p.value or test statistics versus rank or density
#'@export
#'
#' @seealso \code{\link[ggplot2]{ggplot}},
#'   \code{\link[ggplot2]{stat_density}}, \code{\link[ggplot2]{theme_bw}},
#'   \code{\link[ggplot2]{aes}}, \code{\link[ggplot2]{geom_line}}
#' @examples
#'#numerical matrix
#'set.seed(134)
#' S <- matrix(stats::rnbinom(10000, mu = 6, size = 10), 500, 80)
#' dat <- matrix(runif(1600,min =1, max=10 ), 80, 80, byrow = TRUE)
#' A <- dat / rowSums(dat)
#' X <- data.frame(S %*% A)
#' res_run_ica <- run_fastica(X, row.center = TRUE, n.comp = 5, overdecompose = FALSE)
#'
#'#run the funtion selecting wide = FALSE
#' res.ttest <- dist_test_samples(A = res_run_ica$A,
#' sample.names = res_run_ica$samples,
#' X.counts = res_run_ica$log.counts,
#' test.type = "t.test",
#' thr=0.5,
#' isLog = 2,
#' return = "p.value",
#' wide = FALSE)
#'
#'#plot results
#'plot_dist_test(res.ttest, plot.type = "density")
#'plot_dist_test(res.ttest, plot.type = "line")
plot_dist_test <- function(df, plot.type = c("line", "density")) {
  package <- c("ggplot2", "grDevices")
  new.packages <-
    package[!(package %in% utils::installed.packages()[, "Package"])]
  if (length(new.packages))
    utils::install.packages(new.packages)

  value <-  component <-  rank <- NULL
  names(df) <- c("rank", "component", "value")
  grDevices::dev.off()
  switch(
    plot.type,
    density = ggplot2::ggplot(df , ggplot2::aes(x = value, color = component)) +
      ggplot2::stat_density(position = "identity", geom = "line") + ggplot2::theme_bw(),
    line = ggplot2::ggplot(data = df,  ggplot2::aes(
      y = value, x = rank, color = component
    )) +
      ggplot2::geom_line() + ggplot2::theme_bw()
  )
}
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Lolypop plot for correlations
#'
#' Plot correlations between one metagene or known profile and all components in a form of linear plot
#' which is a variant of a signal plot. Wrapper using ggplot2.
#'
#' @details Values are order from highest correlation to lowest correlation. Colors and fonts can be overwritten.
#' To see all correlations simultaneously choose \code{\link{radar_plot_corr}}
#'
#' @param r correlation matrix \code{r} matrix of output \code{\link{correlate_metagenes}}
#' @param col select column either index or column name
#' @param head.size size of the point of correlation
#' @param digits parameter of \code{\link[base]{round}} for the correlation showed on the plot. integer
#' indicating the number of decimal places (round) or significant digits (signif) to be used.
#' @param head.text.size size of the correlation text font
#' @param head.text.color color of the correlation text font
#' @param vertical \code{TRUE} for vertical plot,  \code{FALSE} for horizontal plot
#' @param head.color  by default colored by correlation values, if you want one color provide color name
#'
#' @return
#' returns \code{\link[ggplot2]{ggplot}}
#' @export
#'
#' @seealso \code{\link[ggplot2]{ggplot}},
#'   \code{\link[ggplot2]{aes_string}}, \code{\link[ggplot2]{theme_bw}},
#'   \code{\link[ggplot2]{geom_point}}, \code{\link[ggplot2]{labs}},
#'   \code{\link[ggplot2]{scale_color_distiller}}, \code{\link[ggplot2]{coord_flip}},
#'    \code{\link[ggplot2]{geom_segment}}

#'
#' @examples
#'res_run_ica <- run_fastica (
#'  Example_ds,
#'  overdecompose = FALSE,
#'  n.comp = 20,
#'  with.names = TRUE
#')
#'corr <- correlate_metagenes(
#'    S = res_run_ica$S,
#'    gene.names = res_run_ica$names)
#' #horizontal
#' lolypop_plot_corr(corr$r,2, vertical =FALSE)
#' # vertical
#' lolypop_plot_corr(corr$r,"M8_IMMUNE")
#' #change colors
#' lolypop_plot_corr(corr$r,"M8_IMMUNE",head.color = "black" , head.text.color = "green")
#' #remove title
#' lolypop_plot_corr(corr$r,"M8_IMMUNE")+ ggplot2::labs(title="",subtitle="")
lolypop_plot_corr <-
  function(r,
           col,
           head.size = 10,
           head.color = "value",
           digits = 2,
           head.text.size = 3.5,
           head.text.color = "white",
           vertical = TRUE) {
    package <- c("ggplot2", "grDevices")
    new.packages <-
      package[!(package %in% utils::installed.packages()[, "Package"])]
    if (length(new.packages))
      utils::install.packages(new.packages)

    vec <- r[, col, drop = FALSE]
    r.1 <- data.frame(component = row.names(vec), corr = vec)
    r.2 <- r.1[order(r.1[, 2]), ]
    levels(r.2[, 2]) <- r.2[, 2]
    r.2$component <- factor(r.2$component, levels =  r.2[, 1])
    #input metagene and component (ICs) ordered + level correction
    p <-
      ggplot2::ggplot(
        r.2,
        ggplot2::aes_string(
          x = "component",
          y = colnames(vec),
          label = round(r.2[, colnames(vec)], digits = digits),
          color = colnames(vec)
        )
      )
    if (head.color == "value") {
      p <- p +  ggplot2::geom_point(stat = 'identity', size = head.size)
    } else {
      p <-
        p +  ggplot2::geom_point(stat = 'identity',
                                 size = head.size,
                                 color = head.color)
    }
    p <- p +
      ggplot2::geom_segment(
        ggplot2::aes_string(
          y = min(r.2[, colnames(vec)]) - 0.1,
          x = "component",
          yend = colnames(vec),
          xend = "component"
        ),
        color = "black"
      ) +
      ggplot2::geom_text(color = head.text.color, size = head.text.size) +
      ggplot2::labs(
        title = paste(colnames(vec), " component correlation", sep = ""),
        subtitle = "pearson correlation coefficients",
        y = "Pearson correlation coefficient"
      ) +
      ggplot2::scale_color_distiller(palette = "Spectral") +
      ggplot2::theme_bw()
    if (vertical) {
      (p <- p + ggplot2::coord_flip())
    } else {
      p
    }
    return(p)
  }
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Correlation plot of abundance scores
#'
#' Produces correlation plot of abundance scores estimated versus expected
#'
#' @param x \code{matrix} or \code{data.frame} of abundance scores, samples in rows and cell types in columns
#' @param y \code{matrix} or \code{data.frame} of expected (or to compare) abundance scores, samples in rows and cell types in columns
#' @param ...  additional parameters for \code{method} from \code{\link[corrplot:corrplot]{corrplot}}
#' @details correlation plot between different abundance scores of cell types in samples, correlates both matrices with each other
#' merging two \code{data.frame}s by \code{row.name}s, on \code{\link[corrplot:corrplot]{corrplot}} \code{is.corr} parameter is
#' set to \code{FALSE}
#' @return \describe{\item{correlation plot}{based on \code{\link[corrplot:corrplot]{corrplot}}}
#' \item{\code{corr.full}}{full correlation matrix}
#' \item{\code{corr.filtere}}{correlation without correlation with itself}}
#' @export
#' @seealso  \code{\link[Hmisc:rcorr]{rcorr}}, \code{\link[corrplot:corrplot]{corrplot}}
#' @examples
#' x <- matrix(runif(1000), ncol = 10, nrow = 10)
#' y <- matrix(runif(1000), ncol = 10, nrow = 10)
#' row.names(x) <- row.names(y) <- paste0("S", 1:10)
#' colnames(x) <- paste0("CL_",  1:10, "_estimated")
#' colnames(y) <- paste0("CL_",  1:10, "_expected")
#' scores_corr_plot(x,y, method = "number", tl.col = "black")
#' scores_corr_plot(x,y, method = "square", tl.col = "black")
scores_corr_plot <- function(x, y, ...) {
  package <- c("corrplot", "grDevices")
  new.packages <-
    package[!(package %in% utils::installed.packages()[, "Package"])]
  if (length(new.packages))
    utils::install.packages(new.packages)

  m <- merge(x, y, by = "row.names")
  row.names(m) <- m[, 1]
  m$Row.names <- NULL
  df <- data.frame(m)
  cex.before <- graphics::par("cex")
  graphics::par(cex = 0.7)
  corr <- Hmisc::rcorr(as.matrix(df))
  grDevices::dev.off()
  corrplot::corrplot(
    corr$r[(ncol(x) + 1):nrow(corr$r), 1:ncol(x)],
    is.corr = FALSE,
    tl.cex = 1 / graphics::par("cex"),
    cl.cex = 1 / graphics::par("cex"),
    ...
  )
  return(list(
    corr.full = corr,
    corr.filtered = corr$r[(ncol(x) + 1):nrow(corr$r), 1:ncol(x)]
  ))
}

#
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#' Plot cell proportions
#'
#' Plots scores for all samples as a fraction of one in each samples
#'
#' @param dat scores data.frame with cell types in lines and samples in columns
#'
#' @return
#' a stacked bar plot based on ggplot2
#' @export
#'
#' @examples
#' #random matrix y
#' y <- data.frame(matrix(runif(10000), ncol = 100, nrow = 10))
#' #plot
#' stacked_proportions_plot(y)

stacked_proportions_plot <- function(dat) {
  packages <- c("ggplot2", "reshape")
  new.packages <-
    packages[!(packages %in% utils::installed.packages()[, "Package"])]
  if (length(new.packages))
    utils::install.packages(new.packages)

  cols <- colnames(dat)
  dat <- as.matrix(dat)
  dat <- dat %*% diag(1 / colSums(dat))
  dat <- data.frame(dat)
  colnames(dat) <- cols
  dat$row <- row.names(dat)
  dat2 <- reshape::melt(dat, id.vars = "row")
  colnames(dat2)[1] <- "cell_type"

  dat2[, 1] <- as.factor(dat2[, 1])
  variable <- value <- cell_type <- NULL
  ggplot2::ggplot(dat2, ggplot2::aes(x = variable, y = value, fill = cell_type)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::xlab("\nSample") +
    ggplot2::ylab("Relative proportion\n") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

}
