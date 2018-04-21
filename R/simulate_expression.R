#' Simulate gene expression
#'
#' Function simulating gene expression of mixed cell types with a perturbator (i.e. proliferation, stress)
#'
#' @param x number of cell types
#' @param n number of genes
#' @param p number of samples
#' @param z number of perturbators
#' @param dist.cells distribution and parameters from which cell profiles will
#'   be drawn
#' @param markers number of markers that will distinguish cell types, can be a
#'   number (the same number of marker genes for cell types and perturbator), can
#'   be a vector of length x+z, it will be set to ceiling(n/20) if not provided
#' @param mfold number of fold change between gene markers and other genes
#' @param CLnames column names (cell and perturbator)
#' @param genes gene names
#' @param dist.noise.sources noise that will be added to each column of basis matrix (to each source)
#' @param alpha parameter for the dirichlet distribution from which are drawn the cell proportions, using rdirichlet.
#' @param dist.noise.global distribution and parameters of global noise (added to each sample one mixture is obtained)
#' @param perturb function of distribution
#' @param pargs arguments of perturbation function
#'
#' @return \describe{
#' \item{expression}{mixed expression matrix}
#' \item{marker.genes}{list of marker genes per cell type}
#' \item{basis_matrix}{pure cell type and perturbator profile}
#' \item{prop}{pure cell type and perturbator proportions (from 0 to 1)}
#' }
#'
#'
#' @export
#'
#' @examples
#'res <- simulate_gene_expresssion (3, 30, 10, 2 , markers = c(4,5,5,3,4))
#'#visualise the basis matrix
#'pheatmap::pheatmap(res$basis_matrix)
#'#visualize expression
#'pheatmap::pheatmap(res$expression)
#'#observe distribution of signals
#'par(mfrow=c(2,2))
#'apply(res$basis_matrix, 2, hist)
simulate_gene_expresssion <-
  function(x,
           n,
           p,
           z = 0,
           dist.cells = list(dist = stats::rnbinom,
                             size = 3,
                             mu = 5),
           markers = NULL,
           mfold = 2,
           CLnames = NULL,
           genes = NULL,
           dist.noise.sources = list(dist = stats::rnorm,
                                     mean = 0,
                                     sd = 0.05),
           alpha = 1,
           dist.noise.global = list(dist = stats::rgamma,
                                    shape = 5,
                                    scale = 1),
           perturb = .pos.gaussian,
           pargs = list(
             p = p,
             mean = 0.5,
             sd = 0.2,
             lwr = 0,
             upr = 1
           )) {
    package <- c("analytics", "NMF")
    new.packages <-
      package[!(package %in% utils::installed.packages()[, "Package"])]
    if (length(new.packages))
      utils::install.packages(new.packages)

    # draw cell types from given distribution
    x.cells <-
      do.call(eval(parse(text = "NMF::rmatrix")), c(list(x = n, y = x), dist.cells))

    # add comumn for perturbation
    if (z > 0) {
      x.cells.pert <-
        cbind(x.cells, NMF::rmatrix(n, z, min = 0, max = 0))
    } else {
      x.cells.pert <- x.cells
    }
    if (is.null(CLnames)) {
      if (z > 0) {
        CLnames <- c(paste0("CL_", 1:x), paste0("P_", 1:z))
      }
      else{
        CLnames <- c(paste0("CL_", 1:x))
      }
    }
    if (is.null(genes)) {
      genes <-
        paste0("gene_", 1:nrow(x.cells.pert))
    }
    colnames(x.cells.pert) <- CLnames
    row.names(x.cells.pert) <- genes

    # define number of markers
    if (is.null(markers)) {
      markers <- ceiling(nrow(x.cells.pert) / 20)
      markers_x <- rep(markers, x)
      markers_z <- rep(markers, z)
    } else if (length(markers) == 1L)  {
      markers_x <- rep(markers, x)
      markers_z <- rep(markers, z)
    } else {
      markers_x <- markers[1:x]
      markers_z <- markers[(x + 1):(x + z)]
    }

    # draw random markers for cell types
    sample.genes_x <- sample(genes, sum(markers_x))
    list.mkgenes_x <-
      split(sample.genes_x, rep(1:length(markers_x), times = markers_x))
    names(list.mkgenes_x) <- CLnames[1:x]
    # draw random markers for paerturbators (they can overlap)
    prob_base <- 1
    prob <- rep(prob_base, length(genes))
    times <- 1.1
    prob[genes %in% sample.genes_x] <- prob_base * times
    if (z > 0) {
      list.mkgenes_z <-
        lapply(markers_z, function(m) {
          sample(genes, m, prob = prob)
        })

      names(list.mkgenes_z) <-
        CLnames[(x + 1):(length(markers_x) + length(markers_z))]

      list.mkgenes <- c(list.mkgenes_x, list.mkgenes_z)
      # build matrix of markers
    } else {
      list.mkgenes <- list.mkgenes_x
    }

    marker.genes.matrix_x <-
      do.call('rbind', lapply(1:length(list.mkgenes_x), function(m) {
        marker <- list.mkgenes_x[[m]]
        res <- apply(x.cells.pert[marker, ], 1, function(gene, i) {
          if (gene[i] >= 10) {
            gene[-c(i)] <-
              stats::runif(length(gene[-c(i)]) , min = 0, max = gene[i] / mfold)
          } else {
            gene[i] <- as.numeric(sample(10:18, 1))
            gene[-i] <-
              stats::runif(length(gene[-c(i)]) , min = 0, max = gene[i] / mfold)
          }
          return(gene)
        }, i = m)
        t(res)
      }))
    if (z > 0) {

      marker.genes.matrix_z <-
        do.call('rbind', lapply(1:length(list.mkgenes_z), function(m) {
          marker <- list.mkgenes_z[[m]]
          new_marker <- marker[!(marker %in% sample.genes_x)]
          if(length(new_marker) ==1L ){
            gene = x.cells.pert[new_marker, ]
            i = length(list.mkgenes_x) + m
            if (gene[i] >= 10) {
              gene[-c(i)] <-
                stats::runif(length(gene[-c(i)]) ,
                             min = 0,
                             max = gene[i] / mfold)
            } else {
              gene[i] <- as.numeric(sample(10:18, 1))
              gene[-i] <-
                stats::runif(length(gene[-c(i)]) ,
                             min = 0,
                             max = gene[i] / mfold)
            }
            res <- t(data.frame(gene))
            row.names(res) <-  new_marker
            return(res)
          } else if (length(new_marker) > 1L) {
          res <-
            apply(x.cells.pert[new_marker, ], 1, function(gene, i) {
              if (gene[i] >= 10) {
                gene[-c(i)] <-
                  stats::runif(length(gene[-c(i)]) ,
                               min = 0,
                               max = gene[i] / mfold)
              } else {
                gene[i] <- as.numeric(sample(10:18, 1))
                gene[-i] <-
                  stats::runif(length(gene[-c(i)]) ,
                               min = 0,
                               max = gene[i] / mfold)
              }
              return(gene)
            }, i = length(list.mkgenes_x) + m)
          t(res)
          }
        }  ))
      for (m in 1:length(list.mkgenes_z)) {
        marker <- list.mkgenes_z[[m]]
        old_marker <- marker[marker %in% sample.genes_x]

        if (length(old_marker) == 1) {
          marker.genes.matrix_x[old_marker, length(list.mkgenes_x) + m] <-
            as.numeric(sample(10:18, 1))
        } else if (length(old_marker) > 1L) {
          i = length(list.mkgenes_x) + m
          marker.genes.matrix_x[old_marker, i] <-
            as.numeric(sample(10:18, 1))
        }
      }

      marker.genes.matrix <-
        rbind(marker.genes.matrix_x, marker.genes.matrix_z)

      marker.genes.matrix <-
        analytics::rowmean(marker.genes.matrix, row.names(marker.genes.matrix))
    } else {
      marker.genes.matrix <- marker.genes.matrix_x
    }

    sample.genes <- row.names(marker.genes.matrix)
    x.cells.pert[sample.genes,] <- marker.genes.matrix
    # add noise
    xN <-
      abs(x.cells.pert +  do.call(eval(parse(text = "NMF::rmatrix")), c(
        list(x = x.cells.pert), dist.noise.sources
      )))
    # simulate proportions of cells
    alpha <- rep(alpha, (ncol(xN) - z))
    prop <- t(gtools::rdirichlet(p, alpha = alpha))
    # simualte proportions of perturbator from normal dist
    if (z > 0) {
      prop.all <-
        rbind(prop, t(sapply(1:z, function(z)
          do.call(perturb, pargs))))
    } else{
      prop.all <- prop
    }

    row.names(prop.all) <- CLnames
    colnames(prop.all) <- paste0("sample_", 1:p)

    ## compose matrix

    y <- xN %*% prop.all

    ### add global noise for each sample
    lambda <-
      do.call(eval(parse(text = "NMF::rmatrix")), c(list(x = nrow(y)), dist.noise.global))
    y <- t(sapply(1:nrow(y), function(i) {
      y[i,] + stats::rnorm(y[i,], mean = 0, sd = 1 / lambda[i])
    }))

    colnames(y) <- paste0("sample_", 1:p)
    row.names(y) <- genes
    return(list(
      expression = y,
      marker.genes  = list.mkgenes,
      basis_matrix = xN,
      prop =  prop.all
    ))
  }
