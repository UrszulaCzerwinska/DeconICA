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
#' @param mfold nulber of folds marker genes should be different from other genes
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
#'res <- simulate_gene_expresssion (3, 100, 100, 2 , markers = 5)
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
           dist.cells = list(dist = stats::rnbinom, size = 3, mu = 5),
           markers = NULL,
           mfold = 2,
           CLnames = NULL,
           genes = NULL,
           dist.noise.sources = list(dist = stats::rnorm, mean = 0, sd = 0.05),
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
    # draw cell types from given distribution
    x.cells <-
      do.call(eval(parse(text = "NMF::rmatrix")), c(list(x = n, y = x), dist.cells))

    # add comumn for perturbation
    if (z > 0) {
      x.cells.pert <-
        cbind(x.cells, NMF::rmatrix(n, z, min = 0, max = 0))
    }
    else {
      x.cells.pert <- x.cells
    }
    if (is.null(CLnames)) {
      CLnames <- paste0("CL_", 1:(x + z))
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
      markers <- rep(markers, x + z)
    } else if (length(markers) == 1L)  {
      markers <- rep(markers, x + z)
    }

    # draw random markers
    sample.genes <- sample(genes, sum(markers))
    list.mkgenes <-
      split(sample.genes, rep(1:length(markers), times = markers))
    names(list.mkgenes) <- CLnames

    # build matrix of markers
    marker.genes.matrix <-
      do.call('rbind', lapply(1:length(list.mkgenes), function(m) {
        marker <- list.mkgenes[[m]]
        res <- apply(x.cells.pert[marker,], 1, function(gene, i) {
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
    x.cells.pert[sample.genes, ] <- marker.genes.matrix
    # add noise
    xN <-
      abs(x.cells.pert +  do.call(eval(parse(text = "NMF::rmatrix")), c(
        list(x = x.cells.pert), dist.noise.sources
      )))
    # simulate proportions of cells
    alpha <- rep(alpha, (ncol(xN) - z))
    prop <- t(gtools::rdirichlet(p, alpha = alpha))
    # simualte proportions of pert from normal dist
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
      y[i, ] + stats::rnorm(y[i, ], mean = 0, sd = 1 / lambda[i])
    }))


    return(list(
      expression = y,
      marker.genes  = list.mkgenes,
      basis_matrix = xN,
      prop =  prop.all
    ))
  }

