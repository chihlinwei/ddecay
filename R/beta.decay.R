#' Fit an exponential distance decay model for beta diversity
#'
#' Fit a generalized dissimilarity model to the spatial turnover, nestness or
#' dissmilarity of the ecological assemblages(including zero similarity points)
#' between two sites as a function of their distance apart along an environmental
#' gradient.
#'
#' @param gradient A vector of environmental values of interests
#' @param counts A dataframe of community matrix. Presence/absence if dis.fun="beta.pair"
#' @param coords A dataframe of 2 column with x and y coordinates of the sampling site. The
#' same sites will be remove from bootsrap procedures if like.pairs = TRUE
#' @param nboots A integer number indicating the numbers of bootsrap resampling
#' @param dis.fun Matching function to compute dissimilarity. The default is \code{\link[betapart]{beta.pair}}
#' @param dis An integer of 1, 2 or 3, indicating either the spatial turnover, nestness or total dissmilarity
#' are used in the distance decay modeling.
#' @param like.pairs Logical, whether to remove the like pairs or not. During bootstrap resampling,
#' because the sites are randomyl resampled with replacement and hence there may be multiple copies
#' of some sites. like.pairs = TRUE removes those site pairs having zero separation distance.
#' @param ... Arguments to be passed to the matching function
#' @details This funtion extends \code{\link{dist.decay}} to fit a generalized linear models (GLMs)
#' using quasibinomial distribution (with a log-link function) for the partition or total beta diversity
#' in an ecological assemblage as a function of an environmental gradient. The beta diversity is partitioned
#' into spatial turnover, nestness and totdal dissimilarity using \code{\link[betapart]{beta.pair}},
#' \code{\link[betapart]{bray.part}}, \code{\link[betapart]{functional.beta.pair}}, or \code{\link[betapart]{phylo.beta.pairt}}
#' following the methods decribed by Baselga (2010, 2012). A modified bootstrap method was used to estimate
#' the bootraped mean and standard error of model parameters (alpha, beta), \eqn{s=1-alpha*exp(-beta*d)} for
#' the spatial turnover and dissimilarity and \eqn{s=alpha*exp(-beta*d)} for the nestedness component of betat diversity,
#' where s is assemblage sptail turnover, nestedness or disimilarity and d is environmental distance,
#' and then use these parameters to calculate the turnover, nestedness or disimilarity at zero distance (s0)
#' and the halving distance (halfd) for which the turnover, nestedness or disimilarity between sites
#' decreased/increase by 50\% following the method decribed by Millar et al. (2010).
#' @return A list of data frames, including:
#' \describe{
#'   \item{Results}{A data frame of bootstrap estimate of GLM coefficients, (intercept) and x, as well as beta diversity at zero distance (s0) and halving distance (halfd)}
#'   \item{Summary}{A data frame of bootsrap summary statistics}
#'   \item{Predictions}{A data frame of distance difference, bootstrap predictions and 95\% confidence intervals}
#'   \item{CtrlList}{Other information}
#' }
#' @note This function models the complement of dissimilarity (similarity).
#' @references
#' \itemize{
#'   \item Millar, R.B., Anderson, M.J., Tolimieri, N., 2011. Much ado about nothings: using zero similarity points in distance-decay curves. Ecology 92, 1717-1722. doi:10.1890/11-0029.1
#'   \item Baselga, A. 2010. Partitioning the turnover and nestedness components of beta diversity. Global Ecology and Biogeography 19:134-143
#'   \item Baselga, A. 2012. The relationship between species replacement, dissimilarity derived from nestedness, and nestedness. Global Ecology and Biogeography 21, 1223-1232
#' }
#' @seealso \code{\link{dist.decay}}
#' @author Chih-Lin Wei <chihlinwei@@gmail.com>
#' @export
#' @examples
#' data(os)
#' dd <- beta.decay(gradient=os$dist, counts=decostand(os[, -1:-7], method="pa"),
#'                  coords=os[, c("longitude", "latitude")], nboots=1000,
#'                  dis.fun = "beta.pair", index.family = "sorensen", dis = 1, like.pairs=T)
#' x <- vegdist(os$dist, method = "euclidean")
#' y <- 1 - dist(decostand(os[,-1:-7], method="pa"), index.family = "sorensen")[[1]]
#' plot(x, y)
#' lines(dd$Predictions[, "x"], dd$Predictions[,"mean"], col="red", lwd=2)

beta.decay = function(gradient, counts, coords, nboots = 1000, dis.fun="beta.pair", dis = 1,
                  like.pairs=T, ...) {
  dist = match.fun(dis.fun)
  if(length(gradient) != nrow(as.matrix(counts)))
    stop("\nDifferent number of sites detected in input data")
  nsites = length(gradient)
  Predictions = NULL
  Results = matrix(NA, nrow = nboots, ncol = 2)
  nPairs = rep(NA,nboots)
  for(i in 1:nboots) {
    indices = sample(1:nsites, replace = TRUE)
    ds = as.vector(vegdist(coords[indices,], method = "euclidean"))
    d = as.vector(vegdist(gradient[indices], method = "euclidean"))
    if(dis==2){s = as.vector(dist(counts[indices,], ...)[[dis]])
    }else{s = as.vector(1 - dist(counts[indices,], ...)[[dis]])
    }
    df = data.frame(ds, d, s)
    n = nrow(df)
    #remove like site pairs (zero distance) if requested
    if(!like.pairs) df = subset(df, subset = c(ds > 0))
    n.used = nrow(df)
    if(n.used < 2) stop("Less than 2 site-pairs. Can not fit glm")
    #log-binomial fit
    g = glm(s~d,family = quasibinomial(link=log), data=df)
    Results[i,] = coef(g); nPairs[i] = n.used
    lims = c(0, max(gradient)-min(gradient))
    x = seq(lims[1], lims[2], len = 100)
    Predictions <- cbind(Predictions, predict(g, data.frame("d" = x), type = "response")) }
  #Get parameters of interest
  Results = cbind(Results[,1], Results[,2],
                  ifelse(dis==2, exp(Results[,1]),
                         1-exp(Results[,1])), -log(2)/Results[,2])
  colnames(Results) = c("(intercept)","x","s0","halfd")
  if(!like.pairs) {
    cat("\n",100*(n-mean(nPairs))/n, "% removed due to like pairs")
    cat("\n Bstrap s.e.s have been corrected for the smaller no. of pairs\n\n") }
  Variances = apply(Results, 2, wtd.var, weights = nPairs/n)*mean(nPairs)/n
  Summary = rbind(apply(Results, 2, wtd.mean, weights = nPairs/n),
                  sqrt(Variances),
                  apply(Results, 2, wtd.quantile, weights = nPairs/n, probs = c(.025, .25, .5, .75, .975))
  )
  rownames(Summary) = c("Bstrap mean","Bstrap s.e.", "Bstrap 2.5%", "Bstrap 25%", "Bstrap 50%", "Bstrap 75%", "Bstrap 97.5%")
  if(dis==2){
    ci = t(apply(Predictions, 1, wtd.quantile, weights = nPairs/n, probs = c(.025, .975)))
    mean = apply(Predictions, 1, wtd.mean, weights = nPairs/n)
  }else{
    ci = t(apply(1-Predictions, 1, wtd.quantile, weights = nPairs/n, probs = c(.025, .975)))
    mean = apply(1-Predictions, 1, wtd.mean, weights = nPairs/n)
  }
  return(list(Fits = Results, Summary = Summary, Predictions = cbind(x, mean, ci),
              CtrlList = list(nboots = nboots, like.pairs = like.pairs)))
}
