#' Fit an exponential distance decay model to assemblage similarities
#'
#' Fit a generalized dissimilarity model in ecological similarity of assemblages
#' (including zero similarity points) between two sites as a function of their
#' distance apart along an environmental gradient.
#'
#' @param gradient A vector of environmental values of interests
#' @param counts A dataframe of community matrix of species abundance data
#' @param coords A dataframe of 2 column with x and y coordinates of the sampling site. The
#' same sites will be remove from bootsrap procedures if like.pairs = TRUE
#' @param nboots A integer number indicating the numbers of bootsrap resampling
#' @param dis.fun Matching function to compute dissimilarity. The default is \code{\link[vegan]{vegdist}}
#' @param method Dissimilarity index, partial match to "manhattan", "euclidean", "canberra",
#' "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" ,
#' "binomial", "chao", "cao" or "mahalanobis".
#' @param like.pairs Logical, whether to remove the like pairs or not. During bootstrap resampling,
#' because the sites are randomly resampled with replacement and hence there may be multiple copies
#' of some sites. like.pairs = TRUE removes those site pairs having zero separation distance.
#' @details This funtion fits a generalized linear models (GLMs) using quasibinomial distribution (with
#' a log-link function) for ecological similarity of assemblages between any two sites as a function
#' of their distance apart along a environmental gradient. A modified bootstrap method was used to estimate
#' the bootraped mean and standard error of model parameters (alpha, beta), \eqn{s=alpha*exp(-beta*d)},
#' where s is assemblage similarity and d is environmental distance, and then use these parameters to calculate
#' the similarity at zero distance (s0) and the halving distance (halfd) for which the similarity between sites
#' decreased by 50\% following the method decribed by Millar et al. (2010).
#' @return A list of data frames, including:
#' \describe{
#'   \item{Results}{A data frame of bootstrap estimate of GLM coefficients, (intercept) and x, as well as beta diversity at zero distance (s0) and halving distance (halfd}
#'   \item{Summary}{A data frame of bootsrap summary statistics}
#'   \item{Predictions}{A data frame of distance difference, bootstrap predictions and 95\% confidence intervals}
#'   \item{CtrlList}{Other information}
#' }
#' @note This function models the complement of dissimilarity (similarity).
#' @references
#' \itemize{
#'   \item Millar, R.B., Anderson, M.J., Tolimieri, N., 2011. Much ado about nothings: using zero similarity points in distance-decay curves. Ecology 92, 1717–1722. doi:10.1890/11-0029.1
#' }
#' @seealso \code{\link{beta.decay}}
#' @author Chih-Lin Wei <chihlinwei@@gmail.com>
#' @export
#' @examples
#' # Ostrocode exponential distance decay
#' data(os)
#' dd <- dist.decay(gradient=os$dist, counts=os[, -1:-7], coords=os[, c("longitude", "latitude")], nboots=1000, dis.fun = "vegdist", method = "bray", like.pairs=T)
#' x <- vegdist(os$dist, method = "euclidean")
#' y <- 1-vegdist(os[, -1:-7], method = "bray")
#' plot(x, y)
#' lines(dd$Predictions[, "x"], dd$Predictions[,"mean"], col="red", lwd=2)
#'
dist.decay = function(gradient, counts, coords, nboots = 1000, dis.fun = "vegdist", method = "bray",
                  like.pairs=T) {
  dist = match.fun(dis.fun)
  if(length(gradient) != nrow(as.matrix(counts)))
    stop("\nDifferent number of sites detected in input data")
  nsites = length(gradient)
  Predictions = NULL
  Results = matrix(NA, nrow = nboots, ncol = 2)
  nPairs = rep(NA,nboots)
  for(i in 1:nboots) {
    indices = sample(1:nsites, replace = TRUE)
    ds = as.vector(dist(coords[indices,], method = "euclidean"))
    d = as.vector(dist(gradient[indices], method = "euclidean"))
    s = as.vector(1 - dist(counts[indices,], method = method))
    df = data.frame(ds, d, s)
    n = nrow(df)
    #remove like site pairs (zero distance) if requested
    if(!like.pairs) df = subset(df,subset = c(ds > 0))
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
                  exp(Results[,1]), -log(2)/Results[,2])
  colnames(Results) = c("(intercept)","x","s0","halfd")
  if(!like.pairs) {
    cat("\n",100*(n-mean(nPairs))/n, "% removed due to like pairs")
    cat("\n dist.decay s.e.s have been corrected for the smaller no. of pairs\n\n") }
  Variances = apply(Results, 2, wtd.var, weights = nPairs/n)*mean(nPairs)/n
  Summary = rbind(apply(Results, 2, wtd.mean, weights = nPairs/n),
                  sqrt(Variances),
                  apply(Results, 2, wtd.quantile, weights = nPairs/n, probs = c(.025, .25, .5, .75, .975))
  )
  rownames(Summary) = c("Bstrap mean","Bstrap s.e.", "Bstrap 2.5%", "Bstrap 25%", "Bstrap 50%", "Bstrap 75%", "Bstrap 97.5%")
  ci = t(apply(Predictions, 1, wtd.quantile, weights = nPairs/n, probs = c(.025, .975)))
  mean = apply(Predictions, 1, wtd.mean, weights = nPairs/n)
  return(list(Fits = Results, Summary = Summary, Predictions = cbind(x, mean, ci),
              CtrlList = list(nboots = nboots, like.pairs = like.pairs)))
}
