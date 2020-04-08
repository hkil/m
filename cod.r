

#===============================================================================================================================

trim <- function(X){
  X <- setNames(X, trimws(names(X)))
  y <- sapply(names(X), function(x) is.character(as.vector(X[[x]])))
  X[y] <- lapply(X[y], trimws)
  return(X)
}


#===============================================================================================================================

reget <- function(List, what){
  
  s <- substitute(what)  
  
  if(class(List)[1] != "list") List <- list(List)  
  
  h <- lapply(List, function(x) do.call("subset", list(x, s)))
  
  res <- Filter(NROW, h)
  
  if(length(res) == 0) NULL else res
}

#===============================================================================================================================

odiag <- function(x) x[(n <- nrow(x))^2-(1:n)*(n-1)]

#===============================================================================================================================

get.uni <- function(data, what){
  
  data$study.name <- trimws(data$study.name)
  m <- split(data, data$study.name)
  m <- Filter(NROW, rm.allrowNA2(m)) 
  
  G <- substitute(what)
  E <- quote(x$x)
  E[[3]] <- G[[2]]
  G[[2]] <- E
  
  f <- sapply(m, function(x) sum(eval(G)) == nrow(x))
  
  h <- m[names(f)[f]]
  
  res <- Filter(NROW, h)
  
  if(length(res) == 0) NULL else res
}


#===============================================================================================================================


get.gen <- function(data, what){
  
  s <- substitute(what)  
  data$study.name <- trimws(data$study.name)
  m <- split(data, data$study.name)
  m <- Filter(NROW, rm.allrowNA2(m)) 
  
  h <- lapply(m, function(x) do.call("subset", list(x, s)))
  
  res <- Filter(NROW, h)
  
  if(length(res) == 0) NULL else res
}

#===============================================================================================================================

find.stud <- function(data, what, timevar = TRUE){
  
  s <- substitute(what)
  data$study.name <- trimws(data$study.name)
  if(!timevar) { unique(as.vector(subset(data, eval(s))$study.name))
    
  } else {
    chep <- sort(unique(na.omit(data$time)))
    G <- lapply(chep, function(x) bquote(.(s) & time == .(x)))
    setNames(lapply(seq_along(G), function(j) unique(as.vector(subset(data, eval(G[[j]]))$study.name))), as.character(chep))
  }
}       

#===============================================================================================================================

find.miss <- function(data, space = FALSE, all = FALSE){    
  
  res <- Filter(length, lapply(data, function(x) which(if(!space & !all) is.na(x) else if(space & !all) x == "" else is.na(x) | x == "")))
  if(length(res) == 0) NA else res
}                    

#===============================================================================================================================

mod.level <- function(data, what){
  
  sort(unique(na.omit(unlist(data[paste0(substitute(what))]))))
  
}

#===============================================================================================================================

mods.level <- function(data){
  
  f <- function(data, what) sort(unique(na.omit(unlist(data[what]))))
  
  ar <- c(formalArgs(d.prepos)[-c(21, 22)], c("dint", "SD", "id"))
  
  dot.names <- names(data)[!names(data) %in% ar]
  
  setNames(lapply(seq_along(dot.names), function(i) f(data = data, what = dot.names[i])), dot.names)
}                  

#===============================================================================================================================


cor.mat <- function(r, dim) { 
  
  m <- matrix(r, dim, dim)
  diag(m) <- 1
  m
}              

#===============================================================================================================================

option1 <- function(ds, sds, r = .5){ 
  
  V <- sds^2  
  
  m <- length(V)
  
  r <- cor.mat(r, m)
  
  SD <- sqrt((1/m)^2 * sum(sqrt(outer(V, V)) * r))
  
  D <- mean(ds)
  
  return(c(D, SD))
}

#===============================================================================================================================

decimal <- function(x, k = 3) format(round(x, k), nsmall = k)              

#===============================================================================================================================

roundi <- function(x, digits = 7){
  
  if(!inherits(x, "data.frame")) stop("Only used for a 'data.frame'.", call. = FALSE)
  
  num <- sapply(x, is.numeric)
  
  x[num] <- lapply(x[num], round, digits)
  
  return(x)
}                

#===============================================================================================================================

cov.dint <- function(sds, r = .5, no.names = TRUE){
  
  m <- length(sds) 
  
  r <- cor.mat(r, m)
  
  D <- diag(sds)
  
  m <- D%*%r%*%D
  
  if(!no.names) rownames(m) <- colnames(m) <- paste0("d", 1:length(sds))
  return(m)
}              


#===============================================================================================================================

option2 <- function(ds, sds, r = .5){
  
  d <- matrix(ds)
  
  if(length(d) == 1) { return(c(d, sds)) }
  
  r <- cor.mat(r, length(d))
  
  e <- matrix(rep(1, length(d)))
  
  A <- cov.dint(sds, r)
  
  w <- t((solve(A)%*%e)%*%solve((t(e)%*%solve(A)%*%e)))
  
  se <- as.vector(sqrt(solve(t(e)%*%solve(A)%*%e)))
  
  ave.d <- as.vector(w%*%d)
  
  return(c(ave.d, se))
}              


#===============================================================================================================================

fuse <- function(..., per.study){
  
  ll <- per.study
  
  L <- list(...)
  
  if(all(sapply(list(...), inherits, "list"))){
    
    g <- lapply(1:length(L), function(i) split(L[[i]], rep(seq_along(ll), ll)))
    
    h <- lapply(1:length(L), function(i) lapply(g[[i]], function(x) do.call(rbind, x)))
    
    lapply(1:length(h), function(i) Filter(NROW, h[[i]]))
    
  } else {
    
    g <- split(L, rep(seq_along(ll), ll))
    
    h <- lapply(g, function(x) do.call(rbind, x))
    
    Filter(NROW, h)
  }
}              


#===============================================================================================================================


d.prepos <- function(d = NA, study.name = NA, group.name = NA, n = NA, mdif = NA, stder = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, r = NA, rev.sign = FALSE, rev.group = FALSE, autoreg = FALSE, t.pair = NA, df = NA, sdif = NA, post = NA, control = NA, outcome = NA, time = NA, ...) 
{
}  

#===================================Modern Inter-rater Reliability in Meta-Analysis=====================================================

rm.allrowNA <- function(X) { 
  
  if(inherits(X, "list")){
    
    lapply(X, function(i) i[rowSums(is.na(i) | i == "") != ncol(i), , drop = FALSE])
    
  } else { X[rowSums(is.na(X) | X == "") != ncol(X), , drop = FALSE] }
}

#===============================================================================================================================

rm.allcolNA <- function(X) { 
  
  if(inherits(X, "list")){
    
    lapply(X, function(i) i[, colSums(is.na(i) | i == "") != nrow(i), drop = FALSE])
    
  } else { X[, colSums(is.na(X) | X == "") != nrow(X), drop = FALSE] }
}

#===============================================================================================================================

rm.colrowNA <- function(X){
  
  r <- rm.allrowNA(X)
  rm.allcolNA(r)  
  
}                                      


#================================================================================================================================

drop.col <- function(dat, vec){
  
  vec <- trimws(vec)
  names(dat) <- trimws(names(dat))
  
  f <- function(dat, vec) {
    i1 <- !names(dat) %in% vec
    setNames(dat[i1], names(dat)[i1])
  }
  
  if(inherits(dat, "list")) { lapply(dat, f, vec = vec)
  } else { f(dat, vec) }
}               

#================================================================================================================================              

full.clean <- function(X, omit, all = TRUE, omit.auto.suffix = TRUE)
{
  
  X <- rm.colrowNA(X)
  
  X <- if(inherits(X, "list") & omit.auto.suffix){ lapply(X, function(x) setNames(x, sub("\\.\\d+$", "", names(x)))) 
    
  } else if(inherits(X, "data.frame") & omit.auto.suffix) { setNames(X, sub("\\.\\d+$", "", names(X))) } else { X }
  
  if(all){ X } else { 
    
    drop.col(X, vec = omit)
  }
}              

#================================================================================================================================                                      

kap <- function (x, level = .95)
{
  
  x <- matrix(table(x[[1]], x[[2]]), 2)
  d  <- diag(x)
  n  <- sum(x)
  nc <- ncol(x)
  colFreqs <- colSums(x)/n
  rowFreqs <- rowSums(x)/n
  
  kappa <- function (po, pc) (po - pc) / (1 - pc)
  
  std  <- function (p, pc, kw, W = diag(1, ncol = nc, nrow = nc)) {
    sqrt((sum(p * sweep(sweep(W, 1, W %*% colSums(p) * (1 - kw)), 2, W %*% rowSums(p) * (1 - kw)) ^ 2) - (kw - pc * (1 - kw)) ^ 2) / crossprod(1 - pc) / n)
  }
  
  po <- sum(d) / n
  pc <- crossprod(colFreqs, rowFreqs)[1]
  k <- kappa(po, pc)
  s <- std(x / n, pc, k)
  
  p <- (1 - level) / 2
  q <- qnorm(c(p, 1-p))
  ci <- k + q*s
  
  return(c(
    KAPPA = k,
    lower = ci[1],
    upper = ci[2],
    conf.level = level))
}                                                                             

#===============================================================================================================================                                      

kappa <- function(X, Y, level = .95, raw.sheet = FALSE){
  
  if(raw.sheet){
    
    ar <- head(formalArgs(d.prepos), -1)
    dot.names <- names(X)[!names(X) %in% ar]
    X <- X[dot.names]
  }
  
  L <- Map(table, X, Y[names(X)])
  
  lapply(L, kap, level = level)
}    


#===============================================================================================================================

detail2 <- function(X, useNA = "ifany"){
  
  nr <- nrow(X)
  nc <- ncol(X)
  tab <- table(row(X), unlist(X), useNA = useNA)
  pj <- apply(tab, 2, sum)/(nr * nc)
  pjk <- (apply(tab^2, 2, sum) - nr * nc * pj)/(nr * nc * (nc - 1) * pj)
  K <- (pjk - pj)/(1 - pj)
  h <- names(K)
  h[is.na(h)] <- "NA"
  setNames(K, h)
}         

#===============================================================================================================================

detail <- function(X, useNA = "ifany") {
  X <- as.matrix(X)
  tab <- table(row(X), unlist(X), useNA = useNA)
  w <- diag(ncol(tab))
  rosum <- rowSums(tab)
  obs_oc <- tab * (t(w %*% t(tab)) - 1)
  obs_c <- colSums(obs_oc)
  max_oc <- tab * (rosum - 1)
  max_c <- colSums(max_oc)
  SA <- obs_c / max_c
  h <- names(SA)
  h[is.na(h)] <- "NA"
  setNames(SA, h)
}                                                       

#===============================================================================================================================                                                          

set.margin <- function() 
{
  par(mgp = c(1.5, 0.14, 0), mar = c(2.5, 2.6, 1.8, .5), 
      tck = -0.02)
}                                                         

#===============================================================================================================================                                                          

splot <- function(y, main, lwd = 5, lend = 2, show.sa = FALSE, digits = 3, cex.sa = .9){
  
  ll <- length(y)
  
  x <- seq_len(ll)
  
  plot(x, y, type = "h", main = main, xlim = c(.95, 1.02*max(x)), ylim = 0:1,
       ylab = "SA%", xaxt = "n", xlab = "Category", lend = lend, lwd = lwd,
       col = colorRampPalette(c(4, 2))(ll), font.lab = 2, 
       panel.first = abline(h = 0, col = 8), las = 1, cex.axis = .9, padj = .3)
  
  if(show.sa) text(x[y != 0]-.015, .4, round(y[y != 0], digits), pos = 2, xpd = NA, srt = 90, font = 2, cex = cex.sa)
  
  axis(1, at = x, labels = names(y))
}


#===============================================================================================================================

irr <- int <- function (X, nsim = 1e3, useNA = "ifany", level = .95, digits = 6, raw = TRUE) 
{
  
  if(!inherits(X, c("data.frame", "matrix", "table"))) stop("Ratings must be 'data.frame', 'matrix', and if not raw, a 'table'.", call. = FALSE)
  
  if(raw) X <- table(row(X), unlist(X), useNA = useNA)
  
  X2 <- X * (X - 1)
  sumcol <- colSums(X)
  sumrow <- rowSums(X)
  nc <- ncol(X)
  nr <- nrow(X)
  tot <- sum(X)
  pij <- X2/(sumrow * (sumrow - 1))
  pi <- rowSums(pij)
  p <- mean(pi)
  pj <- sumcol/tot
  pj2 <- pj^2
  pe <- sum(pj2)
  KAPPA <- (p - pe)/(1 - pe)
  s <- (nc * p - 1)/(nc - 1)
  pi.v.boot <- replicate(nsim, pi.boot <- sample(pi, size = nr, replace = TRUE))
  p.boot <- colMeans(pi.v.boot)
  s.boot <- sapply(seq_len(nsim), function(i) (nc * p.boot[i] - 1)/(nc - 1))
  
  p <- (1 - level) / 2
  s.boot.ci <- quantile(s.boot, probs = c(p, 1-p), na.rm = TRUE)
  
  return(round(c(KAPPA = KAPPA, 
                 Sindex = s, 
                 lower = s.boot.ci[[1]], 
                 upper = s.boot.ci[[2]], 
                 conf.level = level), digits))
}                                      


#===============================================================================================================================                   


int2 <- function(X, level = .95, useNA = "ifany", nsim = 1e3, digits = 4, raw = TRUE){ 
  
  X <- table(row(X), unlist(X), useNA = useNA)
  
  agree.mat <- as.matrix(X) 
  n <- nrow(agree.mat) # number of studies or groups within studies
  q <- ncol(agree.mat) # number of categories
  f <- 0               # population correction 
  
  weights.mat <- diag(q)
  
  agree.mat.w <- t(weights.mat%*%t(agree.mat))
  
  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more
  
  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat) * sum(pi.vec*(1-pi.vec)) / (q*(q-1))
  ac1 <- (pa-pe)/(1-pe)
  
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec == 0) # this replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec
  
  pe.r2 <- pe*(ri.vec>=2)
  ac1.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  pe.ivec <- (sum(weights.mat)/(q*(q-1))) * (agree.mat%*%(1-pi.vec))/ri.vec
  ac1.ivec.x <- ac1.ivec - 2*(1-ac1) * (pe.ivec-pe)/(1-pe)
  
  var.ac1 <- ((1-f)/(n*(n-1))) * sum((ac1.ivec.x - ac1)^2)
  stderr <- sqrt(var.ac1)
  p.value <- 2*(1-pt(ac1/stderr,n-1))
  
  lower <- ac1 - stderr*qt(1-(1-level)/2,n-1)
  upper <- min(1,ac1 + stderr*qt(1-(1-level)/2,n-1))
  
  return(round(c(AC = ac1, lower = lower, upper = upper, conf.level = level), digits))
}                   

#===============================================================================================================================

is.constant <- function(x) length(unique(x)) == 1L 

#===============================================================================================================================                   

drop.inner.list <- function(L, what, omit.auto.suffix = TRUE) {
  
  if(omit.auto.suffix) L <- lapply(L, function(x) setNames(x, sub("\\.\\d+$", "", names(x))))
  
  L[!names(L) %in% what]
}

#===============================================================================================================================

is.unique <- function(X, which){
  
  f <- function(X, which) { nrow(unique(X[which])) == nrow(X[which]) }
  
  test <- if(inherits(X, "list")) sapply(X, f, which) else if(inherits(X, "data.frame")) f(X, which) else {
    
    length(unique(X)) == length(X)
  }
  base::all(test)  
}                                   

#===============================================================================================================================


interrate <- function(..., nsim = 1e3, level = .95, useNA = "ifany", na.rm = FALSE, digits = 3, common = FALSE, all = TRUE, drop = NULL, plot = TRUE, lwd = 5, lend = 1, show.sa = TRUE, group.level = NULL, study.level = NULL, file.name = NULL, reset = TRUE, rev.page = FALSE, cex.sa = .9)
{
  
  r <- list(...) 
  
  if(!all(sapply(r, inherits, c("data.frame", "matrix")))) stop("Coding-sheet(s) must be 'Excel CSV' files, 'data.frame' or 'matrix'.", call. = FALSE)
  
  n.df <- length(r)
  
  r <- lapply(r, as.data.frame)
  
  ar <- formalArgs(d.prepos)[-c(2, 22)]
  
  r <- full.clean(r, ar, all)
  
  r <- lapply(r, function(i) setNames(i, trimws(names(i))))
  
  check <- all(sapply(r, function(i) "study.name" %in% names(i)))
  
  if(!check) stop("Add a new column named 'study.name' to the coding sheet(s).", call. = FALSE)
  
  r <- lapply(r, function(i) {i$study.name <- trimws(i$study.name); i})
  
  r <- lapply(r, function(x) do.call(rbind, c(split(x, x$study.name), make.row.names = FALSE)))
  
  drop <- trimws(drop)              
  drop <- drop[!drop %in% "study.name"]
  
  if(length(drop) != 0) r <- drop.col(r, drop)   
  
  r <- unname(r)
  
  if(n.df == 1) tbl <- table(names(r[[1]])[!names(r[[1]]) %in% c("study.name", "group.name")])
  
  com.names <- if(n.df >= 2) { 
    
    if(common) { Reduce(intersect, lapply(r, names)) 
      
    } else {
      
      vec <- names(unlist(r, recursive = FALSE))
      unique(vec[duplicated(vec)])
    }
    
  } else { 
    
    if(common) { 
      
      names(which(tbl == max(tbl)))
      
    } else {
      
      names(which(tbl >= 2))
    }
  }
  
  dot.names <- if(all) com.names else com.names[!com.names %in% ar]
  
  if(length(dot.names) == 0) stop("No 2 raters detected OR no two moderators names match.", call. = FALSE)
  
  if(n.df >= 2) { 
    
    r <- do.call(cbind, r)
    
    tbl <- table(names(r)[!names(r) %in% c("study.name", "group.name")]) 
    
  } else { r <- r[[1]]
  
  }
  
  n.coder <- if(common) { 
    
    tbl[tbl == max(tbl)] 
    
  } else {
    
    tbl[tbl >= 2]
  }
  
  st.level <- c(names(Filter(base::all, aggregate(.~study.name, r, is.constant, na.action = na.pass)[-1])), if(is.null(study.level)) study.level else trimws(study.level))
  
  st.level <- st.level[st.level %in% dot.names]
  
  exclude <- trimws(group.level)
  
  st.level <- st.level[!st.level %in% exclude]
  
  L <- split.default(r[names(r) %in% dot.names], names(r)[names(r) %in% dot.names])
  
  if(length(st.level) != 0) L[st.level] <- lapply(L[st.level], function(x) x[ave(seq_along(x[[1]]), r$study.name, FUN = seq_along) == 1, ]) 
  
  L <- drop.inner.list(L, c("study.name", "group.name"))
  
  if(na.rm) L <- lapply(L, na.omit)
  
  out <- lapply(L, int, nsim = nsim, level = level, digits = digits, useNA = useNA, raw = TRUE)
  
  A <- lapply(L, detail, useNA = useNA)
  
  study.level <- sapply(seq_along(out), function(i) names(out)[[i]] %in% st.level)
  
  d <- data.frame(out)
  
  d[] <- lapply(d, as.list)
  
  if(plot){
    
    n <- length(L)
    
    if(reset){
      graphics.off()
      org.par <- par(no.readonly = TRUE)
      on.exit(par(org.par))
    }
    dev <- if(!rev.page) n2mfrow(n) else rev(n2mfrow(n))
    if(n > 1L) { par(mfrow = dev) ; set.margin() }
    
    invisible(mapply(splot, y = A, main = names(A), lwd = lwd, lend = lend, show.sa = show.sa, digits = digits, cex.sa = cex.sa))
  }
  
  res <- data.frame(t(rbind(d, row.comprd = sapply(L, nrow), min.cat = sapply(A, function(i) if(any(i < 1)) names(i)[which.min(i)] else "--"), 
                            n.coder = n.coder, study.level = ifelse(study.level, "Yes", "No"))))
  
  file.name <- trimws(file.name)
  
  if(length(file.name) != 0){
    output <- data.frame(lapply(res, unlist))
    nm <- paste0(file.name, ".csv")
    ur <- try(write.csv(output, nm), silent = TRUE)
    if(inherits(ur, "try-error")) stop(paste0("\nClose the Excel file '", nm, "' and try again OR pick another file name."), call. = FALSE)
    message(paste0("\nNote: Check folder '", basename(getwd()),"' for the Excel file '", nm, "'.\n"))
  }
  
  return(res)
}

#===============================================================================================================================                   

do.factor <- function(data, exclude = NULL, char = TRUE, drop = NULL){
  
  colnames(data) <- trimws(colnames(data))
  
  if(!is.null(drop)) data <- drop.col(data, drop)
  
  data <- rm.allrowNA(data) 
  
  ar <- c(formalArgs(d.prepos)[-(20:22)], c("SD", "dint", "id"), exclude)
  
  dot.names <- names(data)[!names(data) %in% ar]
  
  data[dot.names] <- lapply(data[dot.names], if(char) as.character else as.factor)
  
  return(data)
}                  

#===============================================================================================================================

tplot <- function(y, main, lwd = 4, lend = 2, cat.level = 0, low = FALSE){
  
  if(!low) low <- NULL
  
  z <- length(y) 
  x <- seq_len(z)
  
  if(cat.level != 0 & z >= cat.level) { main <- bquote(bold(.(main)~symbol(("\326")))) ; col.main <- "magenta"} else { main <- main ; col.main <- 1} 
  plot(x, y, type = "h", main = main, xlim = c(.95, 1.02*max(x)),
       ylab = "Frequency", axes = FALSE, xlab = "Category", lwd = lwd,
       col = colorRampPalette(c(4, 2))(z), font.lab = 2, lend = lend, col.main = col.main)
  box()
  axis(1, at = which(!names(y) %in% names(low)), labels = names(y)[!names(y) %in% names(low)], cex.axis = .9)
  if(!is.null(low)) axis(1, at = which(names(y) %in% names(low)), labels = names(low), cex.axis = .9, col = "magenta", col.axis = "magenta", font = 2)
  if(!is.null(low)) text(which(names(y) %in% names(low)), low, low, font = 2, cex = .9, pos = 3, col = "magenta")
  axis(2, at = pretty(y), cex.axis = .85, las = 1, padj = .3)
}

#================================================================================================================================================================

plot.mods <- function(data, exclude = NULL, lwd = 4, lend = 2, cat.level = 0, code = NULL, low = NULL){
  
  lo <- low
  low <- if(is.null(low)) FALSE else TRUE
  
  names(data) <- trimws(names(data))
  
  data <- rm.colrowNA(data) 
  
  ar <- c(formalArgs(d.prepos)[-(20:22)], c("SD", "dint", "id", "study.name"), exclude)
  
  mods <- names(data)[!names(data) %in% ar]
  
  A <- setNames(lapply(seq_along(mods), function(i) table(data[[mods[i]]], dnn = NULL)), mods)
  Ls <- lapply(A, length)
  
  bad <- Ls < 2
  bad.names <- names(A[bad])
  A <- A[!bad]
  
  A <- A[Ls >= cat.level]
  if(length(A) == 0) stop(paste("No variable with cat.level >=", if(cat.level != 0) cat.level else 2, "found."), call. = FALSE)  
  
  if(!is.null(code)){
    target <- sapply(seq_along(A), function(i) any(names(A[[i]]) == code))
    A <- A[target]
  }
  
  if(low){
    target <- sapply(seq_along(A), function(i) any(A[[i]] <= lo))
    A <- A[target]
    low <- if(length(A) == 0) FALSE else setNames(lapply(seq_along(A), function(i) A[[i]][which(A[[i]] <= lo)]), names(A))
  }
  
  n <- length(A)
  graphics.off()
  org.par <- par(no.readonly = TRUE)
  on.exit(par(org.par))
  if(n > 1L) { par(mfrow = n2mfrow(n)) ; set.margin() }
  
  if(length(bad.names) != 0) message("Note: ", toString(dQuote(bad.names)), "ignored due to insufficient category levels.")
  if(length(A) > 0) { invisible(mapply(tplot, y = A, main = names(A), lwd = lwd, lend = lend, cat.level = cat.level, low = low)) ; return(A) } else NA
}


#================================================================================================================================================================

mix.level <- function(data, exclude = NULL){
  
  ar <- c(formalArgs(d.prepos)[-c(20:22)], c("SD", "dint", "id"), exclude)  
  
  mods <- names(data)[!names(data) %in% ar]
  
  tmp <- do.call(rbind, lapply(mods, function(x){
    d <- setNames(unique(data[c("study.name", x)]), c("study.name", "code"))
    transform(d, mod.name = x)
  }))
  
  res <- tmp[with(tmp, ave(code, code, mod.name, FUN = length) == 1),]
  if(nrow(res) == 0) return(NULL)
  res <- res[order(res$study.name),]
  rownames(res) <- NULL
  res
}                           


#================================================================================================================================================================

group.level <- function(data, exclude = NULL){
  
  ar <- c(formalArgs(d.prepos)[-c(2,20:22)], c("SD", "dint", "id"), exclude)
  
  d <- drop.col(data, ar)
  
  molten <- data.frame(d[, 1, drop = FALSE], stack(d[, -1]))
  
  res <- molten[as.logical(ave(molten[['values']], molten[['ind']], 
                               FUN = function(x) !duplicated(x) & !duplicated(x, fromLast = TRUE))), ]
  
  if(nrow(res) == 0) return(NULL)
  res <- setNames(res[order(res$study.name),], c("study.name", "code", "mod.name"))
  rownames(res) <- NULL
  res
}                               

#================================================================================================================================================================

study.level <- function(data, exclude = NULL){
  
  mix <- mix.level(data = data, exclude = exclude)
  grp <- group.level(data = data, exclude = exclude)
  
  if(!is.null(grp) & !is.null(mix)) { 
    
    res <- subset(mix, !((study.name %in% grp$study.name) & (code %in% grp$code) & (mod.name %in% grp$mod.name)))
    res <- res[order(res$study.name),]    
    rownames(res) <- NULL
    res
  } else if(is.null(grp) & !is.null(mix)){ mix } else { NULL }
}                                                                              

#================================================================================================================================================================

print.labdf <- function(data) {
  dn <- dimnames(data)
  names(dn) <- attr(data, "rclab")
  data <- as.matrix(data)
  dimnames(data) <- dn
  print(data, quote = FALSE)
}                               

#================================================================================================================================================================

suggest <- function(data, exclude = NULL){
  
  ar <- c(formalArgs(d.prepos)[-c(2,20:22)], c("SD", "dint", "id"), exclude)
  data <- drop.col(data, ar)
  
  DF <- data[with(data, ave(as.numeric(study.name), study.name, FUN = length)) >= 3, ]
  long.df <- cbind(DF[1], stack(DF[-1]))
  res <- long.df[with(long.df, ave(values, study.name, ind, values, FUN = length) == 1 &
                        ave(values, study.name, ind, FUN = function(x) length(unique(x))) == 2), ]
  
  if(nrow(res) == 0) return(NULL)
  res <- setNames(res[order(res$study.name),], c("study.name", "code", "mod.name"))   
  rownames(res) <- NULL
  res
}    

#================================================================================================================================================================

rule3 <- function(data, exclude = NULL, low = 4){
  
  ar <- c(formalArgs(d.prepos)[-c(20:22)], c("SD", "dint", "id"), exclude)  
  
  mods <- names(data)[!names(data) %in% ar]
  
  A <- setNames(lapply(seq_along(mods), function(i) table(data[[mods[i]]], dnn = NULL)), mods)
  
  target <- sapply(seq_along(A), function(i) any(A[[i]] <= low))
  A <- A[target]
  
  if(length(A) == 0) return(NULL)
  
  lo <- setNames(lapply(seq_along(A), function(i) A[[i]][which(A[[i]] <= low)]), names(A))
  
  lst <- Filter(length, lapply(split(data[names(lo)], data$study.name), 
                               function(dat) Filter(nrow, Map(function(x, y) 
                                 merge(x, y[setdiff(names(y), "values")], by = "ind"), lapply(dat, 
                                                                                              function(x) stack(table(x))), lapply(lo, stack)))))
  
  res <- do.call(rbind, c(Map(cbind, study.name = names(lst), lapply(lst, 
                                                                     function(x) do.call(rbind, c(Map(cbind, x, mod.name = names(x)),
                                                                                                  make.row.names = FALSE)))), make.row.names = FALSE))
  
  r3 <- setNames(res, c("study.name","code","occurs","mod.name"))[, c(1,2,4,3)]
  r3 <- r3[order(r3$mod.name),]    
  rownames(r3) <- NULL
  
  grp <- group.level(data = data, exclude = exclude)
  
  if(is.null(grp)){ r3 } else {
    res <- subset(r3, !((study.name %in% grp$study.name) & (code %in% grp$code) & (mod.name %in% grp$mod.name)))
    res <- res[order(res$mod.name),]    
    rownames(res) <- NULL
    res
  } 
}                                     

#================================================================================================================================================================

exam.code <- function(data, exclude = NULL, rule = 1, lwd = 4, lend = 2, cat.level = 0, code = NULL, low = NULL, suggest = FALSE, plot = TRUE){
  
  names(data) <- trimws(names(data))
  check <- "study.name" %in% names(data)
  if(!check) stop("Add a new column named 'study.name'.", call. = FALSE)
  
  data$study.name <- trimws(data$study.name)
  data <- rm.colrowNA(data)
  
  if(length(unique(data$study.name)) < 2) stop("At least two coded studies required.", call. = FALSE)
  
  exclude <- trimws(exclude)  
  excl <- setdiff(exclude, "study.name")
  
  exclude <- if(!is.null(excl) & length(excl) != 0) exclude else NULL
  
  if(plot) invisible(plot.mods(data = data, exclude = exclude, lwd = lwd, lend = lend, cat.level = cat.level, code = code, low = low))
  
  h <- if(rule == 1 & !suggest) study.level(data = data, exclude = exclude) else
    if(rule == 2 & !suggest) group.level(data = data, exclude = exclude) else 
      if(rule == 3 & !suggest) rule3(data = data, exclude = exclude, low = if(is.null(low)) 4 else low) else
        if(suggest) suggest(data = data, exclude = exclude)
  
  if(!is.null(h)){    
    
    attr(h, "rclab") <- c("", paste0(if(suggest) "Possible Mistakes" else paste("Violations of Rule", rule), ":"))
    class(h) <- c("labdf", class(h)) 
    h } else { h }
}                                                                           


#================================================================================================================================================================


impute <- function(D, FUN = median){
  y <- sapply(D, is.numeric)
  D[y] <- lapply(D[y], function(x){x[is.na(x)] <- FUN(x, na.rm = TRUE); x})
  return(D)
}                                   


#================================================================================================================================================================


plan.item <- function(margin = .2, S.index = NA){
  
  margin[margin > .99] <- .99  
  margin[margin < .01] <- .01
  
  data.frame(n.item = ceiling(1/margin^2), margin = paste0("+-", margin), S.index = S.index, lower = if(!is.na(S.index)) S.index - margin else NA, upper = if(!is.na(S.index)) S.index + margin else NA)
}

#================================================================================================================================================================

plan.coder <- function(se2irr = .2, S.index = NA){
  
  se2irr[se2irr > .99] <- .99  
  se2irr[se2irr < .01] <- .01 
  
  n <- ceiling(2/se2irr)
  se <- if(!is.na(S.index)) se2irr*S.index else NA
  
  data.frame(n.coder = n, se2irr = se2irr, S.index = S.index, lower = if(!is.na(S.index)) S.index - 2*se else NA, upper = if(!is.na(S.index)) S.index + 2*se else NA, conf.lev = .95)
}

#================================================================================================================================================================

irr.diag <- function(X, useNA = "ifany"){
  
  a <- detail2(X, useNA = useNA)
  b <- detail(X, useNA = useNA)
  
  round(data.frame(KAPPA = a, SA = b), 3)
}

#==========================================================================================================================================

find.irr <- function(X, what){
  
  if(!inherits(X, "data.frame")) stop("Data must be an Excel CSV file or a 'data.frame'.", call. = FALSE)
  
  s <- as.list(substitute(what))  
  
  res <- Filter(NROW, X[rowSums(X[grep(as.character(s[[2]]), names(X))] == s[[3]], na.rm = TRUE) > 0,][c("study.name", "group.name")])
  
  if(length(res) == 0) NULL else res
}           

#===========================# Datasets # ===================================================================================== 

table1 <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/irr1.csv", row.names = 1)
table2 <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/irr2.csv", row.names = 1)          
table3 <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/X.csv", row.names = 1)
c1 <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/c1.csv")
c2 <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/c2.csv")
c3 <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/c3.csv")
c4 <- read.csv("https://raw.githubusercontent.com/rnorouzian/m/master/c4.csv")             

#================================================================================================================================================================


options(warn = -1)
                                      
