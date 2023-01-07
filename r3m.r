# functions for 'formetafor'

# H===============================================================================================================================

trim_ <- function(X){
  X <- setNames(X, trimws(names(X)))
  y <- sapply(names(X), function(x) is.character(as.vector(X[[x]])))
  X[y] <- lapply(X[y], trimws)
  return(X)
}

# H===============================================================================================================================

abber_case <- function(X, abb_except = NULL, abb_length = 10){
  y <- names(Filter(function(i) is.character(i) | is.factor(i), X[setdiff(names(X), abb_except)]))
  X[y] <- lapply(X[y], abbreviate, minlength=abb_length, named=FALSE)
  return(X)
}                    

# H===============================================================================================================================                    

numerize_case <- function(X, num_except = NULL, num_zero = FALSE){
  y <- names(Filter(function(i) is.character(i) | is.factor(i), X[setdiff(names(X), num_except)]))
  X[y] <- lapply(X[y], function(k) as.integer(as.factor(k))-if(num_zero) 1 else 0)
  return(X)
} 

# H================================================================================================================================

rm.allrowNA <- function(X) { 
  
  X[rowSums(is.na(X) | X == "") != ncol(X), , drop = FALSE]  
  
}

# H===============================================================================================================================

rm.allcolNA <- function(X) { 
  
  X[, colSums(is.na(X) | X == "") != nrow(X), drop = FALSE]
  
}

# H===============================================================================================================================

rm.colrowNA <- function(X){
  
  rm.allcolNA(rm.allrowNA(X))  
  
}

# H===============================================================================================================================

get_vars_ <- function(gls_fit, as_fml = TRUE){ 
  
  m <- sub("[^:]+\\(([^)]+).*", "\\1",attr(terms(gls_fit), "term.labels"))
  
  if(as_fml) sapply(paste0("~",m), as.formula) else m
}

# H===============================================================================================================================

get_data_ <- function (object) 
{
  if ("data" %in% names(object)) {
    data <- object$data
  }
  else {
    dat_name <- object$call$data
    envir_names <- sys.frames()
    ind <- sapply(envir_names, function(e) exists(as.character(dat_name), 
                                                  envir = e))
    e <- envir_names[[min(which(ind))]]
    data <- eval(dat_name, envir = e)
  }
  if (is.null(data)) 
    return(data)
  naAct <- object[["na.action"]]
  if (!is.null(naAct)) {
    data <- if (inherits(naAct, "omit")) {
      data[-naAct, ]
    }
    else if (inherits(naAct, "exclude")) {
      data
    }
    else eval(object$call$na.action)(data)
  }
  subset <- object$call$subset
  if (!is.null(subset)) {
    subset <- eval(asOneSidedFormula(subset)[[2]], data)
    data <- data[subset, ]
  }
  data
}                 

# H=============================================================================================================================== 

odds_. <- function(x) subset(x, x %% 2 != 0)

# H===============================================================================================================================

print.post_rma <- function(post_rma_call){
  print(post_rma_call$table)
}                 

# H=============================================================================================================================== 

full_clean <- function(data) rm.colrowNA(trim_(data))           

# H===============================================================================================================================

odiag <- function(x) suppressMessages(x[col(x) != row(x)])                   

# H===============================================================================================================================

shift_rows <- function(data, user_index, up = TRUE){
  
  indx <- seq_len(nrow(data))
  remain <- indx[!indx %in% user_index]
  
  data[if(up) c(user_index, remain) else c(remain, user_index), ]
}                    

# H===============================================================================================================================

get_error_rho <- function(fit){
  
  is_V <- any(odiag(fit$V) != 0)  
  
  if(is_V) {
    
    u <- unique(odiag(round(cov2cor(fit$V), 7)))
    u[u!=0]
    
  } else 0
  
}                    

# H===============================================================================================================================                    

is_crossed <- function(obj){
  
  mod_struct <- clubSandwich:::parse_structure(obj)
  highest_cluster <- names(mod_struct$level_dat)[which.min(mod_struct$level_dat)]
  cluster <- mod_struct$cluster_dat[[highest_cluster]]
  out <- !clubSandwich:::test_nested(cluster, fac = mod_struct$cluster_dat)
  out[names(out) %in% obj$s.names]
}

# H===============================================================================================================================

pluralify_ <- function (x, keep.original = FALSE, 
                        irregular = lexicon::pos_df_irregular_nouns) {
  
  stopifnot(is.data.frame(irregular))
  
  hits <- match(tolower(x), tolower(irregular[[1]]))
  
  ends <- "(sh?|x|z|ch)$"
  plural_ <- ifelse(grepl(ends, x), "es", "s")
  out <- gsub("ys$", "ies", paste0(x, plural_))
  out[which(!is.na(hits))] <- irregular[[2]][hits[which(!is.na(hits))]]
  
  c(if (keep.original) {
    x
  }, out)
}           

# M================================================================================================================================     

cat_overlap <- function(data, study_id, ..., blank_sign = "-"){
  
  data <- full_clean(data)
  study_id <- rlang::ensym(study_id)
  cat_mod <- rlang::ensyms(...)
  cat_nms <- purrr::map_chr(cat_mod, rlang::as_string)
  
  idx <- cat_nms %in% names(data)
  if(!all(idx)) stop(toString(dQuote(cat_nms[!idx]))," not found in the 'data'.", call. = FALSE)
  
  setNames(purrr::map(cat_mod,  ~ {
    
    studies_cats <- 
      data %>%
      dplyr::group_by(!!study_id, !!.x) %>%
      dplyr::summarise(effects = n(), .groups = 'drop')
    nm1 <- rlang::as_string(.x)
    cat_names <- paste0(nm1, c(".x", ".y"))
    
    studies_cats <- 
      studies_cats %>%
      dplyr::inner_join(studies_cats, by = rlang::as_string(study_id)) %>%
      dplyr::group_by(!!!rlang::syms(cat_names)) %>%
      dplyr::summarise(
        studies = n(),
        effects = sum(effects.x), .groups = 'drop') %>% 
      dplyr::mutate(n = paste0(studies, " (", effects, ")") )
    
    out1 <- studies_cats %>%
      dplyr::select(-studies, -effects) %>%        
      tidyr::pivot_wider(names_from = cat_names[2], 
                         values_from = n, names_sort = TRUE) %>%
      dplyr::rename_with(~nm1,  cat_names[1]) %>%
      dplyr::arrange(dplyr::across(tidyselect::all_of(nm1))) %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(nm1), as.character))  
    
    out2 <- out1[-1]
    out2[upper.tri(out2)] <- blank_sign
    dplyr::bind_cols(out1[1],out2)
    
  }), cat_nms)
}           


# H================================================================================================================================ 

data.tree_ <- function(data, toplab = NULL, cex = 1, rowcount = FALSE, cex_top = 1, ...){
  
  toplab <- if(is.null(toplab)) names(data) else toplab
  
  sizetree2(data, toplab = toplab, stacklabels = FALSE, border = 0, base.cex = cex, showcount = rowcount, cex_top = cex_top, ...)
}           

# M===============================================================================================================================

meta_tree <- function(data, ..., effect = TRUE, highest_level_name = NULL,  
                      structure = c("simple","typical","complex"), toplab = NULL, 
                      main = NULL, main_extra = NULL, rowcount = FALSE, 
                      abb_names = FALSE, abb_length = 6, abb_except = NULL, 
                      num_names = FALSE, num_except = NULL, num_zero = FALSE, 
                      panel_label = TRUE, cex = 1, cex_main = 1, rev_order = TRUE, 
                      rev_page = FALSE, reset = TRUE, index = NULL, cex_top = 1) 
{
  
  data <- full_clean(data) %>%
    mutate(effect = row_number())
  
  dot_cols <- rlang::ensyms(...)
  dot_cols <- if(effect) append(dot_cols, rlang::sym("effect")) else dot_cols
  str_cols <- purrr::map_chr(dot_cols, rlang::as_string)
  
  ss <- dot_cols[[1]]
  sss <- str_cols[1]
  
  idx <- str_cols %in% names(data)
  if(!all(idx)) return(message("Error: ",toString(dQuote(str_cols[!idx]))," not found in the data."))
  
  main_org <- main
  
  if(!is.null(highest_level_name) & abb_names) abb_except <- c(sss, abb_except)
  if(!is.null(highest_level_name) & num_names) num_except <- c(sss, num_except)
  
  if(abb_names) data <- abber_case(data, abb_length = abb_length, abb_except = abb_except) 
  if(num_names) data <- numerize_case(data, num_except = num_except, num_zero = num_zero)
  
  if(reset){
    graphics.off()
    org.par <- par(no.readonly = TRUE)
    on.exit(par(org.par))
  }
  
  data <- data %>%
    dplyr::select(!!!dot_cols)
  
  if(is.null(highest_level_name)){
    
    struc <- match.arg(structure) 
    
    hlist <- data %>%
      dplyr::group_by(!!sym(sss)) %>%
      dplyr::mutate(grp = dplyr::across(tidyselect::all_of(str_cols[-1]), ~ {
        tmp <- dplyr::n_distinct(.)
        dplyr::case_when(tmp == 1 ~ 1, tmp == n() ~ 2, tmp > 1 & tmp < n() ~ 3,  TRUE ~ 4)
      }) %>%
        purrr::reduce(stringr::str_c, collapse = "")) %>%
      dplyr::ungroup(.) %>%
      dplyr::group_split(grp, .keep = FALSE)
    
    hlist <- if(rev_order) rev(hlist) else hlist
    
    res <- Filter(NROW, hlist)
    
    LL <- length(res)
    
    if(!is.null(index)){ 
      
      index <- index[index <= LL & index >= 1]
      if(length(index)==0) index <- NULL
    } 
    
    res <- if(!is.null(index)) res[index] else res
    
    main_no. <- sapply(res, function(i) length(unique(i[[sss]])))
    
    typic <- function(vec) vec[ceiling(length(vec)/2)]
    
    nms <- lapply(res, function(i){
      nr <- sapply(split(i, i[[sss]]), nrow);
      study_type <- if(struc == "typical") {typic(as.numeric(names(table(nr))))
      } else if(struc == "simple") {min(as.numeric(names(table(nr))))
      } else {max(as.numeric(names(table(nr))))};
      names(nr)[nr == study_type][1]
    })
    
    list2plot <- lapply(seq_along(res),function(i) filter(res[[i]], !!ss == nms[i]))
    
    LL <- length(list2plot)
    
    if(LL > 1L) { 
      
      dev <- if(!rev_page) n2mfrow(LL) else rev(n2mfrow(LL))
      par(mfrow = dev) 
      
    }
    
    main <- if(is.null(main)) stringr::str_to_title(ifelse(main_no. > 1, pluralify_(sss), sss)) else main
    
    main <- if(is.null(main_org)) paste(main_no., main) else main
    
    if(panel_label) {
      
      pan_lab <- make.unique(rep(LETTERS,1e1))[seq_along(list2plot)]
      
      main <- paste0("(",pan_lab,") ", main)
    }
    
    if(!is.null(main_extra)) main <- paste0(main, " [",main_extra,"]")
    
    invisible(lapply(seq_along(list2plot), function(i) data.tree_(list2plot[[i]], main = main[i], toplab, cex, rowcount, cex.main = cex_main, cex_top = cex_top)))
    
    invisible(if(panel_label) setNames(res, pan_lab) else res)
    
  } else {
    
    highest_level_name <- trimws(highest_level_name)
    highest_level_names <- unique(data[[sss]])
    
    idx <- highest_level_name %in% highest_level_names 
    
    if(!all(idx)) return(message("Error: ",toString(dQuote(highest_level_name[!idx]))," not found in the ", paste0("'",sss,"'", " data column.")))
    
    list2plot <- lapply(highest_level_name, function(i) filter(data, !!ss == i))
    
    LL <- length(list2plot)
    
    if(!is.null(index)){ 
      
      index <- index[index <= LL & index >= 1]
      if(length(index)==0) index <- NULL
    } 
    
    list2plot <- if(!is.null(index)) list2plot[index] else list2plot
    
    LL <- length(list2plot)
    
    if(LL > 1L) { 
      
      dev <- if(!rev_page) n2mfrow(LL) else rev(n2mfrow(LL))
      par(mfrow = dev) 
      
    }
    
    invisible(lapply(list2plot, data.tree_, toplab, cex, rowcount, cex.main = cex_main, main = main, cex_top = cex_top))
  }
}

# M================================================================================================================================================

interactive_outlier <- function(fit, cook = NULL, st_del_res_z = NULL, 
                                cex_add_point = .5, cex_multi_point = 1.2,
                                whisker_coef = 2.5, cex_text_outlier = .6,
                                cex_main = .9, parallel = "no", ncpus = 1, 
                                reestimate = FALSE, save = FALSE, 
                                file_name_cook = "cooks1",
                                file_name_res_z = "rstudent1",
                                view = 1, pos = 2)
{
  
  if(!inherits(fit,c("rma.mv"))) stop("Model is not 'rma.mv'.", call. = FALSE)
  datziola <- get_data_(fit) %>%
    mutate(obsss = fit$slab)
  
  hat <- hatvalues.rma.mv(fit)
  
  if(is.null(cook)){
    
    cook <- cooks.distance.rma.mv(fit, progbar=TRUE,
                                  parallel = parallel, 
                                  ncpus = ncpus, 
                                  reestimate = reestimate)
    
    if(save){
      
      filenm <- paste0(file_name_cook,".rds")
      saveRDS(cook, filenm)
      
      message("\nNote: Check folder '", basename(getwd()),"' for the ", dQuote(filenm)," file.\n") 
    }
  }
  
  if(is.null(st_del_res_z)){
    
    st_del_res_z <- rstudent.rma.mv(fit, progbar=TRUE,
                                    parallel = parallel, 
                                    ncpus = ncpus, 
                                    reestimate = reestimate)$z
    
    if(save){
      
      filenm <- paste0(file_name_res_z,".rds")
      saveRDS(st_del_res_z, filenm)
      
      message("\nNote: Check folder '", basename(getwd()),"' for the ", dQuote(filenm)," file.\n") 
    }
  }
  
  # Make visual size of effects proportional to their hat/cook's distance (estimate influence)
  cex <- cex_add_point+cex_multi_point*sqrt(if(view == 1)cook else hat)
  
  outlier_limits <- qnorm(c(.025,.5,.975))
  ylim <- range(c(outlier_limits, st_del_res_z))
  
  # Plot Leverage against Studentized residuals proportioned on cook's distances
  plot(if(view == 1) hat else cook, st_del_res_z, cex=cex, las=1, mgp=c(1.5,.3,0),
       xlab = if(view == 1) "Leverage (Hat Value)" else "Effect Influence (Cook's Dis.)", 
       ylab = "Outlier (Standardized Del. Value)",pch=19,cex.axis = .9,tcl = -.3,
       col = adjustcolor(1, .5),
       ylim = ylim)
  
  title(if(view == 1) "Size of points denote \nestimate-influencing effects\n (Cook's distances)" 
        else "Size of points denote \nleverage effects\n (Hat value)", 
        cex.main = cex_main, line = .3)
  
  abline(h=outlier_limits, lty=c(3,1,3), lwd=c(1,2,1))
  
  max_hat <- max(mean(range(hat)), boxplot.stats(hat, coef = whisker_coef)$stats[5])
  
  max_cook <- max(mean(range(cook)), boxplot.stats(cook, coef = whisker_coef)$stats[5])
  
  abline(v = if(view == 1) max_hat else max_cook, col=2)
  
  # To be outlier, an estimate must simultaneously (a) be outlying (per studentized value)
  # (b) have high leverage (per hat value), and (c) high model influence (per cook's distance)
  
  i <- abs(st_del_res_z) > outlier_limits[3]  
  j <- hat > max_hat
  k <- cook > max_cook
  L <- which(i & j & k)
  
  if(length(L)==0) { 
    
    message("Note: No interactive outlier detected.") 
    
    return(NA)
  }
  
  u <- par()$usr
  
  if(any(st_del_res_z[L]>0)) rect(if(view == 1) max_hat else max_cook, outlier_limits[3], u[2], u[4], col = adjustcolor(2, .2), border = NA)
  if(any(st_del_res_z[L]<0)) rect(if(view == 1) max_hat else max_cook, outlier_limits[1], u[2], u[3], col = adjustcolor(2, .2), border = NA)
  
  # Show which effects meet all the three conditions
  text(if(view == 1) hat[L] else cook[L], st_del_res_z[L], labels = names(L), pos = pos, col = "magenta", cex = cex_text_outlier, xpd = NA)
  points(if(view == 1) hat[L] else cook[L], st_del_res_z[L], cex=cex[L])
  
  LL <- names(L)
  
  removed <- filter(datziola, obsss %in% LL)
  new_data <- filter(datziola, !obsss %in% LL)
  
  return(list(removed = removed, 
              new_data = new_data))
}   

# H=================================================================================================================================================

fixed_form_rma <- function(fit){ 
  
  a <- fit$formula.yi
  b <- fit$formula.mods
  y <- fit$call$yi
  
  if(!is.null(a) & !is.null(b)) a 
  else if(!is.null(a) & is.null(b)) a
  else if(is.null(a) & !is.null(b)) as.formula(paste(as.character(y), paste(as.character(b), collapse = "")))
  else as.formula(paste(as.character(y), "~ 1"))
}

# H=================================================================================================================================================  

shorten_ <- function(vec, n = 3) { 
  
  gsub("\\s+", "", 
       sub(sprintf("^(([^,]+, ){%s}).*, ([^,]+)$", n), "\\1...,\\3", toString(vec)))  
} 

# M======================================================================================================================================================  

plot_rma <- function(fit, full=TRUE, multiline=TRUE, 
                     confint = list(style="auto"), x.var, 
                     key.args= list(space="top",cex=.7,cex.title=.8), main=NA,
                     index=NULL, xlab, ylab, z.var, colors, cex, lty, 
                     lwd, ylim, xlim, factor.names, band.transparency, 
                     band.colors, grid=TRUE, axes, lattice, rotx, roty,
                     symbols=list(pch = 19), ticks.x, lines=TRUE, robust=FALSE,
                     cluster, plot=TRUE, dots=FALSE, ...) 
{
  
  if(!inherits(fit,c("rma.mv","rma","rma.uni"))) stop("Model is not 'rma()' or 'rma.mv()'.", call. = FALSE)
  
  lm_fit <- lm(fixed_form_rma(fit), data = get_data_(fit), na.action = "na.omit")
  
  is_singular <- anyNA(coef(lm_fit))
  
  fit2 <- rma2gls(fit)
  
  if(is_singular) { 
    
    fit2$coefficients <- replace(lm_fit$coefficients, !is.na(lm_fit$coefficients), fit2$coefficients)
  }
  
  x <- if(!full) allEffects(fit2, ...) else predictorEffects(fit2, ...)
  
  x <- Filter(function(i) !all(is.na(i$fit)), x)
  
  if(!is.null(index)){ 
    
    len <- length(x)  
    index <- index[index <= len & index >= 1] 
    if(length(index)==0) index <- NULL
  }            
  
  x <- if(!is.null(index)) x[index] else x
  
  if(plot){   
    
    if(missing(ylab)) ylab <- paste0("Effect Size (",as.character(fixed_form_rma(fit))[2],")")
    
    xcv <- plot(x, multiline=multiline, main=main, rug=FALSE, dots=dots,
                confint=confint, x.var=x.var, z.var=z.var, key.args=key.args, 
                xlab=xlab, ylab=ylab, colors=colors, cex=cex, lty=lty, 
                lwd=lwd, ylim=ylim, xlim=xlim, factor.names=factor.names, band.transparency=band.transparency, 
                band.colors=band.colors, grid=grid, lattice=lattice, axes=axes, rotx=rotx, roty=roty, symbols=symbols,
                ticks.x=ticks.x, lines=lines)
    
    xcv$x.scales$tck=c(1,0)
    xcv$y.scales$tck=c(1,0)
    
    return(invisible(as.data.frame(x)))
    
    xcv 
    
  } else {
    
    return(invisible(as.data.frame(x)))
  }
}


# H=================================================================================================================================================

random_left <- function(random_fml) {
  
  as.formula(as.character(parse(text = sub("\\|.*", "", random_fml))))
  
}

# H=================================================================================================================================================

rma2gls <- function(fit){
  
  data_. <- get_data_(fit)
  form_. <- fixed_form_rma(fit)
  
  rownames(fit$b)[rownames(fit$b) %in% "intrcpt"] <- "(Intercept)"
  rownames(fit$beta)[rownames(fit$beta) %in% "intrcpt"] <- "(Intercept)"
  rownames(fit$vb)[rownames(fit$vb) %in% "intrcpt"] <- "(Intercept)"
  colnames(fit$vb)[colnames(fit$vb) %in% "intrcpt"] <- "(Intercept)"
  
  fit2 <- nlme::gls(form_., data = data_., method = fit$method, na.action = "na.omit",
                    control = glsControl(singular.ok = TRUE))
  
  fit2$call$model <- form_.
  fit2$call$data <- data_.
  
  fit2$coefficients <- coef.rma(fit) 
  fit2$varBeta <- fit$vb
  
  return(fit2)
} 

# H================================================================================================================================================

clean_reg <- function(fm, nm, uniq = TRUE) {
  
  vars <- vapply(attr(terms(fm), "variables"), deparse, "")[-1L]
  subpat <- paste0(gsub("([()])", "\\\\\\1", vars), collapse = "|")
  l <- rapply(strsplit(nm, ":"), sub, how = "list",
              perl = TRUE,
              pattern = sprintf("^(?!(%1$s)$)(%1$s)(.+)$", subpat),
              replacement = "\\3")
  vec <- vapply(l, paste0, "", collapse = ":")
  if(uniq) vec <- make.unique(vec)
  vec[vec=="intrcpt"] <- "(Intercept)"
  return(vec)
}     

# H=================================================================================================================================================     

any_num_vec <- function(vec){
  any(grepl("^[0-9]{1,}$", gsub(":", "", vec)))
}

# M=================================================================================================================================================    

results_rma <- function(fit, digits = 3, robust = TRUE, blank_sign = "", 
                        cat_shown = 1, shift_up = NULL, shift_down = NULL, 
                        drop_rows = NULL, drop_cols = NULL, QM = TRUE, 
                        QE = FALSE, sig = TRUE, clean_names = NULL, 
                        tidy = FALSE, tol_large = 1e4, random_only = FALSE){
  
  if(!inherits(fit, "rma.mv")) stop("Model is not 'rma.mv()'.", call. = FALSE)
  
  fixed_eff <- is.null(fit$random)
  cr <- if(!fixed_eff) is_crossed(fit) else FALSE
  
  lm_fit <- lm(fixed_form_rma(fit), data = get_data_(fit), na.action = "na.omit")
  
  cl <- clean_reg(lm_fit, names(coef(lm_fit)))
  
  if(is.null(clean_names)){
    
    if(any_num_vec(cl)) {
      
      clean_names <- FALSE
      
    } else { 
      
      clean_names <- TRUE
      
    }
  } 
  
  if(clean_names) names(lm_fit$coefficients) <- cl
  if(clean_names) fit <- clean_reg_names(fit)
  
  lm_coef <- coef(lm_fit)
  
  is_singular <- anyNA(lm_coef)
  
  if(is_singular) message("Note:",dQuote(toString(names(lm_coef)[is.na(lm_coef)])), " dropped due to lack of data.\n")
  
  if(robust & any(cr) || robust & fixed_eff) { 
    
    robust <- FALSE
    message("Note: Robust estimation not available for models with", if(any(cr))" crossed random-" else " only fixed-", "effects.\n")
  }
  
  res <- if(!robust) { 
    
    a <- coef(summary(fit))
    
    nm <- c("Estimate","SE","t","Df","p-value","Lower","Upper")
    
    nm <- if(fit$test == "t") { nm } else { nm[3] <- "z";
    nm[-4] }
    
    setNames(a, nm)
    
  } else {
    
    a <- suppressWarnings(as.data.frame(conf_int(fit, vcov = "CR2"))[-1])
    b <- suppressWarnings(coef_test(fit, vcov = "CR2"))
    a$p <- b$p_Satt
    a$t <- b$tstat
    
    a <- a[c(1:2,7,3,6,4:5)]
    
    setNames(a, c("Estimate","SE","t","Df","p-value","Lower","Upper"))
  }
  
  rn <- rownames(res)[1]
  
  if(rn == "intrcpt") rownames(res)[1] <- "(Intercept)"
  
  res_org <- res
  res <- na.omit(res)
  
  if(robust & nrow(res) != nrow(res_org)) message("Note:",dQuote(toString(setdiff(rownames(res_org),rownames(res)))), " dropped due to inestimablity under Robust estimation.\n")
  
  if(random_only) 
  { QE <- FALSE  
  QM <- FALSE
  sig <- FALSE
  drop_cols <- 1:7
  drop_rows <- 1:(nrow(res)+1)
  }
  
  if(QE){
    qe <- data.frame(Estimate = fit$QE, Df = nobs.rma(fit), 
                     pval = fit$QEp, row.names = "QE") %>%
      dplyr::rename("p-value"="pval") 
    
    res <- bind_rows(res, qe)
  }
  
  if(QM){
    
    mc <- try(clubSandwich::Wald_test(fit, constrain_zero(fit$btt), "CR2"), silent = TRUE)   
    
    bad <- inherits(mc,"try-error")
    
    if(robust){
      
      if(bad || !bad && is.na(mc$p_val)) { 
        robust <- FALSE
        message("Note: Robust QM unavailable,likely: \n1- Some moderators in <2 clusters OR/AND \n2- High # of coefficients vs. # of highest clusters.\nQM results are model-based.\n")
      }
      
      if(!bad && mc$Fstat>tol_large) { message("Note: Robust QM is unreasonably large, likely robust estimation is unfit for the model (use 'robust=FALSE').") }
    }
    
    qm <- if(robust) {
      
      data.frame(Estimate = mc$Fstat, Df = mc$df_num, 
                 pval = mc$p_val, row.names = "QM") %>%
        dplyr::rename("p-value"="pval") 
      
    } else {
      
      data.frame(Estimate = fit$QM, Df = fit$QMdf[1], 
                 pval = fit$QMp, row.names = "QM") %>%
        dplyr::rename("p-value"="pval") 
      
    }
    res <- bind_rows(res, qm)
  }
  
  u <- get_error_rho(fit)
  cte <- length(u) == 1
  
  d6 <- data.frame(r = if(cte) u else mean(u, na.rm = TRUE), 
                   row.names = paste0("Within Corr. (",if(cte) "constant" else "average",")"))
  
  blk <- paste0(paste0(rep(" ",digits-1), collapse=""), "NA", collapse ="")
  
  if(fixed_eff) { 
    
    out <- roundi(dplyr::bind_rows(res, d6), digits = digits)
    out[out== blk] <- blank_sign
    
    if(sig){ 
      
      p.values <- as.numeric(out$"p-value")
      
      Signif <- symnum(p.values, corr = FALSE, 
                       na = FALSE, cutpoints = 
                         c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                       symbols = c("***", "**", "*", ".", " ")) 
      
      out <- add_column(out, Sig. = Signif, .after = "p-value")
    }
    
    if(!is.null(shift_up)) out <- shift_rows(out, shift_up)
    if(!is.null(shift_down)) out <- shift_rows(out, shift_down, up = FALSE)
    if(!is.null(drop_rows)) out <- out[-drop_rows, ]
    
    out <- dplyr::select(out, -tidyselect::all_of(drop_cols))
    
    if(tidy) out <- cbind(Terms = rownames(out), set_rownames_(out, NULL))
    
    return(out)
  }
  
  res <- rbind(res, "|RANDOM|" = NA)
  
  Sys.setlocale(locale = "Greek")
  
  if(fit$withS){
    
    d1 <- data.frame(Sigma = sqrt(fit$sigma2), 
                     row.names = paste0(names(cr), ifelse(cr," (crossed)"," (nested)"))) 
    
    d1 <- setNames(d1, intToUtf8(963))
  } else { d1 <- NULL}
  
  if(fit$withG){
    
    h <- paste(fit$struct[1], "Corr.")
    is_un <- fit$struct[1] == "UN"
    is_gen <- fit$struct[1] == "GEN"
    is_diag <- fit$struct[1] == "DIAG"
    is_simple <- length(fit$tau2) == 1
    
    rnm <- paste("Outer:", tail(fit$g.names,1))
    clnm <- clean_GH_names(fit)
    
    d2 <- data.frame(Tau = sqrt(fit$tau2), 
                     row.names = paste0(if(!is_simple) clnm else fit$g.names[1],
                                        paste0(if(is_diag)" (Uncor. " 
                                               else " (Cor. ",if(!is_simple & !is_gen) paste0(" ", fit$g.names[1]),")")))
    
    d2 <- setNames(d2, intToUtf8(964))
    
    d2 <- rbind(NA, d2)
    rownames(d2)[1] <- rnm
    
    d3 <- data.frame(Rho = fit$rho, 
                     row.names = if(is_un || is_gen) apply(combn(clnm,2),2,paste0, collapse = "~") 
                     else paste0(h," (",shorten_(clnm, cat_shown),")")) 
    
    d3 <- setNames(d3, intToUtf8(961))
    
  } else { d2 <- NULL; d3 <- NULL}
  
  if(fit$withH){
    
    h <- paste(fit$struct[2], "Corr.")
    is_un <- fit$struct[2] == "UN"
    is_gen <- fit$struct[2] == "GEN"
    is_diag <- fit$struct[2] == "DIAG"
    is_simple <- length(fit$gamma2) == 1
    
    rnm <- paste("Outer:", paste0(tail(fit$h.names,1)," "))
    
    clnm <- clean_GH_names(fit, G=FALSE)
    
    d4 <- data.frame(Gamma = sqrt(fit$gamma2), 
                     row.names = paste0(if(!is_simple) clnm else fit$h.names[1],
                                        paste0(if(is_diag)" (Uncor. " 
                                               else " (Cor. ",if(!is_simple) paste0(" ",if(!is_gen)fit$h.names[1]),") "))) 
    
    d4 <- setNames(d4, intToUtf8(933))
    
    d4 <- rbind(NA, d4)
    rownames(d4)[1] <- rnm
    
    d5 <- data.frame(Phi = fit$phi, 
                     row.names = if(is_un || is_gen) apply(combn(clnm,2),2,paste0, collapse = "~ ")
                     else paste0(h," (",shorten_(clnm, cat_shown),") "))
    
    d5 <- setNames(d5, intToUtf8(966))
    
  } else { d4 <- NULL; d5 <- NULL}
  
  out <- roundi(dplyr::bind_rows(res, d1, d2, d3, d4, d5, d6), digits = digits)
  
  
  out[out== blk] <- blank_sign
  
  if(sig){ 
    
    p.values <- as.numeric(out$"p-value")
    
    Signif <- symnum(p.values, corr = FALSE, 
                     na = FALSE, cutpoints = 
                       c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                     symbols = c("***", "**", "*", ".", " ")) 
    
    out <- add_column(out, Sig. = Signif, .after = "p-value")
  }
  
  if(!is.null(shift_up)) out <- shift_rows(out, shift_up)
  if(!is.null(shift_down)) out <- shift_rows(out, shift_down, up = FALSE)
  if(!is.null(drop_rows)) out <- out[-drop_rows, ]
  if(!is.null(drop_cols)) out <- dplyr::select(out, -tidyselect::all_of(drop_cols))
  
  if(tidy) out <- cbind(Terms = rownames(out), set_rownames_(out, NULL))
  
  return(out)
}                      

# H=================================================================================================================================================                     

roundi <- function(x, digits = 7){
  
  if(!inherits(x, c("data.frame","tibble","matrix"))) stop("'x' must be a 'data.frame' or a 'matrix'.", call. = FALSE)
  if(inherits(x,"matrix")) x <- as.data.frame(x)
  
  num <- sapply(x, is.numeric)
  
  x[num] <- lapply(x[num], function(i) formatC(round(i, digits), digits, format = "f"))
  
  return(x)
}                      

# H=================================================================================================================================================

random_GH_form <- function(fit, G = TRUE){
  
  fm <- fit$random
  if(G) fm[[1]] else fm[[2]]
  
}

# H=================================================================================================================================================

clean_GH_names <- function(fit, G = TRUE) {
  
  fmla <- random_left(random_GH_form(fit, G = G))
  vec <- if(G) rownames(fit$G) else rownames(fit$H)
  clean_reg(fmla, vec)
  
}                                     

# H=================================================================================================================================================

clean_reg_names <- function(fit) {
  
  fmla <- fixed_form_rma(fit)
  vec <- rownames(fit$b)
  
  vec <- clean_reg(fmla, vec)
  if(fit$int.only) vec[vec=="(Intercept)"||vec==""] <- "Overall Effect"
  rownames(fit$b) <- vec
  rownames(fit$beta) <- vec
  rownames(fit$vb) <- colnames(fit$vb) <- vec
  return(fit)
}

# H=================================================================================================================================================                       

set_rownames_ <- function (object = nm, nm) 
{
  rownames(object) <- nm
  object
}       

# M=================================================================================================================================================

smooth_vi <- function(data, study, vi, digits = 8, fun = sd, ylab = "Studies", xlab = NULL, max_bar = FALSE, return_list = TRUE,
                      breaks = "Sturges", main = "Closeness of Sampling Variances\n (in each study)"){
  
  dot_cols <- rlang::ensyms(study, vi, fun)
  str_cols <- purrr::map_chr(dot_cols, rlang::as_string)
  
  d_clean <- full_clean(data)
  
  out <- map_dfr(group_split(d_clean, study), 
                 ~data.frame(study = unique(.[[str_cols[[1]]]]), sd = round(fun(.[[str_cols[[2]]]]),digits)))
  
  names(out)[2] <- str_cols[[3]]
  
  res <- hist(out[[2]], plot = FALSE, breaks = breaks)
  graphics.off()
  
  hist(out[[2]], breaks = breaks, xaxt = "n", yaxt = "n", cex.axis = .8, mgp = c(1.5,.4,0), xlab = if(is.null(xlab)) toupper(str_cols[[3]]) else xlab,
       main = main, cex.lab = .9, ylab = ylab, cex.main = .8)
  
  axis(1, at = res$breaks, las=1, cex.axis=.8, mgp=c(1.5,.4,0))
  axis(2, at = c(axTicks(2), max(res$counts)), las=1, cex.axis=.8, mgp=c(1.5,.55,0))
  text(res$mids[res$counts!=0], res$counts[res$counts!=0], res$counts[res$counts!=0], pos = 3, col = 4, font = 2, xpd = NA, cex = .75)
  
  whichbar <- findInterval(out[[2]], res$breaks)
  mm <- sort(unique(whichbar[!is.na(whichbar)]))
  
  out <- if(!max_bar) setNames(lapply(mm, function(i) na.omit(out[whichbar==i,])), 
                               paste0("Bar",seq_along(mm))) else na.omit(out[whichbar==which.max(res$counts),])
  
  return(invisible(if(return_list) out else dplyr::bind_rows(out, .id = "Bar")))
}


# H=================================================================================================================================================

# Testing jointly whether EMMs *across* multiple
# categorical variables in a MUTOS block are equal to each other or not?

comb.facs <- function(object, facs, new.name = paste(facs, collapse = ".")) {
  
  if((length(facs) < 2))
    stop("Block test needs at least two categorical variables.", call. = FALSE)
  levs = object@levels
  if(any(sapply(facs, function(x) !(x %in% names(levs)))))
    stop("Categorical variables not found.", call. = FALSE)
  grid = object@grid
  idx = sapply(facs, function(x) which(names(levs) == x))
  levs[[new.name]] = do.call(paste, c(do.call(expand.grid, levs[idx]), sep = ":"))
  grid[[new.name]] = do.call(paste, c(grid[idx], sep = ":"))
  levs[idx] = grid[idx] = NULL
  object@levels = levs
  object@grid = grid
  object@roles$predictors = c(object@roles$predictors[-idx], new.name)
  ord = emmeans:::.std.order(grid, levs)
  object[ord]
}


# M=================================================================================================================================================

post_rma <- function(fit, specs = NULL, cont_var = NULL, by = NULL, horiz = TRUE, ci = TRUE, block_vars = NULL,
                     adjust = "none", mutos = FALSE, mutos_contrast = FALSE, compare = FALSE, plot_pairwise = FALSE, p_value = TRUE,
                     reverse = FALSE, digits = 3, xlab = "Estimated Effect", shift_up = NULL, shift_down = NULL, 
                     drop_rows = NULL, mutos_name = "(M)UTOS Term", drop_cols = NULL, contrast_contrast = FALSE, 
                     na.rm = TRUE, robust = FALSE, cluster, show0df = FALSE, sig = TRUE, contr,
                     get_rows = NULL, get_cols = NULL, ...){
  
  if(!inherits(fit, "rma.mv")) stop("Model is not 'rma.mv()'.", call. = FALSE)
  
  dot_args <- list(...)
  dot_args_nm <- names(dot_args)
  
  cl <- match.call()
  
  if(robust) { 
    
    fixed_eff <- is.null(fit$random)
    cr <- if(!fixed_eff) is_crossed(fit) else FALSE
    
    if(any(cr) || fixed_eff) { 
      
      robust <- FALSE
      message("Robust estimation not available for models with", if(any(cr))" crossed random-" else " only fixed-", "effects.")
    }
    
    vcov_. <- try(as.matrix(vcovCR(fit, type = "CR2", cluster = cluster)), silent=TRUE)
    
    if(inherits(vcov_., "try-error")) { 
      robust <- FALSE
      message("Robust tests unavailable (likely having <2 clusters for some moderators).\nResults are model-based.")
    }
  }
  
  data_. <- get_data_(fit)
  fml <- fixed_form_rma(fit)
  lm_fit <- lm(fml, data = data_.)
  lm_fit$call$data <- data_.
  lm_fit$call$formula <- fml
  
  is_singular <- anyNA(coef(lm_fit))
  
  rma.mv_fit <- fit
  
  fit <- rma2gls(fit)
  
  specs_org <- specs
  if(is.null(specs)) specs <- as.formula(bquote(~.(terms(fit)[[3]])))
  
  if(robust) { 
    
    rownames(vcov_.)[rownames(vcov_.) %in% "intrcpt"] <- "(Intercept)"
    colnames(vcov_.)[colnames(vcov_.) %in% "intrcpt"] <- "(Intercept)"
    
    fit$varBeta <- vcov_.
  }
  
  
  fit_org <- fit
  
  if(is_singular) { 
    
    lm_coef <- lm_fit$coefficients
    gls_coef <- fit$coefficients
    
    fit$coefficients <- replace(lm_coef, !is.na(lm_coef), gls_coef)
    
    fit <- suppressMessages(emmeans::ref_grid(fit))
    fit@nbasis <- suppressMessages(emmeans::ref_grid(lm_fit)@nbasis)
  }
  
  infer <- c(ci, p_value)
  
  lookup <- c(Contrast="contrast",Estimate="estimate","Mean"="emmean","Response"="response",t="t.ratio",
              Df="df","p-value"="p.value",Lower="lower.CL",Upper="upper.CL",
              Df1="df1", Df2="df2","F"="F.ratio",m="model term")
  
  names(lookup)[12] <- if(is.null(block_vars)) paste(mutos_name, "Contrast") else "Block Contrast"
  
  tran. <- if('tran' %in% dot_args_nm) dot_args$tran else FALSE           
  type. <- if('type' %in% dot_args_nm) dot_args$type else FALSE
  
  is_contr <- !missing(contr)            
  
  ems <- try(if(is.null(cont_var)){
    
    emmeans(object = fit, specs = specs, infer = infer, adjust = adjust, contr = contr, ...)
    
  } else {
    
    if(!is_contr){ 
      
      emtrends(object = fit, specs = specs, var = cont_var, infer = infer, adjust = adjust, ...)
      
    } else {
      
      emtrends(object = fit, specs = specs, var = cont_var, infer = infer, adjust = adjust, contr = contr, ...)
      
    }
    
  }, silent = TRUE)
  
  
  if(inherits(ems,"try-error")) return(message("Error: Wrong specification OR no relavant data for the comparisons found."))
  
  con_methods <- c("pairwise","revpairwise","tukey","consec",
                   "poly","trt.vs.ctrl","trt.vs.ctrlk","trt.vs.ctrl1",
                   "dunnett","mean_chg","eff","del.eff","identity")
  
  is_pair <- any(con_methods %in% as.character(specs))
  
  if(!is.null(cont_var) & is_pair) names(lookup)[2] <- paste0(cont_var,".dif")
  
  out <- if(is_pair){
    
    methd <- as.character(specs[2])
    
    if(plot_pairwise) print(plot(ems, by = by, comparisons = compare, horizontal = horiz, adjust = adjust, xlab = xlab)) 
    
    pp <- contrast(ems, method = methd, each="simple", infer=infer, reverse=reverse, adjust=adjust)[[if(!contrast_contrast) 1 else 2]]
    
    if(!is_contr) pp else contrast(pp, contr)
    
  }
  
  else {
    
    if(mutos_contrast) {  
      
      fit_used <- if(is.null(specs_org)) fit else ems
      
      joint_tests(fit_used, by = by, adjust = adjust, show0df = show0df, ...)
      
    } else if (mutos){
      
      if(is.null(specs_org)){
        
        zz <- cbind(mod=get_vars_(fit_org,as_fml=FALSE),as.data.frame(map_dfr(get_vars_(fit_org),~emmeans::test(emmeans(ems,.),joint=TRUE))))
        names(zz)[1] <- mutos_name
        zz 
        
      } else { 
        
        cbind("(M)UTOS Test"= as.character(specs)[2],as.data.frame(emmeans::test(ems, joint=TRUE)))
        
      }
      
    } else if(!is.null(block_vars)){ 
      
      com <- comb.facs(ems, block_vars)
      
      joint_tests(com)
      
    } else {
      
      ems
    }
  }
  
  out <- as.data.frame(out, adjust = adjust) %>%
    dplyr::rename(tidyselect::any_of(lookup)) %>% 
    dplyr::select(-tidyselect::any_of("note"))
  
  if(!is.null(block_vars)) out <- subset(out,`Block Contrast`==paste0(block_vars, collapse="."))
  
  out <- set_rownames_(out,NULL)
  
  p.values <- as.numeric(out$"p-value")
  
  if(all(is.na(p.values))) { 
    return(message("Error: No relavant data for the comparisons found."))
  }
  
  if(sig){
    Signif <- symnum(p.values, corr = FALSE, 
                     na = FALSE, cutpoints = 
                       c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                     symbols = c("***", "**", "*", ".", " "))
    
    out <- tibble::add_column(out, Sig. = Signif, .after = "p-value")
  }
  
  if(na.rm) out <- na.omit(out)
  
  out <- roundi(out, digits = digits)
  
  if(!is.null(shift_up)) out <- shift_rows(out, shift_up)
  if(!is.null(shift_down)) out <- shift_rows(out, shift_down, up = FALSE)
  if(!is.null(drop_rows)) out <- out[-drop_rows, ]
  if(!is.null(get_rows)) out <- out[get_rows, ]
  
  if(!is.null(drop_cols)) out <- dplyr::select(out, -tidyselect::all_of(drop_cols))
  if(!is.null(get_cols)) out <- dplyr::select(out, tidyselect::all_of(get_cols))
  
  out <- list(table = out, specs = specs, call = cl, fit = fit, rma.mv_fit = rma.mv_fit, ems = ems,
              tran. = tran., type. = type.)
  
  class(out) <- "post_rma"
  return(out)
}                  

# M=================================================================================================================================================

R2_rma <- function(..., robust = TRUE, digits = 3, 
                   model_names = NULL, null_model = NULL,
                   level_names = NULL, blank_sign = "", 
                   null_name = "No (M)UTOS", tol_large = 1e4)
{
  
  LL <- list(...)
  if(!all(sapply(LL,inherits,"rma.mv"))) stop("All models must be 'rma.mv()'.", call. = FALSE)
  
  bad <- sapply(LL, function(i) i$withG || i$withH || i$withR)
  
  if(any(bad)) stop("These models not yet supported.", call. = FALSE)
  
  ok <- length(unique(map(LL,~map_chr(strsplit(.$s.names,"/",fixed=TRUE),tail,1))))==1
  
  if(!ok) stop("Models must have the exact same random-effects.", call.=FALSE)
  
  first <- LL[[1]]
  
  Model <- if(is.null(model_names)) as.character(substitute(...())) else model_names
  
  .zolqui_. <- as.formula(paste0(as.character(fixed_form_rma(first))[2],"~1"))
  
  null_fit <- if(is.null(null_model)) update.rma(first, yi = .zolqui_.) else null_model
  
  lvl_names <- if(is.null(level_names)) sapply(strsplit(null_fit$s.names,"/",fixed=TRUE),tail,1) else level_names
  
  sigmasn <- setNames(sqrt(null_fit$sigma2), lvl_names)
  
  sigma_totaln <- sqrt(sum(sigmasn^2))
  
  null_res <- data.frame(Model = null_name,._A_.= sigma_totaln,._D_.=NA,R2=NA)
  
  null_res <- add_column(null_res, as.data.frame(t(sigmasn)), .after = "._A_.")
  
  z <- function(nm) paste(intToUtf8(963),nm) 
  
  Sys.setlocale(locale = "Greek")
  
  f <- function(fit){
    
    if(robust){
      
      mc <- try(clubSandwich::Wald_test(fit, constrain_zero(fit$btt), "CR2"), silent=TRUE)   
      
      bad <- inherits(mc,"try-error")
      
      if(bad || !bad && is.na(mc$p_val)) { 
        robust <- FALSE
        message("Note: Robust p-value unavailable,likely: \n1- Some moderators in <2 clusters OR/AND \n2- High # of coefficients vs. # of highest clusters.\np-value is model-based.\n")
      }
      
      if(!bad && mc$Fstat>tol_large) { message("Note: Robust estimation seems unfit for the model (use 'robust=FALSE').") } 
    }
    
    p <- if(robust) mc$p_val else fit$QMp
    
    sigmas <- setNames(sqrt(fit$sigma2), lvl_names) 
    
    sigma_total <- sqrt(sum(sigmas^2)) 
    
    R2 <- (sigma_totaln - sigma_total) / sigma_totaln*1e2
    
    if(R2<0) message("Negative R2 was set to 0.")
    
    R2 <- if(R2>0) R2 else 0
    
    c(._A_.=sigma_total,sigmas,._D_.=p, R2=R2)
  }
  
  out <- map_dfr(LL, f) %>% as.data.frame() %>%
    add_column(Model = Model, .before = "._A_.")
  
  res <- roundi(bind_rows(null_res, out), digits = digits)
  
  res[-1,]$R2 <- paste0(res[-1,]$R2,"%")
  
  names(res)[names(res) %in% c("._A_.",lvl_names,"._D_.")] <- c(z(c("total",lvl_names)),"p-value")
  
  blk <- paste0(paste0(rep(" ",digits-1), collapse=""), "NA", collapse ="")
  
  res[res == blk] <- blank_sign
  
  return(res)
}                                                                        


# H=================================================================================================================================================

sizetree2 <- function (x, left = 0, top, right = 1, lastcenter = NA,
                       showval = TRUE, showcount = FALSE, stacklabels = TRUE, firstcall = TRUE,
                       col = NULL, border = NA, toplab = NULL, base.cex = 1, cex_top = 1, ...) 
{
  dimx <- dim(x)
  colname <- names(x)[1]
  if (firstcall) {
    x <- x[do.call(order, x), ]
    oldmar <- par("mar")
    par(mar = c(1, 2, 2, 1))
    top <- sum(!is.na(x[, 1]))
    if (top < dimx[1]) 
      cat(dimx[1] - top, "NA values dropped from first stack.\n")
    plot(0, xlim = c(0, dimx[2]), ylim = c(0, top), type = "n", 
         axes = FALSE, xlab = "", ylab = "", ...)
  }
  xfreq <- table(x[, 1])
  lenxf <- length(xfreq)
  if (firstcall) {
    if (is.null(col)) {
      col <- list()
      for (index in 1:dimx[2]) col[[index]] <- rainbow(length(table(x[, 
                                                                      index])))
    }
    for (index in 1:dimx[2]) if (is.null(names(col[[index]]))) 
      names(col[[index]]) <- names(table(x[, index]))
  }
  if (lenxf) {
    if (is.list(col)) {
      barcol <- col[[1]]
      barcol <- barcol[names(col[[1]]) %in% names(xfreq)]
    }
    else barcol <- col[names(col) %in% names(xfreq)]
    labels <- names(xfreq)
    squeeze <- (right - left)/10
    for (bar in 1:lenxf) {
      if (length(xfreq[bar])) {
        if (!is.na(xfreq[bar])) {
          if (xfreq[bar] > 0) {
            rect(left + squeeze, top - xfreq[bar], right - 
                   squeeze, top, col = barcol[bar], border = border)
            labelheight <- strheight(labels[bar])
            cex <- ifelse((1.5 * labelheight) > xfreq[bar], 
                          base.cex * 0.75 * xfreq[bar]/labelheight, 
                          base.cex)
            if (showval) {
              textcol <- ifelse(colSums(col2rgb(unlist(barcol[bar])) * 
                                          c(1.4, 1.4, 0.5)) < 350, "white", "black")
              bartext <- ifelse(showcount, paste(labels[bar], 
                                                 " (", xfreq[bar], ")", sep = ""), labels[bar])
              text((left + right)/2, top - xfreq[bar]/2, 
                   bartext, cex = cex, col = textcol)
            }
            if (!is.na(lastcenter)) 
              segments(left + squeeze, top - xfreq[bar]/2, 
                       left - squeeze, lastcenter)
            xvalue <- ifelse(is.numeric(x[, 1]), as.numeric(labels[bar]), 
                             labels[bar])
            if (dimx[2] > 1) {
              newcol <- col
              newcol[[1]] <- NULL
              nextx <- subset(x, x[, 1] == xvalue, 2:dimx[2])
              sizetree2(nextx, right, top, right + 1, 
                        lastcenter = top - xfreq[bar]/2, showval = showval,
                        showcount = showcount, stacklabels = stacklabels, firstcall = FALSE,
                        col = newcol, border = border, base.cex = base.cex)
            }
          }
        }
      }
      top <- top - xfreq[bar]
    }
    if (stacklabels) 
      mtext(colname, side = 1, at = (left + right)/2, line = -0.4)
  }
  if (firstcall) {
    if (!is.null(toplab)) {
      par(xpd = TRUE)
      top <- sum(!is.na(x[, 1]))
      text(seq(0.5,dimx[2] - 0.5,by=1), 1.01 * top, toplab, adj = c(0.5,
                                                                    0), cex = cex_top)
      par(xpd = FALSE)
    }
    par(mar = oldmar)
  }
}

#================================================================================================================================================

contr_rma <- function(post_rma_fit, contr_index){
  
  if(!inherits(post_rma_fit, "post_rma")) stop("post_rma_fit is not 'post_rma()'.", call. = FALSE)
  
  post_rma_fit <- post_rma_fit$table
  
  all_emms <- rep(0, nrow(post_rma_fit))
  
  all_emms[abs(contr_index)] <- sign(contr_index)
  
  return(all_emms)
}                

#===================================================================================================================================================

prob_rma <- function(post_rma_fit, target_effect = 0, condition = c("or larger", "or smaller"), gain = FALSE, none_names = NULL, ...){
  
  
  if(!inherits(post_rma_fit, "post_rma")) stop("post_rma_fit is not 'post_rma()'.", call. = FALSE)   
  
  fit <- post_rma_fit$rma.mv_fit
  
  if(fit$withG || fit$withH || fit$withR) stop("These models not yet supported.", call. = FALSE)
  
  specs <- post_rma_fit$specs
  
  post_rma_fit <- type.convert(post_rma_fit$table, as.is=TRUE)
  
  nms <- names(post_rma_fit)
  
  contr <- if("contr" %in% as.character(post_rma_fit$call) || c("pairwise","revpairwise","tukey","consec",
                                                                "poly","trt.vs.ctrl","trt.vs.ctrlk","trt.vs.ctrl1",
                                                                "dunnett","mean_chg","eff","del.eff","identity") %in% as.character(specs)) "Contrast" else NULL
  
  vv <- nms[!nms %in% c("Mean","Response","SE","Df","Lower","Upper","t",      
                        "p-value","Sig.",contr,"F","Df1","Df2",
                        "Estimate","m","Block Contrast","(M)UTOS Term", none_names)]
  
  Term <-sapply(seq_len(nrow(post_rma_fit)), 
                function(i) paste0(as.vector(unlist(post_rma_fit[vv][i,])), collapse = " "))
  
  
  cond <- match.arg(condition)  
  
  lower.tail <- switch(cond, 
                       "or larger" = FALSE, 
                       "or smaller" = TRUE)
  
  lower_ave_eff <- post_rma_fit$Lower
  
  upper_ave_eff <- post_rma_fit$Upper
  
  ave_eff <- if(!is.null(post_rma_fit$Mean)) post_rma_fit$Mean else if(!is.null(post_rma_fit$Response)) post_rma_fit$Response else post_rma_fit$Estimate
  
  ci <- as.data.frame(confint.rma.mv(fit, ...))
  
  all_lvls <- ci[odds_.(seq_len(nrow(ci))), , drop=FALSE]
  
  #total_sd: estimate, lower, upper
  total_sd <- if(!gain) sqrt(colSums(all_lvls)) else sqrt(2)*ci[nrow(ci),] #sqrt(2) * sqrt(colSums(all_lvls[-1, ,drop=FALSE]))
  
  # Probability at the estimates
  Probability <- paste0(formatC(round(pnorm(target_effect, ave_eff, total_sd[1], lower.tail=lower.tail), 4)*1e2,digits = 2, format = "f"),"%")
  
  # Probability at the lowest ave_eff and highest variability
  min_Probability <- paste0(formatC(round(pnorm(target_effect, lower_ave_eff, total_sd[3], lower.tail=lower.tail), 4)*1e2,digits = 2, format = "f"),"%")
  
  # Probability at the highest ave_eff and lowest variability
  max_Probability <- paste0(formatC(round(pnorm(target_effect, upper_ave_eff, total_sd[2], lower.tail=lower.tail), 4)*1e2, digits = 2, format = "f"),"%")
  
  data.frame(Term=Term, Target_Effect = paste(target_effect, cond, collapse = " "), Probability = Probability, 
             Min = min_Probability, Max = max_Probability)
  
}

#M==============================================================================================================================================

sense_rma <- function(post_rma_fit = NULL, var_name, fit = NULL, 
                      r = (3:7)*.1, cluster = NULL, clean_names = NULL,
                      regression = NULL, label_lines = TRUE, none_names=NULL,
                      cex_labels = .55, plot = TRUE, digits = 3, ...){
  
  
  if(is.null(fit) & is.null(post_rma_fit)) stop("Provide either 'fit=' or 'post_rma_fit='.", call. = TRUE)
  if(!is.null(fit) & !inherits(fit, "rma.mv")) stop("Model is not 'rma.mv()'.", call. = FALSE)
  if(is.null(fit)) fit <- post_rma_fit$rma.mv_fit
  if(!is.null(post_rma_fit) & !inherits(post_rma_fit, "post_rma")) stop("post_rma_fit is not 'post_rma()'.", call. = FALSE) 
  if(fit$withG || fit$withH || fit$withR) stop("These models not yet supported.", call. = FALSE)
  
  tran. <- post_rma_fit$tran.
  type. <- post_rma_fit$type.
  
  dat <- get_data_(fit)
  
  regression <- if(is.null(regression) & !is.null(post_rma_fit)) {
    
    if("cont_var" %in% as.character(post_rma_fit$call)) TRUE else FALSE
    
  } else if(is.null(regression) & is.null(post_rma_fit)) TRUE else regression 
  
  
  if(regression){
    
    lm_fit <- lm(fixed_form_rma(fit), data = dat, na.action = "na.omit")
    
    cl <- clean_reg(lm_fit, names(coef(lm_fit)))
    
    if(is.null(clean_names)){
      
      if(any_num_vec(cl)) {
        
        clean_names <- FALSE
        
      } else { 
        
        clean_names <- TRUE
        
      }
    } 
    
    if(clean_names) names(lm_fit$coefficients) <- cl
    if(clean_names) fit <- clean_reg_names(fit)
  }   
  
  if(!is.null(post_rma_fit)){
    
    specs <- post_rma_fit$specs
    
    post_rma_fit <- post_rma_fit$table
  }
  
  cluster_name <- if(is.null(cluster)) strsplit(fit$s.names,"/",fixed=TRUE)[[1]] else cluster
  
  V_list <- lapply(r, function(i) impute_covariance_matrix(vi = dat[[var_name]], cluster=dat[[cluster_name]], r=i))
  
  model_list <- lapply(V_list, function(i) suppressWarnings(update.rma(fit, V = i)))
  
  xaxis_lab <- paste0("r=",r)
  
  total_hetros <- sapply(model_list, function(i) sqrt(sum(i$sigma2)))
  
  if(plot){
    
    graphics.off()
    org.par <- par(no.readonly = TRUE)
    on.exit(par(org.par))
    
    par(mfrow = c(2,1), mgp = c(1.5, 0.14, 0), mar = c(1.5, 2.6, 1.5, .5), 
        tck = -0.02)
  }  
  
  output <- if(regression){
    
    fixed_eff_list <- lapply(model_list, function(i) setNames(coef(summary(i))$estimate, rownames(fit$b)))
    if(plot){    
      
      matplot(t(as.data.frame(fixed_eff_list)), type = "l", xaxt = "n", ylab = "Estimates", ...)
      
      axis(1, at = axTicks(1), labels = xaxis_lab,...)
      
      mn <- mean(seq_len(length(fixed_eff_list)))
      
      if(label_lines) text(mn, as.data.frame(fixed_eff_list)[,mn], rownames(fit$b),
                           cex = cex_labels)
    }
    
    output <- as.data.frame(t(do.call(rbind, fixed_eff_list)))
    
    setNames(output, xaxis_lab)
    
  } else {
    
    if(!is.null(post_rma_fit)){
      
      nms <- names(post_rma_fit)
      
      vv <- nms[!nms %in% c("Mean","Response","SE","Df","Lower","Upper","t",      
                            "p-value","Sig.","Contrast","F","Df1","Df2",
                            "Estimate","m","Block Contrast","(M)UTOS Term", none_names=none_names)]
      
      Term <-sapply(seq_len(nrow(post_rma_fit)), 
                    function(i) paste0(as.vector(unlist(post_rma_fit[vv][i,])), collapse = " "))
      
      ave_col <- if(!is.null(post_rma_fit$Mean)) "Mean" else 
        if(!is.null(post_rma_fit$Estimate)) "Estimate" else "Response"
      
      post_rma_list <- lapply(model_list, function(i) 
        setNames(as.numeric(post_rma(i, specs, tran = tran., type = type.)$table[[ave_col]]),Term))
      
    } else {
      
      stop("Please provide a 'post_rma_fit' or use 'regression=TRUE'.", call. = FALSE)
      
    }
    
    if(plot){    
      
      matplot(t(as.data.frame(post_rma_list)), type = 'l', xaxt = "n", ylab = "Mean Effect", xlab = NA,...)
      
      axis(1, at = axTicks(1), labels = xaxis_lab,...)
      
      mn <- mean(seq_len(length(post_rma_list)))
      
      if(label_lines) text(mn, as.data.frame(post_rma_list)[,mn], Term,
                           cex = cex_labels)
    }    
    output <- as.data.frame(t(do.call(rbind, post_rma_list)))
    
    setNames(output, xaxis_lab)
    
  }
  
  if(plot){
    
    rng <- range(total_hetros)
    mrng <- mean(rng)
    plot(total_hetros, type = "l", ylim = rng+c(-mrng,mrng), xaxt = "n", xlab = NA, ylab = "Total Variation (SD)",...)
    axis(1, at = axTicks(1), labels = xaxis_lab, ...)
  }  
  out <- rbind(output, Total_variation_in_SD = total_hetros)
  out <- cbind(out, sd = sapply(1:nrow(out), function(i) sd(out[i,])))
  roundi(rownames_to_column(out, "Term"), digits = digits)
}                

#M================================================================================================================================================

plot_post_rma <- function(post_rma_fit, formula, ylab, CIs = FALSE, CIarg = list(lwd = .5, alpha = 1), ...){
  
  if(!inherits(post_rma_fit, "post_rma")) stop("post_rma_fit is not 'post_rma()'.", call. = FALSE)
  
  if(missing(ylab)) ylab <- paste0("Effect Size (",as.character(fixed_form_rma(post_rma_fit$rma.mv_fit))[2],")")
  
  emmip(object=post_rma_fit$ems, formula=formula, ylab=ylab, CIs=CIs, CIarg=CIarg, ...)
  
}                               

#M================================================================================================================================================

contrast_rma <- function(post_rma_fit, method, type,
                         digits = 2, ci=TRUE, 
                         p_value=TRUE, adjust="none",
                         na.rm = TRUE, sig=TRUE, ...){
  
  if(!inherits(post_rma_fit, "post_rma")) stop("post_rma_fit is not 'post_rma()'.", call. = FALSE)
  
  
  infer <- c(ci, p_value)
  
  lookup <- c(Contrast="contrast",Estimate="estimate","Mean"="emmean","Response"="response",t="t.ratio",
              Df="df","p-value"="p.value",Lower="lower.CL",Upper="upper.CL",
              Df1="df1", Df2="df2","F"="F.ratio",m="model term")
  
  con <- contrast(post_rma_fit$ems, method = method, type = type, infer = infer, ...)
  
  out <- as.data.frame(con, adjust = adjust) %>% 
    dplyr::rename(tidyselect::any_of(lookup)) %>% 
    dplyr::select(-tidyselect::any_of("note"))
  
  p.values <- as.numeric(out$"p-value")
  
  if(all(is.na(p.values))) { 
    return(message("Error: No relavant data for the comparisons found."))
  }
  
  if(sig){
    Signif <- symnum(p.values, corr = FALSE, 
                     na = FALSE, cutpoints = 
                       c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                     symbols = c("***", "**", "*", ".", " "))
    
    out <- tibble::add_column(out, Sig. = Signif, .after = "p-value")
  }
  
  if(na.rm) out <- na.omit(out)
  
  roundi(out, digits = digits)
}                                


#M================================================================================================================================================

r2z_tran <- list(
  linkfun = function(mu) atanh(mu),
  linkinv = function(eta) tanh(eta),
  mu.eta = function(eta) 1/cosh(eta)^2,
  valideta = function (eta) 
    all(is.finite(eta)) && all(eta > -1) && all(eta < 1),
  name = "r2z"
)                                

#======================== WCF Meta Dataset ======================================================================================================                

wcf <- read.csv("https://raw.githubusercontent.com/hkil/m/master/wcf.csv")

#=================================================================================================================================================

needzzsf <- c('metafor', 'clubSandwich', 'nlme', 'effects', 'lexicon', 'plotrix', 'rlang', 'emmeans','tidyverse')      

not.have23 <- needzzsf[!(needzzsf %in% installed.packages()[,"Package"])]
if(length(not.have23)) install.packages(not.have23)

suppressWarnings(
  suppressMessages({ 
    
    invisible(lapply(needzzsf, base::require, character.only = TRUE))
    
  }))

options(dplyr.summarise.inform = FALSE)                        
