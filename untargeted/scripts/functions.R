# Functions called by Control.Rmd or scripts within
# Should be sourced every time!

cv <- function(x, na.rm=FALSE, inf.rm=FALSE){
  if(inf.rm){
    x[is.infinite(x)] <- NA
  }
  sd(x = x, na.rm = na.rm)/mean(x = x, na.rm = na.rm)
}
speedyQscoreCalculator <- function(file_peaks, grabSingleFileData, qscoreCalculator){
  library(data.table)
  file_data <- grabSingleFileData(unique(file_peaks$filename))
  file_data_table <- as.data.table(file_data)
  peak_splits <- split(file_peaks, ceiling(seq_len(nrow(file_peaks))/100))
  file_qscores <- lapply(peak_splits, function(i){
    eic_many <- file_data_table[mz%between%c(min(i$mzmin), max(i$mzmax))]
    individual_peaks <- split(i, seq_len(nrow(i)))
    output <- sapply(individual_peaks, function(peak_row){
      eic <- eic_many[rt%between%c(peak_row$rtmin, peak_row$rtmax) & 
                        mz%between%c(peak_row$mzmin, peak_row$mzmax)]
      qscoreCalculator(eic)
    }, USE.NAMES = FALSE)
  })
  file_peaks$sn <- unlist(file_qscores)
  return(file_peaks)
}
qscoreCalculator <- function(eic){
  #Check for bogus EICs
  if(nrow(eic)<5){
    return(0)
  }
  #Calculate where each rt would fall on a beta dist (accounts for missed scans)
  scaled_rts <- (eic$rt-min(eic$rt))/(max(eic$rt)-min(eic$rt))
  
  # Create a couple different skews and test fit
  maybe_skews <- c(2.5,3,4,5) #Add 7 to catch more multipeaks and more noise
  #Add 2 to catch very slopey peaks and more noise
  best_skew <- maybe_skews[which.max(sapply(maybe_skews, function(x){
    cor(dbeta(scaled_rts, shape1 = x, shape2 = 5), eic$int)
  }))]
  perf_peak <- dbeta(scaled_rts, shape1 = best_skew, shape2 = 5)
  peak_cor <- cor(perf_peak, eic$int)
  
  
  #Calculate the normalized residuals
  residuals <- eic$int/max(eic$int)-perf_peak/max(perf_peak)
  #Calculate the minimum SD, after normalizing for any shape discrepancy
  old_res_sd <- sd(residuals)
  norm_residuals <- diff(residuals)
  new_res_sd <- sd(norm_residuals)
  while(new_res_sd<old_res_sd){
    old_res_sd <- new_res_sd
    norm_residuals <- diff(residuals)
    new_res_sd <- sd(residuals)
  }
  #Calculate SNR
  SNR <- (max(eic$int)-min(eic$int))/sd(norm_residuals*max(eic$int))
  #Return the quality score
  return(round(SNR*peak_cor^4*log10(max(eic$int))))
}
qscoreCalculator <- compiler::cmpfun(qscoreCalculator)
checkPeakCor <- function(mass, rtmin, rtmax, init_eic, 
                         file_data_table, pmppm, trapz){
  given_eic <- file_data_table[mz%between%pmppm(mass, ppm = 5) & 
                                 rt%between%c(rtmin, rtmax)]
  if(nrow(given_eic)<5){
    return(c(0, 0))
  }
  merged_eic <- merge(init_eic, given_eic, by="rt")
  if(nrow(merged_eic)<5){
    return(c(0, 0))
  }
  
  # plot(merged_eic$rt, scales::rescale(merged_eic$int.y))
  # lines(init_eic$rt, scales::rescale(merged_eic$int.x))
  
  peak_match <- cor(merged_eic$int.x, merged_eic$int.y)
  peak_area <- trapz(merged_eic$rt, merged_eic$int.y)
  return(c(peak_match, peak_area))
}
checkPeakCor <- compiler::cmpfun(checkPeakCor)
trapz <- function(x, y) {
  y <- y[order(x)]
  x <- x[order(x)]
  m <- length(x)
  xp <- c(x, x[m:1])
  yp <- c(numeric(m), y[m:1])
  n <- 2*m
  p1 <- sum(xp[1:(n-1)]*yp[2:n]) + xp[n]*yp[1]
  p2 <- sum(xp[2:n]*yp[1:(n-1)]) + xp[1]*yp[n]
  
  return(0.5*(p1-p2))
}
trapz <- compiler::cmpfun(trapz)

fillChromPeaks_wkumler <- function(object, param){
  msLevel <- 1L
  BPPARAM <- bpparam()
  
  if (length(msLevel) != 1) 
    stop("Can only perform peak filling for one MS level at a time")
  if (!hasFeatures(object, msLevel = msLevel)) 
    stop("No feature definitions for MS level ", 
         msLevel, " present. Please run 'groupChromPeaks' first.")
  if (xcms:::.hasFilledPeaks(object))
    message("Filled peaks already present, adding still missing", 
            " peaks.")
  if (xcms:::hasChromPeaks(object) & !xcms:::.has_chrom_peak_data(object)) 
    object <- updateObject(object)
  startDate <- date()
  expandMz <- expandMz(param)
  expandRt <- expandRt(param)
  fixedMz <- fixedMz(param)
  fixedRt <- fixedRt(param)
  ppm <- ppm(param)
  message("Defining peak areas for filling-in .", 
          appendLF = FALSE)
  fdef <- featureDefinitions(object, msLevel = msLevel)
  aggFunLow <- median
  aggFunHigh <- median
  tmp_pks <- chromPeaks(object)[, c("rtmin", "rtmax", 
                                    "mzmin", "mzmax")]
  pkArea <- do.call(rbind, lapply(fdef$peakidx, function(z) {
    pa <- c(aggFunLow(tmp_pks[z, 1]), aggFunHigh(tmp_pks[z, 2]), 
            aggFunLow(tmp_pks[z, 3]), aggFunHigh(tmp_pks[z, 4]))
    if (ppm != 0) {
      mzmean <- mean(pa[3:4])
      tittle <- mzmean * (ppm/2)/1e+06
      if ((pa[4] - pa[3]) < (tittle * 2)) {
        pa[3] <- mzmean - tittle
        pa[4] <- mzmean + tittle
      }
    }
    if (expandRt != 0) {
      diffRt <- (pa[2] - pa[1]) * expandRt/2
      pa[1] <- pa[1] - diffRt
      pa[2] <- pa[2] + diffRt
    }
    if (expandMz != 0) {
      diffMz <- (pa[4] - pa[3]) * expandMz/2
      pa[3] <- pa[3] - diffMz
      pa[4] <- pa[4] + diffMz
    }
    if (fixedMz != 0) {
      pa[3] <- pa[3] - fixedMz
      pa[4] <- pa[4] + fixedMz
    }
    if (fixedRt != 0) {
      pa[1] <- pa[1] - fixedRt
      pa[2] <- pa[2] + fixedRt
    }
    pa
  }))
  rm(tmp_pks)
  message(".", appendLF = FALSE)
  colnames(pkArea) <- c("rtmin", "rtmax", "mzmin", "mzmax")
  pkArea <- cbind(group_idx = 1:nrow(pkArea), pkArea, mzmed = as.numeric(fdef$mzmed))
  pkGrpVal <- featureValues(object, value = "index", msLevel = msLevel)
  message(".", appendLF = FALSE)
  if (!any(is.na(rowSums(pkGrpVal)))) {
    message("No missing peaks present.")
    return(object)
  }
  message(".", appendLF = FALSE)
  objectL <- vector("list", length(fileNames(object)))
  pkAreaL <- objectL
  req_fcol <- requiredFvarLabels("OnDiskMSnExp")
  min_fdata <- fData(object)[, req_fcol]
  rt_range <- range(pkArea[, c("rtmin", "rtmax")])
  if (hasAdjustedRtime(object)) 
    min_fdata$retentionTime <- adjustedRtime(object)
  for (i in 1:length(fileNames(object))) {
    fd <- min_fdata[min_fdata$fileIdx == i, ]
    fd$fileIdx <- 1L
    objectL[[i]] <- new("OnDiskMSnExp", 
                        processingData = new("MSnProcess", 
                                             files = fileNames(object)[i]), 
                        featureData = new("AnnotatedDataFrame", fd), 
                        phenoData = new("NAnnotatedDataFrame", 
                                        data.frame(sampleNames = "1")), 
                        experimentData = new("MIAPE", 
                                             instrumentManufacturer = "a", 
                                             instrumentModel = "a", 
                                             ionSource = "a", 
                                             analyser = "a", 
                                             detectorType = "a"))
    pkAreaL[[i]] <- pkArea[is.na(pkGrpVal[, i]), , drop = FALSE]
  }
  rm(pkGrpVal)
  rm(pkArea)
  rm(min_fdata)
  message(" OK\nStart integrating peak areas from original files")
  ph <- processHistory(object, type = xcms:::.PROCSTEP.PEAK.DETECTION)
  findPeakMethod <- "unknown"
  mzCenterFun <- "wMean"
  if (length(ph)) {
    if (is(ph[[1]], "XProcessHistory")) {
      prm <- ph[[1]]@param
      findPeakMethod <- .param2string(prm)
      if (.hasSlot(prm, "mzCenterFun")) 
        mzCenterFun <- prm@mzCenterFun
    }
  }
  cp_colnames <- colnames(chromPeaks(object))
  mzCenterFun <- paste("mzCenter", gsub(mzCenterFun, 
                                        pattern = "mzCenter.", replacement = "", 
                                        fixed = TRUE), sep = ".")
  if (findPeakMethod == "MSW") {
    rts <- rtime(object, bySample = TRUE)
    if (any(lengths(rts) > 1)) 
      stop("The data is supposed to be direct injection data, ", 
           "but I got files with more than one spectrum/", 
           "retention time!")
    res <- bpmapply(FUN = xcms:::.getMSWPeakData, objectL, pkAreaL, 
                    as.list(1:length(objectL)), MoreArgs = list(cn = cp_colnames), 
                    BPPARAM = BPPARAM, SIMPLIFY = FALSE)
  } else if (findPeakMethod == "matchedFilter") {
    res <- bpmapply(FUN = xcms:::.getChromPeakData_matchedFilter, 
                    objectL, pkAreaL, as.list(1:length(objectL)), 
                    MoreArgs = list(cn = cp_colnames, param = prm, 
                                    msLevel = msLevel), 
                    BPPARAM = BPPARAM, 
                    SIMPLIFY = FALSE)
  } else {
    res <- bpmapply(FUN = xcms:::.getChromPeakData, objectL, 
                    pkAreaL, as.list(1:length(objectL)), 
                    MoreArgs = list(cn = cp_colnames, mzCenterFun = mzCenterFun, 
                                    msLevel = msLevel), 
                    BPPARAM = BPPARAM, SIMPLIFY = FALSE)
  }
  rm(objectL)
  res <- do.call(rbind, res)
  res <- cbind(res, group_idx = unlist(lapply(pkAreaL, function(z){
    z[, "group_idx"]}), use.names = FALSE))
  res <- res[!is.na(res[, "into"]), , drop = FALSE]
  if (nrow(res) == 0) {
    warning("Could not integrate any signal for the missing ", 
            "peaks! Consider increasing 'expandMz' and 'expandRt'.")
    return(object)
  }
  rm(pkAreaL)
  gc()
  newFd <- new("MsFeatureData")
  newFd@.xData <- xcms:::.copy_env(object@msFeatureData)
  object@msFeatureData <- new("MsFeatureData")
  incr <- nrow(chromPeaks(newFd))
  for (i in unique(res[, "group_idx"])) {
    fdef$peakidx[[i]] <- c(fdef$peakidx[[i]], (which(res[, "group_idx"] == i) + incr))
  }
  fdef <- rbind(fdef, featureDefinitions(newFd)[featureDefinitions(newFd)$ms_level != 
                                                  msLevel, , drop = FALSE])
  if (!any(colnames(fdef) == "ms_level")) {
    fdef$ms_level <- 1L
  } else {fdef <- fdef[order(fdef$ms_level), ]}
  reg_extract <- rownames(chromPeaks(newFd))
  maxId <- max(as.numeric(regmatches(reg_extract,regexpr("\\d+$", reg_extract))))
  if(is.na(maxId)){
    message("Grepping for max compound number dun fuxed, plz fix")
    browser()
  }
  if (maxId < 1)
    stop("chromPeaks matrix lacks rownames; please update ", 
         "'object' with the 'updateObject' function.")
  toId <- maxId + nrow(res)
  rownames(res) <- sprintf(paste0("CP", "%0", ceiling(log10(toId + 1L)), "d"), 
                           (maxId + 1L):toId)
  chromPeaks(newFd) <- rbind(chromPeaks(newFd), res[, -ncol(res)])
  cpd <- chromPeakData(newFd)[rep(1L, nrow(res)), , drop = FALSE]
  cpd[, ] <- NA
  cpd$ms_level <- as.integer(msLevel)
  cpd$is_filled <- TRUE
  if (!any(colnames(chromPeakData(newFd)) == "is_filled")) 
    chromPeakData(newFd)$is_filled <- FALSE
  chromPeakData(newFd) <- rbind(chromPeakData(newFd), cpd)
  rownames(chromPeakData(newFd)) <- rownames(chromPeaks(newFd))
  featureDefinitions(newFd) <- fdef
  lockEnvironment(newFd, bindings = TRUE)
  object@msFeatureData <- newFd
  ph <- xcms:::XProcessHistory(param = param, date. = startDate, 
                               type. = xcms:::.PROCSTEP.PEAK.FILLING, 
                               fileIndex = 1:length(fileNames(object)), 
                               msLevel = msLevel)
  object <- xcms:::addProcessHistory(object, ph)
  object
}
formula2elements <- function(formula_vec){
  split_formulas <- formula_vec %>%
    gregexpr(pattern = "[A-Z][a-z]*[0-9]*") %>%
    regmatches(x = formula_vec)
  elements_in <- lapply(split_formulas, gsub, pattern = "[0-9]", replacement = "")
  element_counts <- lapply(split_formulas, gsub, pattern = "[A-z]", replacement = "")
  element_counts <- lapply(element_counts, function(x){
    x[!nchar(x)]<-1
    return(as.numeric(x))
  })
  mapply(`names<-`, element_counts, elements_in, SIMPLIFY = FALSE)
}
mgf_maker <- function(feature_msdata, ms1, ms2, output_file){
  ms1_to_write <- paste("BEGIN IONS",
                         paste0("PEPMASS=", feature_msdata$mzmed),
                         "MSLEVEL=1",
                         "CHARGE=1+",
                         paste0(apply(ms1, 1, paste, collapse=" "), collapse = "\n"),
                         "END IONS",
                         "", sep = "\n")
  
  if(!nrow(ms2)){
    outtext <- ms1_to_write
  } else {
    ms2_to_write <- lapply(unique(ms2$filename), function(volt_i){
      paste0(c("BEGIN IONS",
             paste0("PEPMASS=", feature_msdata$mzmed),
             "MSLEVEL=2",
             "CHARGE=1+",
             paste0(apply(ms2[ms2$filename==volt_i, c("fragmz", "int")], 
                   1, paste, collapse=" "), collapse="\n"),
             "END IONS",
             ""), collapse="\n")
    })
    
    outtext <- paste0(ms1_to_write, "\n", paste0(unlist(ms2_to_write), collapse = "\n"))
  }
  writeLines(outtext, con = output_file)
}

addH <- function(formula_vec){
  gsub(x = formula_vec, pattern = "(?<=H)\\d+", as.numeric(str_extract(formula_vec, "(?<=H)\\d+"))+1, perl = TRUE)
}
minH <- function(formula_vec){
  gsub(x = formula_vec, pattern = "(?<=H)\\d+", as.numeric(str_extract(formula_vec, "(?<=H)\\d+"))-1, perl = TRUE)
}

grabSingleFileData <- function(filename){
  msdata <- mzR:::openMSfile(filename)
  fullhd <- mzR::header(msdata)
  spectra_list <- lapply(seq_len(nrow(fullhd)), function(x){
    given_peaks <- mzR::peaks(msdata, x)
    rtime <- fullhd[x, "retentionTime"]
    return(cbind(rtime, given_peaks))
  })
  all_data <- `names<-`(as.data.frame(do.call(rbind, spectra_list)), 
                        c("rt", "mz", "int"))
  return(all_data)
}
pmppm <- RaMS::pmppm
