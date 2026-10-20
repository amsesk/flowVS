
###############################################################################
###############################################################################
## Identify denasity peaks in 1D distribution of a marker
## Input:
##    dat: a vector representing a column of a flowFrame 
## output:
##       A matrix of size px2 where p is the number of identified peaks
##       each row of the matrix indicates the location of the modes of the regions 
##       in the density estimates.
###############################################################################
###############################################################################

densityPeaks <- function(y, stain, bwFac = 2, plot=FALSE, save.plot=FALSE, populationQuant=.01, borderQuant=0.001, 
                      peak.density.thr=0.05, peak.distance.thr=0.05, annotation="",save.folder="", diagEnv=NULL, ...)
{
  ## some type checking first
  checkClass(y, "numeric")
  checkClass(stain, "character", 1)
  #if(!stain %in% colnames(x))
  #  stop("'", stain,"' is not a valid parameter in this flowFrame")
  checkClass(plot, "logical", 1)
  checkClass(borderQuant, "numeric", 1)
  
  ################################################################
  y = y[y>min(y) & y<max(y)]
  #min.range = quantile(y,.25) - 4*IQR(y)
  #max.range = quantile(y,.75) + 4*IQR(y)
  min.range = quantile(y,borderQuant) 
  max.range = quantile(y,1-borderQuant)
  y = y[y>min.range & y<max.range]
  mat = as.matrix(y)
  colnames(mat) = stain
  ff = flowFrame(mat)
  # find the picks
  curv1.filter <- filter(ff, curv1Filter(stain,bwFac=bwFac))
  bnds <- curvPeaks(curv1.filter, y)
  peaks = bnds$peaks
  ## DIAGNOSTIC stage 1: raw peaks from curvPeaks
  if(!is.null(diagEnv))
  {
    diagEnv$n_events_trimmed = length(y)
    diagEnv$stage1_rawPeaks = peaks
  }
  

  ## remove peaks with low densities , follow Ryan Brinkman's paper standard
  
#   max.peak.density = max(peaks[,'y'], na.rm = TRUE)
#   peak.rm.idx = which(peaks[,'y']/max.peak.density < peak.density.thr)
#   if(length(peak.rm.idx)>0)
#   {
#     peaks=peaks[-peak.rm.idx,]
#     if(!is(peaks, 'matrix'))
#     {
#       peaks = matrix(peaks, ncol=2)
#       colnames(peaks) = colnames(bnds$peaks)
#     }
#   }
  

  
  # remove peaks close to each other
  peak.rm.idx = vector()
  mrange = max(y, na.rm=TRUE) - min(y, na.rm=TRUE)
  #peaks = bnds$peaks
  npeaks = nrow(peaks)
  if(npeaks>=2)
  {
    p1 = 1
    p2 = 2
    while(p2<=npeaks)
    {
      peak1 = peaks[p1,]
      peak2 = peaks[p2,]
      if (abs(peak1[1]-peak2[1]) < mrange * peak.distance.thr)
      {
        peak.rm.idx = ifelse(peak1[2]<peak2[2], c(peak.rm.idx,p1), c(peak.rm.idx,p2))
        p1 = ifelse(peak.rm.idx==p2, p1, p2)
        p2 = p2 + 1
      }else
      {
        p1 = p2
        p2 = p2+1
      }
    }
  }
  if(length(peak.rm.idx)>0)
  {
    peaks=peaks[-peak.rm.idx,]
    if(!is(peaks, 'matrix'))
    {
      peaks = matrix(peaks, ncol=2)
      colnames(peaks) = colnames(bnds$peaks)
    }
  }
  ###########################################################
  
  
  dens <- density(y)
  peaksStats = matrix(nrow=0, ncol=4)
  colnames(peaksStats) = c('n', 'mean', 'median', 'variance')

  ## DIAGNOSTIC stage 2: peaks after removal of close peaks
  if(!is.null(diagEnv))
  {
    diagEnv$removed_close_peaks = peak.rm.idx
    diagEnv$stage2_filteredPeaks = peaks
    diagEnv$populations = list()
  }
  if(nrow(peaks)==0)
  {
    if(!is.null(diagEnv)) diagEnv$returned = peaksStats
    return (peaksStats)
  }
  rec = NULL
  saveRec = function() if(!is.null(diagEnv)) diagEnv$populations[[i]] = rec
  # find the minima between every consecutive pairs of peaks 
  left = min(y, na.rm=TRUE)
  bounds = vector()
  separators = vector()
  populations = list()
  for(i in 1:nrow(peaks))
  {
    p = peaks[i,]
    
    if(i==nrow(peaks))
    {
      right = max(y, na.rm=TRUE) 
    }else
    {
      sel <- (dens$x > peaks[i,1] & dens$x < peaks[i+1,1])
      right <- dens$x[sel][which.min(dens$y[sel])]
    }
    
    
    y.sel = y[y>=left & y<=right]
    min.range = quantile(y.sel,populationQuant) 
    max.range = quantile(y.sel,1-populationQuant) 
    #cat(left, ' ', right, p[1],' ',min.range,' ', max.range, ' \n')
    
    
    #population = y.sel[ y.sel>=min.range & y.sel<= max.range]
    ## DIAGNOSTIC stage 3: per-population values (saved before anything that can fail)
    rec = list(population=i, npeaks=nrow(peaks), peak=p, left=left, right=right,
               n_ysel=length(y.sel), min.range=min.range, max.range=max.range,
               r=NA_real_, branch=NA_character_, n_population=NA_integer_,
               row_appended=FALSE, rare_added=FALSE, skipped_na_peak=is.na(p[1]))
    saveRec()
    if(!is.na(p[1]))
    {
      ### check skewness in the population 
      r = (max.range-p[1])/(p[1]-min.range)
      rec$r = r
      saveRec()
      #print(r)
      if(r >=3 && i==nrow(peaks)) ##rightmost popualtion 
      {
        right = min(max.range, p[1]+3*(p[1]-min.range)  )
        rec$branch = "rightmost_skew"
        y.sel = y[y>=left & y<=right]
        min.range = quantile(y.sel,populationQuant) 
        max.range = quantile(y.sel,1-populationQuant)
      }else if(r <= 1/3 && i==1) ## leftmost population 
      {
        #print('hhh')
        left = max(min.range, p[1]-3*(max.range-p[1])  )
        rec$branch = "leftmost_skew"
        y.sel = y[y>=left & y<=right]
        min.range = quantile(y.sel,populationQuant) 
        max.range = quantile(y.sel,1-populationQuant)
      }
      population = y.sel[ y.sel>=min.range & y.sel<= max.range]
      rec$left_final = left; rec$right_final = right
      rec$min.range_final = min.range; rec$max.range_final = max.range
      rec$n_population = length(population)
      saveRec()
      if(length(population)!=0)
      {
        peaksStats = rbind(peaksStats, c(length(population), mean(population), median(population), var(population))) 
        populations[[length(populations)+1]] = population
        rec$row_appended = TRUE; rec$row = peaksStats[nrow(peaksStats),]
        saveRec()
        bounds = c(bounds, min.range,max.range)
        if(i< nrow(peaks)) separators = c(separators, right)   
      }
      
      left = right
      
      ## look for potential rare populations at the right side 
      if(r >=5 && i==nrow(peaks) && length(separators)==0) ## look for rare populaiton in the right 
      {
        rare = y[y>=left]
        min.range = quantile(rare,populationQuant) 
        max.range = quantile(rare,1-populationQuant)
        
        rare = rare[ rare>=min.range & rare<= max.range]
        if(length(rare)!=0)
        {
          populations[[length(populations)+1]] = rare
          rec$rare_added = TRUE; rec$rare_n = length(rare)
          saveRec()
          separators = c(separators, left)
          bounds = c(bounds, min.range,max.range)
        }
      }
    }
    
  }
  
  if(!is.null(diagEnv))
  {
    diagEnv$separators = separators
    diagEnv$bounds = bounds
    diagEnv$returned = peaksStats
  }
  if(plot && nrow(peaksStats)>0)
  {
    main.text = paste("Peaks for", stain)
    if(annotation!="")
      main.text = paste("Sample=", annotation," : " ,main.text, sep="")
    plot(dens, main=main.text, cex.main=1,...)
    
    ramp <- colorRamp(c('blue', "green", "red"))
    cols = rgb( ramp(seq(0, 1, length = length(bounds)/2)), alpha=50, maxColorValue = 255)
    
    for(i in seq(1, length(bounds), 2) ){
      sel <- dens$x >= bounds[i] & dens$x <= bounds[i+1]
      polygon(c(dens$x[min(which(sel))], dens$x[sel],
                dens$x[max(which(sel))]), c(0, dens$y[sel], 0),
              col=cols[(i+1)/2], border=NA)
    }
    lines(dens, ...)
    for(sep in separators)
      abline(v = sep, col = 2, lwd = 2)
    
    #legend("topright", c("breakpoints", "dens regions"),
    #       fill=c("red", "lightgray"), bty="n")
    if(save.plot==TRUE && annotation!="")
    {
      
      file.name = paste(stain, "_", annotation, ".pdf", sep="")
      if(save.folder!="")
      {
        dir.create(save.folder, showWarnings=FALSE)
        file.name = paste(save.folder, file.name, sep="/")
      }
      dev.copy2pdf(file=file.name)
    }
  }
  return(peaksStats)
}



# copied from flowStats to remove accessing via flowStats:::curvPeaks
curvPeaks <- function(x, dat, borderQuant=0.01, n=201, from, to, densities=NULL)
{
    ## Some type-checking first
    checkClass(x, "multipleFilterResult")
    checkClass(dat, c("numeric","NULL"))
    checkClass(borderQuant, "numeric", 1)
    checkClass(n, "numeric", 1)
    if(missing(from))
    from <- min(dat)
    checkClass(from, "numeric", 1)
    if(missing(to))
    to <- max(dat)
    checkClass(to, "numeric", 1)
    ## extract boundaries
    bound <- attr(x@subSet, "boundaries")
    border_lower <- quantile(c(from,to), borderQuant)
    border_upper <- quantile(c(from,to), 1-borderQuant)
    
    dens <- if(!is.null(densities)) densities else
    density(dat, n=n, from=from, to=to, na.rm=TRUE)$y
    ## iterate over regions
    regPoints <- list()
    peaks <- midpoints <- regions <- densFuns <- NULL
    i <- 1
    if(!all(is.na(bound[[1]]))){
        #oo <- options(warn=-1)
        #on.exit(options(oo))
        for(b in bound){
            ## discard regions that lie entirely within a margin
            if(b[2] > border_lower && b[1] < border_upper){
                ## clip regions that extend past the density grid [from, to]; approxfun returns
                ## NA outside the grid, which makes optimize() return an NA objective and an
                ## out-of-range maximum
                b[1] <- max(b[1], from)
                b[2] <- min(b[2], to)
                ## approximate density by function
                afun <- approxfun(seq(from, to, len=n), dens)
                sel <- seq(b[1], b[2], len=50)
                regPoints[[i]] <- cbind(x=sel, y=afun(sel))
                ## compute maximum of function
                m <- optimize(afun, b, maximum=TRUE)
                peaks <- rbind(peaks, cbind(x=m$maximum, y= m$objective))
                regions <- rbind(regions, cbind(left=b[1], right=b[2]))
                midpoints <- c(midpoints, mean(b))
                densFuns <- c(densFuns, afun)
                i <- i+1
            }
        }
    }
    if(i==1)
    return(list(peaks=cbind(x=NA, y=NA), regions=cbind(left=NA, right=NA),
    midpoints=NA, regPoints=list(cbind(x=NA, y=NA)),
    densFuns=NA))
    return(list(peaks=peaks, regions=regions, midpoints=midpoints,
    regPoints=regPoints, densFuns=densFuns))
}



## ===========================================================================
##  Copied from flowCore package
## ---------------------------------------------------------------------------
## Check for the class of object x and its length and cast error if wrong
checkClass <- function(x, class, length=NULL, verbose=FALSE,
mandatory=TRUE)
{
    if(mandatory && missing(x))
    stop("Argument '", substitute(x), "' missing with no default",
    call.=verbose)
    msg <- paste("'", substitute(x), "' must be object of class ",
    paste("'", class, "'", sep="", collapse=" or "), sep="")
    fail <- !any(sapply(class, function(c, y) is(y, c), x))
    if(!is.null(length) && length(x) != length)
    {
        if(!is.null(x))
        {
            fail <- TRUE
            msg <- paste(msg, "of length", length)
        }
    }
    if(fail) stop(msg, call.=verbose) else invisible(NULL)
}
