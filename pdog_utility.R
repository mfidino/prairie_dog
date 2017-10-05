###################
# utility functions
###################

# loads packages

package_load<-function(packages = NA, quiet=TRUE, verbose=FALSE, 
  warn.conflicts=FALSE){
  
  # download required packages if they're not already
  pkgsToDownload<- packages[!(packages  %in% installed.packages()[,"Package"])]
  if(length(pkgsToDownload)>0)
    install.packages(pkgsToDownload, repos="http://cran.us.r-project.org", quiet=quiet, verbose=verbose)
  
  # then load them
  for(i in 1:length(packages))
    require(packages[i], character.only=T, quietly=quiet, warn.conflicts=warn.conflicts)
}

pull_year <- function(x){
  strsplit(as.character(x$variable), "\\.") %>% 
    sapply(., "[[", 2) %>% 
    as.numeric(.) %>% 
    return(.)
}

# calculates pc from saura and rubio

lose1 <- function(network, PC = NULL, areas = NULL, tail_distance = NULL, 
  sq.m = NULL){
  y <- V(network)
  
  cores <- detectCores()-2
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  alpha <- log(0.05)/tail_distance
  pk <- rep(0, length(y))
  
  foreach(k = 1:length(y), .combine = 'c') %do% {
    tmp <- network - y[k]
    tmp2 <- distances(tmp, V(tmp), 
                     weights = E(tmp)$weight, to = V(tmp))
    diag(tmp2) <- 1
    pcnum <- exp(-abs(alpha)*tmp2) * tcrossprod(areas[-k])
    pcnum <- sum(pcnum) / (sq.m^2)
    pk[k] <- 100 * ((PC - pcnum)/PC)
    pk[k]
    
  }
  stopCluster(cl)
  return(pk)
  
}

calc_pck <- function(fragments = NULL, sq.m = 347000000, 
  cut_connections_at = 2000, tail_distance = 5000 ){
  
  ppd <- fragments %>% select(one_of(c("easting", "northing"))) %>% 
    dist(, diag = TRUE, upper = TRUE) %>% as.matrix
  ppd[ppd>cut_connections_at] <- 0 
  colnames(ppd) <- fragments$FRAG.ID
  
  # make the fragment network
  frag_network <- graph_from_adjacency_matrix(ppd, weighted = TRUE, 
    mode = "undirected")
  
  #V(frag_network)$size = pk*5
  #node_size = setNames(pk * 5, fragments$FRAG.ID)
  #plot(frag_network, layout = fragments[,7:8], xlim = range(fragments$easting),
  #  ylim = range(fragments$northing), rescale = FALSE, vertex.label = NA, 
  #  edge.arrow.size = 0.0)
  
  # calculate topological distances between sites
  topo_distance <- distances(frag_network, V(frag_network), 
    weights = E(frag_network)$weight, to = V(frag_network))
  # sites are connected to themselves
  diag(topo_distance) <- 1
  
  # calculate tail distnace
  td <- log(0.05)/tail_distance
  # calculate PCnum as in saura and rubio 2010
  PCnum <- exp(-abs(td)*topo_distance) * tcrossprod(fragments$area)
  
  # calculate PC from PCnum
  PC <- (sum(PCnum)) / (sq.m^2)
  
  # calculate pck
  pck <- lose1(frag_network, PC = PC, tail_distance = tail_distance,
    areas = fragments$area, sq.m = sq.m)
  
  return(pck)
}


calc_distances <- function(past_fragments = NULL, 
  past_pdogs = NULL, current_year = NULL){
  
  nearest_dogs <- nearest_frag <- aw_dogs <- aw_frag <- rep(0, 
    nrow(current_year))
  
  for(j in 1:nrow(current_year)){
    # make matrix to pdogs
    to_dist_pd <- rbind(current_year[j,], past_pdogs) %>% 
      filter(!duplicated(FRAG.ID))
    # make matrix to fragments
    to_dist_frag <- rbind(current_year[j,], past_fragments) %>% 
      filter(!duplicated(FRAG.ID))
    # distance matrix pdogs
    my_dist_pd <- to_dist_pd %>% select(one_of(c("easting", "northing"))) %>% 
      dist(, diag = TRUE, upper = TRUE) %>% as.matrix
    # distance matrix fragments
    my_dist_frag <- to_dist_frag %>% select(one_of(c("easting", "northing"))) %>% 
      dist(, diag = TRUE, upper = TRUE) %>% as.matrix
    # nearest pdog colony
    nearest_dogs[j] <- sort(as.numeric(my_dist_pd[,1]))[2]
    # area weigthed nearest pdog colony
    aw_dogs[j] <- (prod(to_dist_pd$area[c(1,order(my_dist_pd[,1])[2])])^0.7)/
      (nearest_dogs[j]^1.7)
    # nearest fragment
    nearest_frag[j] <- sort(as.numeric(my_dist_frag[,1]))[2]
    # area weighted nearest fragment
    aw_frag[j] <- (prod(to_dist_frag$area[c(1,order(my_dist_frag[,1])[2])])^0.7)/
      (nearest_frag[j]^1.7)
  }
  
  return(data.frame(nearest_dogs, aw_dogs, nearest_frag, aw_frag))
  
}

calc_waic <- function(x){
  p_prob <- x[,grep("p_prob", colnames(x))]
  lppd <- apply(p_prob, 2, mean) %>% log() %>% sum
  pwaic <- apply(log(p_prob), 2, var) %>% sum
  
  WAIC <- -2 * (lppd - pwaic)
  return(WAIC)
}


HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}


diagMCMC = function( codaObject , parName=varnames(codaObject)[1] ,
  saveName=NULL , saveType="jpg" ) {
  DBDAplColors = c("skyblue","black","royalblue","steelblue")
  openGraph(height=5,width=7)
  par( mar=0.5+c(3,4,1,0) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
    cex.lab=1.5 )
  layout(matrix(1:4,nrow=2))
  # traceplot and gelman.plot are from CODA package:
  require(coda)
  coda::traceplot( codaObject[,c(parName)] , main="" , ylab="Param. Value" ,
    col=DBDAplColors ) 
  tryVal = try(
    coda::gelman.plot( codaObject[,c(parName)] , main="" , auto.layout=FALSE , 
      col=DBDAplColors )
  )  
  # if it runs, gelman.plot returns a list with finite shrink values:
  if ( class(tryVal)=="try-error" ) {
    plot.new() 
    print(paste0("Warning: coda::gelman.plot fails for ",parName))
  } else { 
    if ( class(tryVal)=="list" & !is.finite(tryVal$shrink[1]) ) {
      plot.new() 
      print(paste0("Warning: coda::gelman.plot fails for ",parName))
    }
  }
  DbdaAcfPlot(codaObject,parName,plColors=DBDAplColors)
  DbdaDensPlot(codaObject,parName,plColors=DBDAplColors)
  mtext( text=parName , outer=TRUE , adj=c(0.5,0.5) , cex=2.0 )
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,"Diag",parName), type=saveType)
  }
}


DbdaAcfPlot = function( codaObject , parName=varnames(codaObject)[1] , plColors=NULL ) {
  if ( all( parName != varnames(codaObject) ) ) { 
    stop("parName must be a column name of coda object")
  }
  nChain = length(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  xMat = NULL
  yMat = NULL
  for ( cIdx in 1:nChain ) {
    acfInfo = acf(codaObject[,c(parName)][[cIdx]],plot=FALSE) 
    xMat = cbind(xMat,acfInfo$lag)
    yMat = cbind(yMat,acfInfo$acf)
  }
  matplot( xMat , yMat , type="o" , pch=20 , col=plColors , ylim=c(0,1) ,
    main="" , xlab="Lag" , ylab="Autocorrelation" )
  abline(h=0,lty="dashed")
  EffChnLngth = effectiveSize(codaObject[,c(parName)])
  text( x=max(xMat) , y=max(yMat) , adj=c(1.0,1.0) , cex=1.25 ,
    labels=paste("ESS =",round(EffChnLngth,1)) )
}

DbdaDensPlot = function( codaObject , parName=varnames(codaObject)[1] , plColors=NULL ) {
  if ( all( parName != varnames(codaObject) ) ) { 
    stop("parName must be a column name of coda object")
  }
  nChain = length(codaObject) # or nchain(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  xMat = NULL
  yMat = NULL
  hdiLims = NULL
  for ( cIdx in 1:nChain ) {
    densInfo = density(codaObject[,c(parName)][[cIdx]]) 
    xMat = cbind(xMat,densInfo$x)
    yMat = cbind(yMat,densInfo$y)
    hdiLims = cbind(hdiLims,HDIofMCMC(codaObject[,c(parName)][[cIdx]]))
  }
  matplot( xMat , yMat , type="l" , col=plColors , 
    main="" , xlab="Param. Value" , ylab="Density" )
  abline(h=0)
  points( hdiLims[1,] , rep(0,nChain) , col=plColors , pch="|" )
  points( hdiLims[2,] , rep(0,nChain) , col=plColors , pch="|" )
  text( mean(hdiLims) , 0 , "95% HDI" , adj=c(0.5,-0.2) )
  EffChnLngth = effectiveSize(codaObject[,c(parName)])
  MCSE = sd(as.matrix(codaObject[,c(parName)]))/sqrt(EffChnLngth) 
  text( max(xMat) , max(yMat) , adj=c(1.0,1.0) , cex=1.25 ,
    paste("MCSE =\n",signif(MCSE,3)) )
}


openGraph = function( width=7 , height=7 , mag=1.0 , ... ) {
  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
    tryInfo = try( X11( width=width*mag , height=height*mag , type="cairo" , 
      ... ) )
    if ( class(tryInfo)=="try-error" ) {
      lineInput = readline("WARNING: Previous graphics windows will be closed because of too many open windows.\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
      graphics.off() 
      X11( width=width*mag , height=height*mag , type="cairo" , ... )
    }
  } else { # Windows OS
    tryInfo = try( windows( width=width*mag , height=height*mag , ... ) )
    if ( class(tryInfo)=="try-error" ) {
      lineInput = readline("WARNING: Previous graphics windows will be closed because of too many open windows.\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
      graphics.off() 
      windows( width=width*mag , height=height*mag , ... )    
    }
  }
}

saveGraph = function( file="saveGraphOutput" , type="pdf" , ... ) {
  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
    if ( any( type == c("png","jpeg","jpg","tiff","bmp")) ) {
      sptype = type
      if ( type == "jpg" ) { sptype = "jpeg" }
      savePlot( file=paste0(file,".",type) , type=sptype , ... )     
    }
    if ( type == "pdf" ) {
      dev.copy2pdf(file=paste0(file,".",type) , ... )
    }
    if ( type == "eps" ) {
      dev.copy2eps(file=paste0(file,".",type) , ... )
    }
  } else { # Windows OS
    file=paste0(file,".",type) 
    savePlot( file=file , type=type , ... )
  }
}
