removeVertexTerm <- function(GG,NAME){
  
  if( !is.null(get.vertex.attribute(GG,NAME)) ){
    GG <- remove.vertex.attribute(GG,name=NAME)
  }
  
  if( !is.null(get.vertex.attribute(GG,gsub("_","",NAME))) ){
    GG <- remove.vertex.attribute(GG,name=gsub("_","",NAME))
  }
  
  return(GG)
  
}

getClustering<-function(gg,alg=c('lec','wt','fc','infomap','louvain','sgG1','sgG2','sgG5','spectral')){
  alg <- match.arg(alg)
  lec<-function(gg){
    lec     <- igraph::leading.eigenvector.community(gg)
    ll      <- igraph::leading.eigenvector.community(gg, start=membership(lec))
  }
  cl<-switch(alg,
             lec=lec(gg),
             wt=igraph::walktrap.community(gg),
             fc=igraph::fastgreedy.community(gg),
             infomap=igraph::cluster_infomap(gg),
             louvain=igraph::cluster_louvain(gg),
             sgG1=igraph::spinglass.community(gg,
                                              spins=as.numeric(500),gamma=1),
             sgG2=igraph::spinglass.community(gg,
                                              spins=as.numeric(500),gamma=2),
             sgG5=igraph::spinglass.community(gg,
                                              spins=as.numeric(500),gamma=5),
             spectral=rSpectral::spectral_igraph_communities(gg,
                                                             Cn_min=5,fix_neig = 1)
  )
  return(cl)
}

fSemilocal<-function(gg){
  N    <- length(V(gg)$name)
  meas <- matrix(0, nrow=N, ncol=3)
  meas[,1]<-ego_size(gg,order = 2,mode='all')-1
  neigSum<-function(i,graph,vec){
    neig <- igraph::neighbors(graph,v=i,mode="all")
    return(sum(vec[neig]))
  }
  meas[,2]<-sapply(1:N,neigSum,graph=gg,vec=meas[,1])
  meas[,3]<-sapply(1:N,neigSum,graph=gg,vec=meas[,2])
  return(as.numeric(meas[,3]))
}


applpMatrixToGraph<-function(gg,m){
  ggm<-gg
  measures<-colnames(m)
  id.col<-which(measures=='ID')
  meas.col<-which(measures!='ID')
  for(i in meas.col){
    #remove previous annotation of that name
    #check does it really needed
    ggm<-removeVertexTerm(ggm,measures[i])
    idx<-match(V(gg)$name,m[,id.col])
    naid<-which(is.na(idx))
    if(length(naid)==0){
      ggm<-set.vertex.attribute(graph=ggm,
                                name=measures[i],
                                index = V(ggm),
                                value = m[idx,i])
    }else{
      gindex<-which(is.na(idx))
      ggm<-set.vertex.attribute(graph=ggm,
                                name=measures[i],
                                index = gindex,
                                value = m[idx[gindex],i])
    }
  }
  return(ggm)
}

calShorestPaths <- function(gg){
  
  N    <- length(V(gg)$name)
  meas <- matrix(0, nrow=N, ncol=3)
  
  for( i in 1:N ){
    sp <- as.numeric(igraph::shortest.paths(gg,i))
    sp <- sp[-i]
    sp <- sp[!sp == Inf]
    meas[i,1] <- min(sp)
    meas[i,2] <- round(mean(sp),3)
    meas[i,3] <- round(sd(sp),3)
  }
  
  return(meas)
  
}

getCentralityMatrix<-function(gg){
  ID <- V(gg)$name
  N  <- length(ID)
  
  CN  <- c("ID","DEG","BET","CC","SL","mnSP","PR","sdSP")
  
  tmp <- matrix(0,nrow=N,ncol=length(CN))
  colnames(tmp) <- CN
  
  tmp[,1] <- ID
  tmp[,2] <- as.vector(igraph::degree(graph=gg))
  tmp[,3] <- as.character(round(betweenness(gg),3))
  tmp[,4] <- as.character(round(transitivity(gg,"local"),3))
  sl<- fSemilocal(gg)
  tmp[,5] <- as.character(round(sl,3))
  
  res <- calShorestPaths(gg)
  tmp[,6]  <- as.character(res[,2])
  tmp[,7]  <- as.character(round(as.vector(page.rank(graph=gg,vids=V(gg),directed=F,options=igraph.arpack.default)$vector),6))
  tmp[,8]  <- as.character(res[,3])
  
  return(tmp)
  
}

calcCentrality<-function(gg){
  m<-getCentralityMatrix(gg)
  ggm<-applpMatrixToGraph(gg,m)
  return(ggm)
}


buildConsensusMatrix<-function(lcc){
  N        = NULL
  I        = NULL
  M        = NULL
  C        = NULL
  initM    = TRUE
  
  NJobs    = 0
  max_com  = 0
  min_com  = 500
  for(i in 1:length(lcc)){
    tb<-lcc[[i]]
    ## make sure node id == -1 if node com == -1
    indx = tb[,3] == -1
    tb[indx,2] = -1
    
    if(initM){
      N     = dim(tb)[1]
      temp  = matrix(0,nrow=N,ncol=N)
      I     = temp
      M     = temp
      initM = FALSE
      rm(temp)
    }
    
    # k.coms    = tb[tb[,3] != -1,3]
    # k.max     = max(k.coms,na.rm=T)
    # k.min     = min(k.coms,na.rm=T)
    #
    # if( k.max > max_com   ){ max_com   = k.max; }
    # if( k.min < min_com   ){ min_com   = k.min; }
    
    study = calculateConsensusMat( data=tb )
    
    I = I + study$I
    M = M + study$M
    
    NJobs = NJobs + 1;
    
    rm(tb,study)
  }
  if( !is.null(N) ){
    C = do.call( cbind, lapply(1:N, function(s) matrixDiv(M[,s],I[,s])))
  }
  return(C)
}

matrixDiv <- function(x,y){
  
  N    = length(x)
  res  = rep(0,N)
  indx = y != 0
  if( sum(indx) != 0 ){
    res[indx] = x[indx]/y[indx]
  }
  
  return(res)
}

sampleGraphClust<-function(gg,mask,alg,type,reclust=FALSE,Cnmin=-1,Cnmax=10){
  IDS <- V(gg)$name;
  ids <- V(gg)$name;
  
  #---subsampling scheme
  if( type == 1 ){
    
    nr  <- ceiling( length(E(gg))*(mask/100) )
    ggM <- delete_edges(gg,sample(E(gg),nr))
    
  }
  
  if( type == 2 ){
    
    nr  <- ceiling( length(V(gg))*(mask/100) )
    ggM <- delete_vertices(gg,sample(V(gg),nr))
    
  }
  
  
  #---Find Largest CC
  ggLCC <- findLCC(ggM)
  #---
  
  
  #---build consensus file
  cc       <- matrix(-1, ncol=3, nrow=length(V(gg)))
  cc[,1]   <- V(gg)$name
  cc[,2]   <- ifelse(cc[,1] %in% V(ggLCC)$name,cc[,1],-1)
  
  
  if( Cnmin > 0 ){
    Cnmin = floor( (Cnmin*length(V(ggLCC)))/100 )
  } else {
    Cnmin = 1;
  }
  if( Cnmax > 0 ){
    Cnmax = floor( (Cnmax*length(V(ggLCC)))/100 )
  } else {
    #---default is 10% of network size
    Cnmax = floor( (10*length(V(ggLCC)))/100 )
  }
  cl<-getClustering(ggLCC,alg)
  if(reclust){
    ggLCC = igraph::set.vertex.attribute(ggLCC,alg,V(ggLCC),cl$membership)
    
    oo = recluster( ggLCC, alg, Cnmax )
    
    if( !is.null(oo) ){
      cc[,3]   <- ifelse(cc[,2] %in% oo[,1],oo[,4],-1)
    }
    
    
  }else{
    cc[,3]   <- ifelse(cc[,2] %in% cl$names,cl$membership,-1)
  }
  return(cc)
}

findLCC <- function(GG){
  
  dec <- decompose.graph(GG)
  d=1
  CC=length(V(dec[[1]]))
  for( i in 1:length(dec) ){
    if(length(V(dec[[i]])) > CC){
      d=i
      CC=length(V(dec[[i]]))
    }
  }
  
  GG  <- decompose.graph(GG)[[d]]
  return(GG)
  
}

calculateConsensusMat <- function( data=NULL ){
  
  I = NULL
  M = NULL
  
  if( !is.null(data) ){
    
    N     = dim(data)[1]
    tempI = matrix(0,nrow=N,ncol=N)
    tempM = matrix(0,nrow=N,ncol=N)
    
    for( i in 1:N ){
      comi = as.numeric(data[i,3])
      keyi = as.numeric(data[i,2])
      jj   = seq(i,N,1)
      if( keyi != -1 ){
        
        comj = as.numeric(data[jj,3])
        keyj = as.numeric(data[jj,2])
        
        ## I
        indxJ = jj[keyj!=-1]
        Nindx = length(indxJ)
        
        tempI[i,indxJ] = as.numeric(rep(1,Nindx))
        tempI[i,i]     = as.numeric(0.5)
        
        ## M
        indxC = jj[comj != -1 & comi == comj]
        Nindx = length(indxC)
        
        tempM[i,indxC] = as.numeric(rep(1,Nindx))
        tempM[i,i]     = as.numeric(0.5)
        
      }
    }
    
    M = tempM + t(tempM)
    I = tempI + t(tempI)
    
    rm(tempM,tempI)
    
  }
  
  return(list(M=M,I=I))
  
}


makeConsensusMatrix<-function(gg,N,mask,alg,type,
                              reclust=FALSE,Cnmin=-1,Cnmax=10){
  lcc<-lapply(1:N, function(.x)sampleGraphClust(gg=gg,
                                                mask=mask,
                                                alg=alg,
                                                type=type,
                                                reclust=reclust,
                                                Cnmin=Cnmin,
                                                Cnmax=Cnmax))
  mm<-buildConsensusMatrix(lcc)
  colnames(mm)<-V(gg)$name
  rownames(mm)<-V(gg)$name
  return(mm)
}


getBridgeness <- function(gg, alg,conmat) {
  #---number of vertices/genes
  N    <- length(V(gg))
  
  #---number of edges/PPIs
  M    <- length(E(gg))
  #---container column names
  if ("GeneName" %in% names(vertex.attributes(gg))) {
    CN   <- c('ENTREZ.ID', 'GENE.NAME', sprintf("BRIDGENESS.%s", alg))
    FROM <- 3
  } else{
    CN   <- c('ENTREZ.ID',  sprintf("BRIDGENESS.%s", alg))
    FROM <- 2
  }
  #---container to store Bridgeness for algorithm 'alg'
  meas <- matrix(0, nrow = N, ncol = length(CN))
  colnames(meas) <- CN
  
  meas[, 1] <- as.character(V(gg)$name)
  if ("GeneName" %in% names(vertex.attributes(gg))) {
    meas[, 2] <- as.character(V(gg)$GeneName)
  }
  ##START filling meas after PageRank column
  
  cat("calculating Bridgeness for: ", alg, "\n")
  
  ##format consensus matrix
  ##the consensus matrix you may have made (as a numeric matrix)
  # cm           <- data.frame(conmat);
  # names(cm)    <- rownames(rm);
  # rownames(cm) <- rownames(rm);
  # cm           <- as.matrix(cm);
  #cm<-conmat
  
  ##get consensus matrix indices for each edge in edgelist
  indA <- match(igraph::get.edgelist(gg)[,1],rownames(conmat))
  indB <- match(igraph::get.edgelist(gg)[,2],rownames(conmat))
  
  dat  <- data.frame(indA,indB)
  
  ##get community assigned to each vertex in edgelist from the algorithm 'alg'
  elA    <- igraph::get.vertex.attribute(gg,alg,V(gg))[match(igraph::get.edgelist(gg)[,1],V(gg)$name)]
  elB    <- igraph::get.vertex.attribute(gg,alg,V(gg))[match(igraph::get.edgelist(gg)[,2],V(gg)$name)]
  
  ##for each edge record the community assigned to each vertex and it's consensus matrix value
  ed      <- matrix(ncol=6,nrow=length(E(gg)))
  ed[,1]  <- igraph::get.edgelist(gg)[,1]
  ed[,2]  <- igraph::get.edgelist(gg)[,2]
  ed[,3]  <- elA
  ed[,4]  <- elB
  ed[,5]  <- apply(dat,1,function(x,mat) mat[x[1],x[2]], mat=conmat)
  ed[,6]  <- (as.numeric(elA)-as.numeric(elB))
  
  ##maximum number of communities found by clustering algorithm
  Cmax  <- max(as.numeric(igraph::get.vertex.attribute(gg,alg,V(gg))))
  
  ##loop over each vertex in the graph
  for( i in 1:length(V(gg)) ){
    
    ##get edges belonging to the i'th veretx
    ind <- which(ed[,1] == V(gg)$name[i] | ed[,2] == V(gg)$name[i])
    
    ##get community belonging to the i'th vertex
    c <- igraph::get.vertex.attribute(gg,alg,V(gg))[i]
    
    ##reorder edge communities, so ed[,3] equals current community no: 'c'
    for( k in 1:length(ind) ){
      if( ed[ind[k],6] != 0 && ed[ind[k],4] == c ){
        ed[ind[k],4] <- ed[ind[k],3]
        ed[ind[k],3] <- c
      }
    }
    
    ##number of communities i'th vertex is connected too (via it's edges)
    cc <- unique(ed[ind,4])
    
    ##use sum of consensus values to calculate the likelihood of i'th
    ##vertex beloning to to k'th community.
    prob <- vector(length=length(cc))
    for( k in 1:length(cc) ){
      prob[k] = sum(as.numeric(ed[which(ed[ind,4]==cc[k]),5]))/length(ind)
    }
    
    ##normalise
    prob <- prob/sum(prob)
    
    ##calculate bridgeness of i'th vertex
    ##Fuzzy communities and the concept of bridgeness in complex networks, T. Nepusz, arXiv, 2007
    b    <- sum( (prob - 1/Cmax) * (prob - 1/Cmax))
    
    Kzero <- Cmax - length(cc)
    b = b + sum(rep((1/(Cmax*Cmax)),times=Kzero))
    
    ##store values
    ##BRIDGENESS.
    meas[i,(FROM)]  <- 1-sqrt( Cmax/(Cmax-1) * b )
    
  }
  
  return(meas)
}



# helper function for graphing 2c
scale <- function(x, VALUE=NULL){
  
  x = as.numeric(as.vector(x))
  
  xmin <- min(x,na.rm=T)
  xmax <- max(x,na.rm=T)
  
  if( is.null(VALUE) ){
    
    x  <- x-xmin
    x  <- ifelse(!is.na(x), x/(xmax-xmin), NA) 
    
    return(x)
  }
  
  value = as.numeric(as.vector(value)[1])
  value = value-xmin
  value = ifelse(!is.na(value), value/(xmax-xmin), NA) 
  return(value)
}