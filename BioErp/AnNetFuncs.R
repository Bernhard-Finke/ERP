# includes functions from the package AnNet (available at https://github.com/lptolik/AnNet/ and
# developed by Anatoly Sorokin, Colin Mclean and Oksana Sorokina) used to perform some of the analysis in Figure 3.



removeVertexTerm <- function(GG,NAME){
  
  if( !is.null(get.vertex.attribute(GG,NAME)) ){
    GG <- remove.vertex.attribute(GG,name=NAME)
  }
  
  if( !is.null(get.vertex.attribute(GG,gsub("_","",NAME))) ){
    GG <- remove.vertex.attribute(GG,name=gsub("_","",NAME))
  }
  
  return(GG)
  
}

getClustering<-function(gg,alg=c('lec','wt','fc','infomap','louvain','sgG1','sgG2','sgG5')){
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
                                              spins=as.numeric(500),gamma=5)
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

calcDiseasePairs<-function(gg,name,diseases=NULL,permute=c('none','random','binned')){
  permute<-match.arg(permute)
  gda<-prepareGDA(gg,name)
  NN  <- length(which(gda!=""))
  if(is.null(diseases)){
    diseases<-getAnnotationList(gda,sort='freq')
  }else{
    remove<-c()
    diseases<-escapeAnnotation(diseases)
    for(d in 1:length(diseases)){
      if(!any(grepl(diseases[d],gda))){
        remove<-c(remove,d)
      }
    }
    if(length(remove)>0){
      diseases<-diseases[-remove]
    }
  }
  if(permute!='none'){
    rgda <- matrix(NA,nrow=vcount(gg),ncol=(length(diseases)))
    colnames(rgda) <- c(diseases)
    if(permute=='binned'){
      map <- degree.binned.GDAs(gg,gda,diseases)
    }
    
    for( d in 1:length(diseases) ){
      
      IDS <- V(gg)$name[grepl(diseases[d],gda,fixed=T)]
      N   <- length(IDS)
      
      if(permute=='random'){
        ## permute the N GDA's relative to all gene ids
        IDS <- permute( V(gg)$name, N) #case
      }else if(permute=='binned'){
        IDS <- sample.deg.binned.GDA(map,diseases[d])
      }
      rgda[match(IDS,V(gg)$name),d]<-1
    }
    gda<-apply(rgda,1,function(.x)paste(diseases[!is.na(.x)],collapse = COLLAPSE))
  }
  res           <- matrix(0 ,ncol=4, nrow=length(diseases))
  colnames(res) <- c("Disease","N","mean_ds","SD_ds")
  res[,1]       <- unescapeAnnotation(diseases)
  
  
  #--- store minimum shorest paths for each gda, and each disease
  oo <- matrix(".",nrow=length(gda),ncol=(length(diseases)+2))
  colnames(oo) <- c("Gene.ID","Gene.Name",diseases)
  oo[,1]       <- V(gg)$name#[gda !=""]
  oo[,2]       <- V(gg)$GeneName#[gda !=""]
  
  ##--- loop over each disease
  for( d in 1:length(diseases) ){
    
    IDS <- V(gg)$name[grepl(diseases[d],gda,fixed=T)]
    N   <- length(IDS)
    
    ## for each gda, find the minimum shortest path to next gda (of the same disease)
    XX=igraph::shortest.paths(gg,IDS,IDS,weights=NA)
    diag(XX)       = NA
    ds             = apply(XX,1,min,na.rm=T)
    indX           = match(names(ds),oo[,1])
    oo[indX,(2+d)] = as.vector(ds)
    
    res[d,2]       = as.character(N)
    res[d,3]       = as.character(mean(ds))
    res[d,4]       = as.character(sd(ds))
    
  }
  
  DAB <- matrix(NA,ncol=length(diseases),nrow=length(diseases))
  colnames(DAB) <- diseases
  rownames(DAB) <- diseases
  
  #--- NOTE ---#
  # DAB is bound by -dmax <= DAB <= dmax
  # where dmax denotes the diameter of the network
  # dmax <- diameter(gg,directed=F)
  #------------#
  
  ##--- calculate disease-disease overlap
  for( i in 1:length(diseases) ){
    for( j in i:length(diseases) ){
      
      DAB[i,j] <- 0
      
      if( i != j ){
        DAB[i,j] <- diseaseOverlap(gg,gda,rownames(DAB)[i],colnames(DAB)[j],oo)
      }
      
    }
  }
  colnames(DAB) <- unescapeAnnotation(diseases)
  rownames(DAB) <- unescapeAnnotation(diseases)
  colnames(oo) <- c("Gene.ID","Gene.Name",unescapeAnnotation(diseases))
  
  if(permute=='none'){
    oo<-oo[gda !="",]
  }
  return(list(disease_separation=DAB,gene_disease_separation=oo,disease_localisation=res))
}

annotate_topOnto_ovg<-function(gg,dis){
  ids = V(gg)$name
  gg <- removeVertexTerm(gg,"TopOnto_OVG")
  gg <- removeVertexTerm(gg,"TopOnto_OVG_HDO_ID")
  
  #--- Set Disease (geneRIF db) attributes in .gml graph
  set.vertex.attribute(gg,"TopOnto_OVG",V(gg),"")
  set.vertex.attribute(gg,"TopOnto_OVG_HDO_ID",V(gg),"")
  
  #par    <- read.table("flatfile_human_gene2HDO.parentTerm.csv",sep="\t",skip=1,strip.white=T,quote="")
  #dis    <- read.table("flatfile_human_gene2HDO.csv",sep="\t",skip=1,header=F,strip.white=T,quote="")
  
  #dis    <- rbind(dis,par)
  
  disIDS <- dis[,3]
  
  disn <-getDiseases()
  dtype<- getDType()
  
  for( i in 1:length(ids) ){
    
    ind1 = which(disIDS==ids[i])
    
    Str1 <- "";
    Str2 <- "";
    
    if( length(ind1) != 0 ){
      #TDOD: refactor this code to work without disn
      disv <- as.vector(dis[ind1,1]);
      
      indx <- match(disv,disn)
      
      for( j in 1:length(disv) ){
        
        if( !is.na(indx[j]) ){
          
          if( Str1 == "" ) { Str1 <- as.character(dtype[indx[j]]) }
          else {
            Str1 <- paste(c(Str1,as.character(dtype[indx[j]])),collapse=COLLAPSE) }
          
          if( Str2 == "" ) { Str2 <- as.character(disn[indx[j]]) }
          else {
            Str2 <- paste(c(Str2,as.character(disn[indx[j]])),collapse=COLLAPSE) }
        }
        
      }
    }
    
    Str1 = paste(unique(strsplit(Str1,COLLAPSE)[[1]]),collapse=COLLAPSE)
    Str2 = paste(unique(strsplit(Str2,COLLAPSE)[[1]]),collapse=COLLAPSE)
    
    V(gg)[i]$TopOnto_OVG = as.character(Str1);
    V(gg)[i]$TopOnto_OVG_HDO_ID = as.character(Str2);
    
    
  }
  return(gg)
}

getAnnotationList<-function(annVec,col=COLLAPSE,sort=c('none','string','frequency')){
  sort <- match.arg(sort)
  res=switch (sort,
              none = unique(unlist(strsplit(annVec,';'))),
              string = sort(unique(unlist(strsplit(annVec,';')))),
              frequency = names(sort(table(unlist(strsplit(annVec,';'))),decreasing = TRUE))
  )
  return(res)
}


getDType<-function(){
  #---HDO Disease short names
  dtype  <- vector(length=12);
  dtype[1]   = "AD";
  dtype[2]   = "BD";
  dtype[3]   = "AUT";
  dtype[4]   = "SCH";
  dtype[5]   = "ASD";
  dtype[6]   = "Epi";
  dtype[7]   = "ID";
  dtype[8]   = "HTN";
  dtype[9]   = "HD";
  dtype[10]  = "PD";
  dtype[11]  = "FTD";
  dtype[12]  = "MS";
  #dtype[12]  = "DMH";
  #dtype[13]  = "CNSD";
  
  return(dtype)
}


getDiseases<-function(){
  #---HDO ID DISEASES of INTEREST
  disn    <- vector(length=12);
  disn[1]  <- "DOID:10652"#Alzheimer's_disease"
  disn[2]  <- "DOID:3312"#bipolar_disorder"
  disn[3]  <- "DOID:12849"#autistic_disorder"
  disn[4]  <- "DOID:5419"#schizophrenia"
  disn[5]  <- "DOID:0060041"#autism_spectrum_disorder
  disn[6]  <- "DOID:1826"#epilepsy_syndrome
  disn[7]  <- "DOID:1059"
  disn[8]  <- "DOID:10763"
  disn[9]  <- "DOID:12858"
  disn[10] <- "DOID:14330"
  disn[11] <- "DOID:9255"
  disn[12] <- "DOID:2377"
  return(disn)
}



prepareGDA<-function(gg,name){
  gda<-get.vertex.attribute(gg,name)
  gda<-escapeAnnotation(gda)
  return(gda)
}

ESC      <- "|"
COLLAPSE <- ";"


unescapeAnnotation<-function(annVec,col=COLLAPSE,esc=ESC){
  res<-gsub(esc,'',annVec,fixed = TRUE)
  return(res)
}

escapeAnnotation<-function(annVec,col=COLLAPSE,esc=ESC){
  if(any(grepl(esc,annVec,fixed = TRUE))){
    stop("Either already escaped or escape charecter found in annotation\n")
  }
  annList<-strsplit(annVec,col,fixed = TRUE)
  escFun<-function(.x){
    if(length(.x)>0){
      return(paste0(esc,.x,esc,collapse = ';'))
    }else{
      return("")
    }
  }
  
  res<-sapply(annList,escFun)
  return(res)
}

diseaseOverlap <- function(GG, GDA, disA, disB, OO){
  
  #disease A genes
  IDS1  <- V(GG)$name[grepl(disA,GDA,fixed=T)]
  NIDS1 <- length(IDS1)
  
  #disease B genes
  IDS2  <- V(GG)$name[grepl(disB,GDA,fixed=T)]
  NIDS2 <- length(IDS2)
  
  #disease A given B
  paths  <- igraph::shortest.paths(GG,IDS1,IDS2,weights=NA)
  dsA    <- as.numeric(as.vector(apply(paths,1,min)))
  
  #disease B given A
  paths  <- igraph::shortest.paths(GG,IDS2,IDS1,weights=NA)
  dsB    <- as.numeric(as.vector(apply(paths,1,min)))
  
  #network-based separation between disease A and B
  dAB <- (sum(dsA)+sum(dsB))/(NIDS1+NIDS2)
  
  #network-based localisation of disease A
  indA <- which(colnames(OO)==disA)
  dA   <- mean(as.numeric(as.vector(OO[OO[,indA[1]]!=".",indA[1]])))
  
  #network-based localisation of disease B
  indB <- which(colnames(OO)==disB)
  dB   <- mean(as.numeric(as.vector(OO[OO[,indB[1]]!=".",indB[1]])))
  
  #overlap between disease A and B
  sAB = as.numeric(dAB) - (as.numeric(dA)+as.numeric(dB))/2
  
  return(sAB)
  
}
permute <- function(GNS, N){
  
  temp <- sample(GNS,N,replace=F)
  
  return(temp)
  
}


zeroNA<-function(x){
  x[is.na(x)]<-0
  return(x)
}

stars    <- c("*","**","***")


runPermDisease<-function(gg,name,diseases=NULL,Nperm=100,alpha=c(0.05,0.01,0.001)){
  mean0<-function(x){
    return(mean(zeroNA(x)))
  }
  sd0<-function(x){
    return(sd(zeroNA(x)))
  }
  print("0")
  resD<-calcDiseasePairs(gg=gg,name=name,diseases=diseases,permute = 'none')
  ds<-resD$gene_disease_separation
  loc<-resD$disease_localisation
  resL<-lapply(1:Nperm,function(.x)calcDiseasePairs(gg=gg,name=name,diseases=diseases,permute='random'))
  resGDS<-sapply(resL,function(.x)apply(.x$gene_disease_separation[,3:dim(.x$gene_disease_separation)[2]],c(1,2),as.numeric),simplify = "array")
  m<-apply(resGDS,c(1,2),mean0)
  RANds<-cbind(as.data.frame(resL[[1]]$gene_disease_separation[,1:2]),as.data.frame(m))
  ##comment out for moment
  disn <- colnames(ds)[3:length(ds[1,])]
  print("1")
  ##--- output results file comparing observed disease pairs against randomised distribution.
  disease_location_sig           <- matrix(0 ,ncol=7, nrow=length(disn))
  colnames(disease_location_sig) <- c("HDO.ID","N","mean_ds","SD_ds","Ran_mean_ds","Ran_SD_ds","Utest.pvalue")
  disease_location_sig[,1]       <- disn
  disease_location_sig[,2]<-loc[match(disease_location_sig[,1],loc[,1]),2]
  
  ## significance of ds for each disease
  for( i in 1:length(disn) ){
    
    ## gda matching indices
    indx <- ds[,(2+i)]!="."
    
    ## gene ids
    ids <- ds[indx,1]
    
    ## observed ds values
    DS       <- as.numeric(as.vector(ds[indx,(2+i)]))
    disease_location_sig[i,3] <- as.numeric(mean(DS))
    disease_location_sig[i,4] <- as.numeric(sd(DS))
    
    indy <- match(ids,RANds[,1])
    
    ## random ds values
    RDS <- as.numeric(as.vector(RANds[indy,(2+i)]))
    disease_location_sig[i,5] <- as.numeric(mean0(RDS))
    disease_location_sig[i,6] <- as.numeric(sd0(RDS))
    disease_location_sig[i,7] <- 1.0
    
    ## compute wilcox test between observable ds and random ds, and store p.values,
    ## see (Menche et al., 2015).
    if( !is.infinite(DS) && !is.nan(DS) && !is.na(DS) &&  !is.infinite(RDS) && !is.nan(RDS) && !is.na(RDS) ){
      if( length(DS) != 0 && length(RDS) != 0 ){
        wt       <- wilcox.test(DS,RDS)
        disease_location_sig[i,7] <- as.numeric(wt$p.value)
      }
    }
  }
  print("2")
  sAB<-resD$disease_separation
  RAW_sAB<-sapply(resL,function(.x).x$disease_separation,simplify = "array")
  RAN_sAB_mean<-apply(RAW_sAB,c(1,2),mean0)
  RAN_sAB_sd<-apply(RAW_sAB,c(1,2),sd0)
  perms <- dim(RAW_sAB)[3]
  Nn    <- length(disn)
  NELE  <- Nn*(Nn+1)/2
  
  ##---no: of levels for Bonferroni correction
  Nlevels = NELE;
  print("3")
  ##--- Output file for disease-disease separation/overlap
  #CN <-  c("HDO.ID","Disease.long","Disease","N","HDO.ID","Disease.long","Disease","N","sAB","Separated","Overlap","zScore","pvalue","Separation/Overlap.than.chance","Bonferroni","p.adjusted","q-value")
  CN <-  c("HDO.ID","N","HDO.ID","N","sAB","Separated","Overlap","zScore","pvalue","Separation/Overlap.than.chance","Bonferroni","p.adjusted","q-value")
  zs <- matrix(".", nrow=NELE, ncol=length(CN))
  colnames(zs) <- CN
  tests <- matrix(0, nrow=NELE,ncol=perms)
  for( k in 0:(NELE-1) ){
    
    ##--- linear indexing for symmetric matrix
    i = floor( (2*Nn+1 - sqrt( (2*Nn+1)*(2*Nn+1) - 8*k ))/2 );
    j = k - Nn*i + i*(i-1)/2;
    
    i = i + 1;
    j = j + i;
    
    zScore = 0
    
    if( !is.nan(as.numeric(RAN_sAB_sd[i,j])) ){
      
      ## compute z-score, i.e. separation of mean sAB, against a randomised model (of the mean of sAB),
      ## see (Menche et al., 2015).
      if( as.numeric(RAN_sAB_sd[i,j]) != 0){
        zScore = (as.numeric(as.vector(sAB[i,j])) - as.numeric(as.vector(RAN_sAB_mean[i,j])))/(as.numeric(as.vector(RAN_sAB_sd[i,j])))
      }
      
      ## compute p.value from the normal distribution
      ## See also http://www.cyclismo.org/tutorial/R/pValues.html
      pval <- pnorm(-abs(zScore))
      pval <- 2 * pval
      
      zs[(k+1),1] <- disn[i]
      zs[(k+1),2] <- as.character(loc[which(loc[,1]==disn[i]),2])
      
      zs[(k+1),3] <- disn[j]
      zs[(k+1),4] <- as.character(loc[which(loc[,1]==disn[j]),2])
      
      ## sAB, the disease-disease separation/overlap measure, on the interactome
      zs[(k+1),5] <- as.character(sAB[i,j])
      
      ## sAB > 0, implies separation
      zs[(k+1),6] <- ifelse((as.numeric(zs[(k+1),5]) > 0), "YES", ".")
      
      ## sAB < 0, implies overlap
      zs[(k+1),7] <- ifelse((as.numeric(zs[(k+1),5]) < 0), "YES", ".")
      
      ## save z-score and p.value
      zs[(k+1),8] <- as.character(zScore)
      zs[(k+1),9] <- as.character(pval)
      
      ## z-scores < 0 (>0), implies separation/overlap smaller (larger) than by chance
      zs[(k+1),10] <- ifelse((as.numeric(zs[(k+1),8]) < 0), "Smaller", "larger")
      
      ## Bonferroni correction for p.value ('stars' can be found in 'setUp.R')
      temp <- "."
      for( x in 1:length(alpha) ){
        if(as.numeric(zs[(k+1),9]) < as.numeric(alpha[x]/Nlevels)){ temp <- stars[x] }
      }
      
      ## save the Bonerroni correction, repersented by stars, here.
      zs[(k+1),11] <- temp
      ## default fill of output container
    } else {
      zs[(k+1),1] <- disn[i]
      zs[(k+1),2] <- as.character(loc[which(loc[,1]==disn[i]),2])
      
      zs[(k+1),3] <- disn[j]
      zs[(k+1),4] <- as.character(loc[which(loc[,1]==disn[j]),2])
    }
    # if( zs[(k+1),1] != zs[(k+1),3] ){
    #
    #   test <- 0
    #
    #   Mn   <- as.numeric(RAN_sAB_mean[i,j])
    #   Sd   <- as.numeric(RAN_sAB_sd[i,j])
    #
    #   RsAB <- as.numeric(RAW_sAB[i,j,])
    #
    #   ## store the random z-score
    #   tests[(k+1),] <- (as.numeric(RsAB - Mn))/Sd
    #
    # }
    
  }
  print("4")
  ## save p.adjusted value in output container
  zs[,12] <- p.adjust(as.numeric(zs[,9]),method="BY")
  
  ## if calFDR is FALSE, we'll use WGCNA's qvalue calculation for FDR
  zs[,13] <- qvalue(as.numeric(zs[,9]))$qvalue
  
  return(list(Disease_overlap_sig=zs,Disease_location_sig=disease_location_sig))
  
}

layoutByCluster<-function(gg,mem,layout=layout_with_kk){
  Cn<-table(mem$membership)
  sgraphs<-lapply(names(Cn),getClusterSubgraphByID,gg=gg,mem=mem$membership)
  layouts <- lapply(sgraphs, layout)
  lay <- merge_coords(sgraphs, layouts)
  ug <- disjoint_union(sgraphs)
  idx<-match(V(gg)$name,V(ug)$name)
  lay<-lay[idx,]
  return(lay)
}


getClusterSubgraphByID<-function(clID,gg,mem){
  idx<-which(mem==clID)
  sg<-induced_subgraph(gg,V(gg)[idx],impl = "auto")
  return(sg)
}

getCommunityGraph<-function(gg,membership){
  g<-gg
  V(g)$composition<-V(gg)$name
  cgg<-simplify(contract(g, membership,vertex.attr.comb =list(composition='concat', 'ignore')))
  V(cgg)$name<-as.character(V(cgg))
  V(cgg)$size<-sapply(V(cgg)$composition,length)
  return(cgg)
}







# This function is taken taken from the WGCNA package available at https://github.com/cran/WGCNA. It is used in the function runPermDisease.



qvalue <- function(p, lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=NULL, robust=FALSE, 
                   smooth.df = 3, smooth.log.pi0 = FALSE) {
  #Input
  #=============================================================================
  #p: a vector of p-values (only necessary input)
  #fdr.level: a level at which to control the FDR (optional)
  #lambda: the value of the tuning parameter to estimate pi0 (optional)
  #pi0.method: either "smoother" or "bootstrap"; the method for automatically
  #           choosing tuning parameter in the estimation of pi0, the proportion
  #           of true null hypotheses
  #robust: an indicator of whether it is desired to make the estimate more robust
  #        for small p-values and a direct finite sample estimate of pFDR (optional)
  #gui: A flag to indicate to 'qvalue' that it should communicate with the gui.  ## change by Alan
  #     Should not be specified on command line.
  #smooth.df: degrees of freedom to use in smoother (optional)
  #smooth.log.pi0: should smoothing be done on log scale? (optional)
  #
  #Output
  #=============================================================================
  #call: gives the function call
  #pi0: an estimate of the proportion of null p-values
  #qvalues: a vector of the estimated q-values (the main quantity of interest)
  #pvalues: a vector of the original p-values
  #significant: if fdr.level is specified, an indicator of whether the q-value
  #    fell below fdr.level (taking all such q-values to be significant controls
  #    FDR at level fdr.level)
  
  #This is just some pre-processing
  
  if(min(p)<0 || max(p)>1) 
    stop("qvalue: p-values not in valid range.")
  if(length(lambda)>1 && length(lambda)<4) 
    stop("qvalue: If length of lambda greater than 1, you need at least 4 values.")
  if(length(lambda)>1 && (min(lambda) < 0 || max(lambda) >= 1)) 
    stop("qvalue: Lambda must be within [0, 1).")
  m <- length(p)
  #These next few functions are the various ways to estimate pi0
  if(length(lambda)==1) {
    if(lambda<0 || lambda>=1) 
      stop("qvalue: Lambda must be within [0, 1).")
    
    pi0 <- mean(p >= lambda)/(1-lambda)
    pi0 <- min(pi0,1)
  }
  else {
    pi0 <- rep(0,length(lambda))
    for(i in 1:length(lambda)) {
      pi0[i] <- mean(p >= lambda[i])/(1-lambda[i])
    }
    
    if(pi0.method=="smoother") {
      if(smooth.log.pi0)
        pi0 <- log(pi0)
      
      spi0 <- smooth.spline(lambda,pi0,df=smooth.df)
      pi0 <- predict(spi0,x=max(lambda))$y
      
      if(smooth.log.pi0)
        pi0 <- exp(pi0)
      pi0 <- min(pi0,1)
    }
    else if(pi0.method=="bootstrap") {
      minpi0 <- min(pi0)
      mse <- rep(0,length(lambda))
      pi0.boot <- rep(0,length(lambda))
      for(i in 1:100) {
        p.boot <- sample(p,size=m,replace=TRUE)
        for(i in 1:length(lambda)) {
          pi0.boot[i] <- mean(p.boot>lambda[i])/(1-lambda[i])
        }
        mse <- mse + (pi0.boot-minpi0)^2
      }
      pi0 <- min(pi0[mse==min(mse)])
      pi0 <- min(pi0,1)
    }
    else {  ## change by Alan: check for valid choice of 'pi0.method' (only necessary on command line)
      stop("qvalue:: 'pi0.method' must be one of 'smoother' or 'bootstrap'.")
      return(0)
    }
  }
  if(pi0 <= 0) 
    stop("qvalue:: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.")
  if(!is.null(fdr.level) && (fdr.level<=0 || fdr.level>1))  ## change by Alan:  check for valid fdr.level
    stop("qvalue:: 'fdr.level' must be within (0, 1].")
  #The estimated q-values calculated here
  u <- order(p)
  
  # change by Alan
  # ranking function which returns number of observations less than or equal
  qvalue.rank <- function(x) {
    idx <- sort.list(x)
    
    fc <- factor(x)
    nl <- length(levels(fc))
    bin <- as.integer(fc)
    tbl <- tabulate(bin)
    cs <- cumsum(tbl)
    
    tbl <- rep(cs, tbl)
    tbl[idx] <- tbl
    
    return(tbl)
  }
  
  v <- qvalue.rank(p)
  
  qvalue <- pi0*m*p/v
  if(robust) {
    qvalue <- pi0*m*p/(v*(1-(1-p)^m))
  }
  qvalue[u[m]] <- min(qvalue[u[m]],1)
  for(i in (m-1):1) {
    qvalue[u[i]] <- min(qvalue[u[i]],qvalue[u[i+1]],1)
  }
  #The results are returned
  if(!is.null(fdr.level)) {
    retval <- list(call=match.call(), pi0=pi0, qvalues=qvalue, pvalues=p, fdr.level=fdr.level, ## change by Alan
                   significant=(qvalue <= fdr.level), lambda=lambda)
  }
  else {
    retval <- list(call=match.call(), pi0=pi0, qvalues=qvalue, pvalues=p, lambda=lambda)
  }
  class(retval) <- "qvalue"
  return(retval)
}
