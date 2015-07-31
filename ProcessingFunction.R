readfiles <- function()
{ 
  if(enable==0)
  {
    winDialog("ok", "A valid Network File must be selected!")
    stop()
  }
  
  if(exists('partitionDialogWindow')==TRUE)
  {
    dispose(partitionDialogWindow)
  }  
  
  #Reads .net file to determine number of nodes in Set1 and Set2 if network is Bipartite
  numAtoms <- read.table(file.path(nPath,svalue(netFile)), nrows = 1)
  NodesinSet1 <- ifelse(is.null(numAtoms$V3)==TRUE, numAtoms$V2, numAtoms$V3)
  NodesinSet2 <- ifelse(is.null(numAtoms$V3)==TRUE,0, numAtoms$V2 - numAtoms$V3)
  TotalNodes <- numAtoms$V2
  
  #Reads .net, .clu, and .vec files to fill respective variables
  Atoms <- read.table(file.path(nPath,svalue(netFile)), skip = 1, nrows = TotalNodes, col.names = c("AtomNum", "AtomLabel", "Xpos", "Ypos", "Zpos"))
  Edges <- read.table(file.path(nPath,svalue(netFile)), skip = TotalNodes+2, col.names = c("Point1", "Point2", "EdgeWeight"))
  
  if(svalue(cluFile, index=TRUE)>1)
  {
    Clusters <- read.table(file.path(cPath,svalue(cluFile)), col.names = "ClusterNum", skip=1)
  } else {
    Clusters <- data.frame(ClusterNum=rep(1,nrow(Atoms)))
  }
  
  if(svalue(vecFile, index=TRUE)>1)
  {
    NodeWeights <- read.table(file.path(vPath,svalue(vecFile)), col.names = "NodeWeights", skip=1)
  } else {
    NodeWeights <- data.frame(NodeWeights=rep(0.5,nrow(Atoms)))
  }
  
  if(nrow(Atoms) != nrow(Clusters))
  {
    winDialog("ok", "Data mismatch!\n\nCluster file does not match Network file!")
    stop()
  }
  
  if(nrow(Atoms) != nrow(NodeWeights))
  {
    winDialog("ok", "Data mismatch!\n\nVector file does not match Network file!")
    stop()
  }
  
  #Numbers Set1 and Set2 if Bipartite
  nodenumbering <- data.frame(V1=1:NodesinSet1)
  if(is.null(numAtoms$V3)==FALSE) {
    nodenumbering2 <- data.frame(V1=1:NodesinSet2)
    nodenumbering <- rbind(nodenumbering,nodenumbering2)}
  
  #Normalizes EdgeWeights, Applies function(5x) to entire column
  trans <- function(x) x/8
  Edges$EdgeWeight <- format(lapply(Edges$EdgeWeight, trans), digits=4)
  
  #Builds color transformation table
  uClusters <- unique(Clusters)
  VMDcolors <- c("black", "green", "red", "yellow", "blue", "orange", "purple", "cyan","gray", "pink", "magenta", "white")
  ColorTransform <- data.frame("Colors" = c(16,7,1,4,0,3,11,10,2,9,27, 8), row.names=VMDcolors)
  cb1 <- vector("list", 50)
  cb2 <- vector("list", 50)
  
  #Determines unique clusters and builds GUI. Queries user for color and label choices
  winAdjust1 <- 0
  winAdjust2 <- 0  
  check.integer <- function(N){
    !length(grep("[^[:digit:]]", as.character(N)))
  }  
  partitionDialogWindow <<- gwindow("Pajekto3DStereo")
  tbl <- glayout(cont=partitionDialogWindow, horizontal=TRUE)
  tbl[1,1] <- "Cluster:" 
  tbl[1,2] <- "   Select Color:"
  tbl[1,3] <- "Show Label?:"  
  for(t in 1:nrow(unique(Clusters)))
  {
    if (check.integer((t-1)/10)==TRUE && t > 10)
    {
      winAdjust1 <- winAdjust1 +3
      winAdjust2 <- winAdjust2 +10
    }
    tbl[t+1-winAdjust2,winAdjust1+1] <- glabel(sub("","      ",uClusters$ClusterNum[t]), cont=tbl)
    tbl[t+1-winAdjust2,winAdjust1+2] <- cb1[uClusters$ClusterNum[t]] <- gcombobox(VMDcolors, cont=tbl, selected=1)
    tbl[t+1-winAdjust2,winAdjust1+3] <- cb2[uClusters$ClusterNum[t]] <-gcombobox(c("No","Yes"), cont=tbl, selected=1)
  }
  tbl[t+2-winAdjust2,2+winAdjust1] <- "  Edges to Write:"
  tbl[t+2-winAdjust2,3+winAdjust1] <- thresholdEdges <- gcombobox(c(nrow(Edges),round(nrow(Edges)*.9),round(nrow(Edges)*.8),round(nrow(Edges)*.7),round(nrow(Edges)*.6),round(nrow(Edges)*.5),round(nrow(Edges)*.4),round(nrow(Edges)*.3),round(nrow(Edges)*.2),round(nrow(Edges)*.1)), cont=tbl, selected=1)
  tbl[t+3-winAdjust2,2+winAdjust1] <- "    Node Scaling:"
  tbl[t+3-winAdjust2,3+winAdjust1] <- scaleNodes <- gcombobox(c(1,1.5,2,2.5,3), cont=tbl, selected=1)
  tbl[t+4-winAdjust2,2:3+winAdjust1] <- helpButton2 <- gbutton("             Info              ", cont=tbl, handler=infohandler2)  
  tbl[t+5-winAdjust2,2:3+winAdjust1] <- calcButton <- gbutton("  Save and Convert  ", cont=tbl, handler = convert) 
}