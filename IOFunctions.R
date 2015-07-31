# Functions for loading Pajek network files
netHandler <- function(netFile) {
  nPath <<- dlgDir(title = "Select Directory containing Network File", dir = nPath)$res
  if (length(dir(path=nPath, pattern=".net")) < 1) {
    winDialog("ok", "Selected Directory has no Network Files!")
    nPath <<- snPath
  } else {
    netFile[] <<- dir(path=nPath, pattern=".net")
    snPath <<- nPath
    enable <<- 1
  }
}

cluHandler <- function(cluFile) {
  cPath <<- dlgDir(title = "Select Directory containing Network File", dir = nPath)$res
  if (length(dir(path=cPath, pattern=".clu")) < 1) {
    winDialog("ok", "Selected Directory has no Cluster Files!")
    cPath <<- scPath
  } else {
    cluFile[] <<- c("No Cluster File", dir(path=cPath, pattern=".clu"))
    scPath <<- cPath
  }
}

vecHandler <- function(vecFile) {
  vPath <<- dlgDir(title = "Select Directory containing Vector File", dir = nPath)$res
  if (length(dir(path=vPath, pattern=".vec")) < 1) {
    winDialog("ok", "Selected Directory has no Vector Files!")
    vPath <<- svPath
  } else {
    vecFile[] <<- c("No Vector File", dir(path=vPath, pattern=".vec"))
    svPath <<- vPath
  }
}

infohandler1 <- function(h,...) {
  winDialog("ok", "Not Available")
}

infohandler2 <- function(h,...) {
  winDialog("ok", "Not Available")
}