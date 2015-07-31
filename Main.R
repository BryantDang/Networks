## Author: Bryant Dang
## 
## The purpose of this application is to convert Pajek network files into VMD visualization files, in order to utilize the 3D stereo features in VMD
##
## Required Pacakges:
##   gdata, svDialogs, gWidgets, gWidgetstcltk
##
## Note: This application uses the package gWidgets with TCL/TK to create the main GUI windows. 
##       If using RStudio, the calls to gWidgets will produce warnings, which may be ignored.

#Clear variables in the workspace 
rm(list = ls())

#Loads libraries
library(gdata)
library(svDialogs)
library(gWidgets)
library(gWidgetstcltk)
options(guiToolkit='tcltk')

#Source Function Files
source("IOFunctions.R")
source("ProcessingFunction.R")

nPath <- getwd()
cPath <- 0
vPath <- 0
enable <- 0


FileWindow <- gwindow("Pajekto3DStereo - File Selection Menu")
fileTbl <- glayout(cont=FileWindow, horizontal=TRUE)
fileTbl[1,1] <- "Network File (Req):"
fileTbl[1,2] <- netFile <- gcombobox("Use 'Select Directory' button to locate Network Files", cont=fileTbl, selected=1)
fileTbl[1,3] <- netButton <- gbutton("Select Directory", cont=fileTbl, handler = netHandler)
fileTbl[2,1] <- "  Cluster File (Opt):"
fileTbl[2,2] <- cluFile <- gcombobox("Use 'Select Directory' button to locate Cluster Files", cont=fileTbl, selected=1)
fileTbl[2,3] <- cluButton <- gbutton("Select Directory", cont=fileTbl, handler = cluHandler)
fileTbl[3,1] <- "   Vector File (Opt):"
fileTbl[3,2] <- vecFile <- gcombobox("Use 'Select Directory' button to locate Vector Files ", cont=fileTbl, selected=1)
fileTbl[3,3] <- vecButton <- gbutton("Select Directory", cont=fileTbl, handler = vecHandler)
fileTbl[4,3] <- fileButton <- gbutton("     Read Files     ", cont=fileTbl, handler = readfiles)
fileTbl[4,2] <- helpButton1 <- gbutton("             Info             ", cont=fileTbl, handler=infohandler1)

winDialog("ok", "Welcome to the 3D Pajek to VMD Converter!\n\nUse the 'Select Directory' button on the top right of the File Selection Menu to choose the folder the Network File (.net) is in.")
