#Creates conversion function from Pajek to VMD
convert <- function(){
  
  pdbPath <<- guiDlgSave(title = "Save Pdb file as...", 
                         defaultDir = snPath, 
                         defaultFile = sub(".net$", "", svalue(netFile)),
                         defaultExt = ".pdb",
                         filters = c("Pdb Files", ".pdb"))
  
  if(nchar(pdbPath)<1)
  {
    winDialog("ok", "Save file was not selected!")
    stop()
  }
  
  vmdPath <<- guiDlgSave(title = "Save VMD file as...", 
                         defaultDir = snPath, 
                         defaultFile = sub(".net$", "", svalue(netFile)),
                         defaultExt = ".vmd",
                         filters = c("VMD Files", ".vmd"))
  if(nchar(vmdPath)<1)
  {
    winDialog("ok", "Save file was not selected!")
    stop()
  }
  
  pb <- winProgressBar("Conversion in progress...", "Writing PDB Table",
                       0, 100, 0)
  Sys.sleep(0.5)
  #Builds data.frame for pdb output file
  pdbtable <- data.frame(V1="ATOM",Atoms$AtomNum)
  #pdbtable$V3 <- with(pdbtable, ifelse(Atoms$AtomNum<=NodesinSet1, "P", "C"))
  pdbtable$V3 <- formatC(nodenumbering$V1, width=4, format = "d", flag="0")
  pdbtable$V4 <- with(pdbtable, ifelse(Atoms$AtomNum<=NodesinSet1, "P", "C"))
  pdbtable$V5 <- formatC(Clusters$ClusterNum, width=2, format = "d", flag="0")
  pdbtable$V6 <- with(pdbtable, ifelse(Atoms$AtomNum<=NodesinSet1, "P", "C"))
  pdbtable$V7 <- pdbtable$Atoms.AtomNum
  pdbtable$V8 <- sprintf("%.3f", Atoms$Xpos*100)
  pdbtable$V9 <- sprintf("%.3f",Atoms$Ypos*100)
  pdbtable$V10 <- sprintf("%.3f",Atoms$Zpos*100)
  
  #Wrties output.pdb in proper fixed-width formated table
  write.fwf(pdbtable,file=pdbPath,width=c(4,7,5,2,2,2,4,12,8,8),justify="right",quote=F,colnames=F,rownames=F, sep="")
  cat("END", file=pdbPath, append=TRUE)
  
  setWinProgressBar(pb, 5, label = "Creating VMD File")
  Sys.sleep(0.5)
  #Writes intro to vmd script
  fileConn<-file(vmdPath)
  writeLines(c("set viewplist {}
               set fixedlist {}
               proc vmdrestoremymaterials {} {
               set mlist { Opaque Transparent BrushedMetal Diffuse Ghost Glass1 Glass2 Glass3 Glossy HardPlastic MetallicPastel Steel Translucent Edgy EdgyShiny EdgyGlass Goodsell AOShiny AOChalky AOEdgy }
               set mymlist [material list]
               foreach mat $mlist {
               if { [lsearch $mymlist $mat] == -1 } { 
               material add $mat
               }
               }         
               material change ambient AOChalky 0.000000
               material change diffuse AOChalky 0.850000
               material change specular AOChalky 0.000000
               material change shininess AOChalky 0.530000
               material change opacity AOChalky 1.000000
               material change outline AOChalky 0.000000
               material change outlinewidth AOChalky 0.000000
               }
               vmdrestoremymaterials
               atomselect macro unparametrized beta<1
               atomselect macro addedmolefacture {occupancy 0.8}
               # Display settings
               display eyesep       0.250000
               display focallength  3.620000
               display height       5.100000
               display distance     96.000000
               display projection   Orthographic
               display nearclip set 0.010000
               display farclip  set 8.950000
               display depthcue   off
               display cuestart   0.500000
               display cueend     10.000000
               display cuedensity 0.500000
               display cuemode    Exp2
               graphics top color 2
               material on
               material Transparent"), fileConn)
  close(fileConn)
  
  setWinProgressBar(pb, 10, label = "Writing Edges to VMD File")
  
  #Appends output.vmd file with a formated draw command for each edge
  edgeCount <- 1
  #  for(i in nrow(Edges):1)
  #  {
  #    if(edgeCount<=as.numeric(svalue(thresholdEdges)))
  #    {
  #      Co1 <- Edges[i,1]
  #      Co2 <- Edges[i,2]
  #      #if statement represents lower bound cutoff for edge weights, max = 6. (higher number = less lines)  
  #      #cat("graphics top cylinder {", pdbtable[Co1,8],"000 ",pdbtable[Co1,9],"000 ",pdbtable[Co1,10],"000} {",pdbtable[Co2,8],"000 ",pdbtable[Co2,9],"000 ",pdbtable[Co2,10],"000} resolution 4 filled no radius ", ifelse(as.numeric(Edges$EdgeWeight[i])>=0.001, sprintf("%.3f", as.numeric(Edges$EdgeWeight[i])), "0.001"), "\n", sep="", file=vmdPath, append = TRUE)
  #      cat("graphics top line {", pdbtable[Co1,8],"000 ",pdbtable[Co1,9],"000 ",pdbtable[Co1,10],"000} {",pdbtable[Co2,8],"000 ",pdbtable[Co2,9],"000 ",pdbtable[Co2,10],"000} style solid width 1", "\n", sep="", file=vmdPath, append = TRUE)
  #      edgeCount<-edgeCount+1
  #    }
  #  }
  
  if(svalue(vecFile, index=TRUE)>1)
  {
    #Min-Max Normalization on NodeWeights, and then multiplied for resizing effect in VMD
    NodeMin <- min(NodeWeights)
    NodeMax <- max(NodeWeights)
    NodeMMN <- function(x) ((x-NodeMin)/(NodeMax-NodeMin) * as.numeric(svalue(scaleNodes)) + 0.1)
    if(svalue(vecFile, index=TRUE)>1) 
    {
      NodeWeights$NodeWeights <- lapply(NodeWeights$NodeWeights, NodeMMN)
    }
  }
  
  setWinProgressBar(pb, 50, label = "Writing Nodes to VMD File")
  Sys.sleep(0.5)
  #Writes row-stacked objects for each atom. First 2 lines represent nodesize and color respectively.
  
  apply(pdbtable[], 1, function(x) 
    cat("draw color ", ColorTransform[svalue(cb1[[Clusters$ClusterNum[as.numeric(x[7])]]]),], "\n",
        "draw sphere {", x[8], " ", x[9], " ", x[10], "} radius ", format(NodeWeights$NodeWeights[as.numeric(x[7])], digits=3), "\n",
        sep="", file=vmdPath, append=TRUE ))
  
  #for(j in 1:TotalNodes)
  #{
  #  cat("draw sphere {"
  #    "mol representation VDW ", format(NodeWeights$NodeWeights[j], digits=3), " 26.0 \n",
  #      "mol color ColorID ", ColorTransform[svalue(cb1[[Clusters$ClusterNum[j]]]),], "\n",
  #      "mol selection {resid ", j, "} \n", 
  #      "mol material AOChalky \n", 
  #      "mol addrep top \n",
  #      "mol selupdate ", j, " top 0 \n",
  #      "mol colupdate ", j, " top 0 \n",
  #      "mol scaleminmax top ", j, " 0.000000 0.000000 \n",
  #      "mol smoothrep top ", j, " 0 \n",
  #      "mol drawframes top ", j, " {now} \n",
  #      sep="", file=vmdPath, append=TRUE)
  #}
  
  setWinProgressBar(pb, 85, label = "Writing Options to VMD File")
  Sys.sleep(0.5)
  #Writes formatting code for stage presentation.
  sink(vmdPath, append=TRUE)
  cat(c("set viewpoints([molinfo top]) {{{1 0 0 -56.502} {0 1 0 -62.4307} {0 0 1 -58.0376} {0 0 0 1}} {{0.328035 -0.755738 0.56772 0} {-0.46883 0.392011 0.792307 0} {-0.821266 -0.526269 -0.225861 0} {0 0 0 1}} {{0.116153 0 0 0} {0 0.116153 0 0} {0 0 0.116153 0} {0 0 0 1}} {{1 0 0 -0.189999} {0 1 0 0.06} {0 0 1 -3.75} {0 0 0 1}}}
        lappend viewplist [molinfo top]
        set topmol [molinfo top]
        done with molecule 0
        foreach v $viewplist {
        molinfo $v set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($v)
        }
        foreach v $fixedlist {
        molinfo $v set fixed 1
        }
        unset viewplist
        unset fixedlist
        mol top $topmol
        unset topmol
        proc vmdrestoremycolors {} {
        color scale colors RWB {1.0 0.0 0.0} {1.0 1.0 1.0} {0.0 0.0 1.0}
        color scale colors BWR {0.0 0.0 1.0} {1.0 1.0 1.0} {1.0 0.0 0.0}
        color scale colors RGryB {1.0 0.0 0.0} {0.5 0.5 0.5} {0.0 0.0 1.0}
        color scale colors BGryR {0.0 0.0 1.0} {0.5 0.5 0.5} {1.0 0.0 0.0}
        color scale colors RGB {1.0 0.0 0.0} {0.0 1.0 0.0} {0.0 0.0 1.0}
        color scale colors BGR {0.0 0.0 1.0} {0.0 1.0 0.0} {1.0 0.0 0.0}
        color scale colors RWG {1.0 0.0 0.0} {1.0 1.0 1.0} {0.0 1.0 0.0}
        color scale colors GWR {0.0 1.0 0.0} {1.0 1.0 1.0} {1.0 0.0 0.0}
        color scale colors GWB {0.0 1.0 0.0} {1.0 1.0 1.0} {0.0 0.0 1.0}
        color scale colors BWG {0.0 0.0 1.0} {1.0 1.0 1.0} {0.0 1.0 0.0}
        color scale colors BlkW {0.0 0.0 0.0} {0.5 0.5 0.5} {1.0 1.0 1.0}
        color scale colors WBlk {1.0 1.0 1.0} {0.5 0.5 0.5} {0.0 0.0 0.0}
        color scale method GWB
        set colorcmds {
{color Display {Background} white}
{color Display {BackgroundTop} black}
{color Display {BackgroundBot} blue2}
{color Display {FPS} white}
{color Element {X} cyan}
{color Chain {P} blue}
{color Chain {C} red}
{color Segname {} blue}
{color Conformation {all} blue}
{color Molecule {0} blue}
{color Structure {3_10_Helix} blue}
{color Surface {Grasp} gray}
{color Labels {Atoms} black}
{color Labels {Springs} orange}
{color Stage {Even} gray}
{color Stage {Odd} silver}
        }
        foreach colcmd $colorcmds {
        set val [catch {eval $colcmd}]
        }
        color change rgb 0 0.0 0.0 1.0
        color change rgb 2 0.9700000286102295 0.9700000286102295 0.9700000286102295
        color change rgb 3 1.0 0.5299999713897705 0.0
        color change rgb 4 1.0 1.0 0.0
        color change rgb 5 0.5 0.5 0.20000000298023224
        color change rgb 6 0.6000000238418579 0.6000000238418579 0.6000000238418579
        color change rgb 7 0.0 1.0 0.0
        color change rgb 8 0.800000011920929 0.800000011920929 0.800000011920929
        color change rgb 9 1.0 0.6000000238418579 0.6000000238418579
        color change rgb 11 0.6499999761581421 0.0 0.6499999761581421
        color change rgb 12 0.5 0.8999999761581421 0.4000000059604645
        color change rgb 13 0.8999999761581421 0.4000000059604645 0.699999988079071
        color change rgb 14 0.5 0.30000001192092896 0.0
        color change rgb 15 0.5 0.5 0.5
        color change rgb 17 0.8799999952316284 0.9700000286102295 0.019999999552965164
        color change rgb 18 0.550000011920929 0.8999999761581421 0.019999999552965164
        color change rgb 19 0.0 0.8999999761581421 0.03999999910593033
        color change rgb 20 0.0 0.8999999761581421 0.5
        color change rgb 21 0.0 0.8799999952316284 1.0
        color change rgb 22 0.0 0.7599999904632568 1.0
        color change rgb 23 0.019999999552965164 0.3799999952316284 0.6700000166893005
        color change rgb 24 0.009999999776482582 0.03999999910593033 0.9300000071525574
        color change rgb 25 0.27000001072883606 0.0 0.9800000190734863
        color change rgb 26 0.44999998807907104 0.0 0.8999999761581421
        color change rgb 27 0.8999999761581421 0.0 0.8999999761581421
        color change rgb 28 1.0 0.0 0.6600000262260437
        color change rgb 29 0.9800000190734863 0.0 0.23000000417232513
        color change rgb 30 0.8100000023841858 0.0 0.0
        color change rgb 31 0.8899999856948853 0.3499999940395355 0.0
        color change rgb 32 0.33000001311302185 0.33000001311302185 0.33000001311302185
        }
        vmdrestoremycolors
        axes location off
        label textsize 0.9142857193946838
        display resetview \n"
    ))
    sink()
    
    setWinProgressBar(pb, 90, label = "Writing Labels to VMD File")
    Sys.sleep(0.5)
    nLabels<-0 
    #Draws labels for nodes. Labels drawn are determined by user input through GUI
    for(k in 1:TotalNodes)
    { 
      if(svalue(cb2[[Clusters$ClusterNum[k]]])=="Yes"){
        l<-as.numeric(NodeWeights$NodeWeights[k])/30
        cat(
          "label add Atoms ", "0/",k-1,"\n",
          "label textoffset Atoms ", nLabels, " {",l," ", l,"}\n",
          "label textformat Atoms ", nLabels, " ", toString(Atoms$AtomLabel[k]), "\n",
          file=vmdPath, append=TRUE, sep="")
        nLabels<-nLabels+1}
    }
    setWinProgressBar(pb, 100, label = "Closing VMD File")
    Sys.sleep(1)
    close(pb)
    winDialog("ok","Conversion Completed!")
  }