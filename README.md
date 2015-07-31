#PajekTo3DStereo

PajekTo3DStereo is a file converter with a Tcl/Tk GUI built in R. This project was published and presented as a demonstration at the American Medical Informatics Association Annual Symposium 2014.

Dang, B.S., Bhavnani S.K. Pajekto3DStereo: Enabling Generation and Interaction with 3D Stereo Networks. Proceedings of AMIA (2014).
[PDF Link](http://www.skbhavnani.com/DIVA/papers/Dang-Bhavnani-3D-Demo-AMIA-2014.pdf)

Code Author: Bryant Dang


##What can it be used for?

The purpose of this application is to convert Pajek network files into VMD visualization files, in order to utilize the 3D stereo features in VMD.


##How to use

The application is launched from the main.R script. 

The main application consists of two GUI windows. The "File Selection Window" is used to select the Pajek data files that contain the network to be converted. Once the files have been loaded, pushing the *Read Files* button will generate a second window, the "Visualization Menu", which allows the user to select different visualization parameters to be added to the VMD output file. Once selections have been made, pushing the *Save and Convert* button will promt the user for the location to save the output files, and then generate and save the output VMD files.

Note: This application uses the package gWidgets with TCL/TK to create the main GUI windows. If this application is launched using RStudio, the calls to gWidgets will produce warnings, which may be ignored.