#=======================================#
### INSTALLATION INSTRUCTIONS FOR MAC ###
#=======================================#

#######Begin by installing xcode and xcode command line tools#######

##Install Xcode using the appstore (free)

##Open a terminal and instal Xcode command line tools
xcode-select --install

##Agree to Xcode license in Terminal
sudo xcodebuild -license



#######The following steps may be necessary in order to install some of chromstaR dependencies#######

##In a terminal window, check if libcurl, libxml and ccache are installed
curl-config --version
xmllint --version
ccache --version

##if not installed, install (or update) them through macports
##Install macports from this website https://www.macports.org/install.php
##Go back to terminal and selfupdate macports
sudo port -v selfupdate
##install libcurl, libxml and ccache
sudo port install curl libxml libxml2 ccache


##Install the gfortran libraries
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /



#######As of June 5th, 2015, the following steps are needed to install chromstaR on OS X Yosemite (10.10) and below#######

##Install Xcode using the appstore (free)

##Open a terminal and instal Xcode command line tools
xcode-select --install

##Agree to Xcode license in Terminal
sudo xcodebuild -license

##OSX doesn't use the latest version of the GCC compiler allowing to use multithreading in R. To make the installation of chromstaR
##successful, a recent version of GCC needs to be installed. This can be done using Homebrew
##Note: you may have to chown the following directories in order for homebrew to install the previous packages and library
sudo chown {your-user-name} /usr/local/bin
sudo chown {your-user-name} /usr/local/include
sudo chown {your-user-name} /usr/local/lib
sudo chown {your-user-name} /usr/local/share

##Install homebrew and the lastest versions of GCC (very long, takes 30 mins to 2 hours to compile) and CCache
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew install gcc --without-multilib
brew install clang-omp
brew install ccache
##Note: be careful in warning messages provided by homebrew and follow it's instruction/suggestions


#Modify ~/.R/Makvars in order to tell R which compiler to use. Replace everything in that file by what is below
###########################################################
CFLAGS +=             -O3 -Wall -pipe -pedantic -std=gnu99
CXXFLAGS +=           -O3 -Wall -pipe -Wno-unused -pedantic 

CC=ccache gcc
CXX=ccache g++
SHLIB_CXXLD=g++

SHLIB_OPENMP_CFLAGS = -fopenmp
SHLIB_OPENMP_CXXFLAGS = -fopenmp
SHLIB_OPENMP_FCFLAGS = -fopenmp
SHLIB_OPENMP_FFLAGS = -fopenmp


FC=ccache gfortran
F77=ccache gfortran
MAKE=make -j8

###########################################################


##Download chromstaR packages locally, as you may have to modify some of the files
##In chromstaR Package, modifies /scr/R_interface.cpp and src/scalehmm.h
##to replace #include <omp.h>  by #include <libiomp/omp.h> 
##if these corrections were already present in the files, please ignore and install directly from bitbucket (see below)


#=====================================================#
### INSTALLATION INSTRUCTIONS FOR LINUX AND WINDOWS ###
#=====================================================#

##In R, install the following packages
install.packages('devtools')
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
biocLite("GenomicAlignments")
biocLite("BSgenome")
library(devtools)
##If you have modified R_interface.cpp and scalhmm.h in chromstaR
install.packages("path-to-chromstaR-package", type='source', repos=NULL)
##If you have not modified chromstaR
install_bitbucket('chakalakka/chromstaR', username='chakalakka', password='chromstaR2015')
install_github('chakalakka/chromstaRExampleData')




