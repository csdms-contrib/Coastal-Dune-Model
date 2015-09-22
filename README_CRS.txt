Notes on building running Coastal Dune Model on Ubuntu with gfortran

csherwood@usgs.gov
21 Sept 2015

Required library:
   Needed to download, build, and install fftw
   http://www.fftw.org/
   
     ./configure
     make
     sudo make install

     ...worked fine

Changes to Makefile:
     Remove (or change) the line that copies exectutable to bin directory

Changes to code:
     My version of gfortran required that in globals.cc, I add:
     # include <cstring>

To do a test run:
     mkdir test
     cp Dune test/.
     cp data/input/*.par test/default.par
     cd test
     mkdir DATA
     ./Dune < default.par

     
     
