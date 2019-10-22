#########################################################################################################
                                                   NemesisPY
#########################################################################################################

NemesisPY is the python-written version of Nemesis

In order to install it, follow the steps:
  
    (1) Compile the fortran code as regular for using Nemesis
	It is not really necessary to compile the programs in the nemesis/ directory, but it is required to
	compile the programs in radtran/, as CIRSdrv_wavePY is required. All libraries must be compiled
 
    (2) Put the nemesisSO.comp file under the nemesis/ directory
        It usually compiles better using gfortran. In order to set up correctly the compiler environment type in the terminal:
          setenv F90 gfortran
          setenv F77 gfortran
	  setenv CC gcc
        Go to the nemesisPY/ directory and type: source nemesisf.comp 

    (3) Look for the function check_arraysize_nemesis() in the nemesis.py library, and change the numbers as used in the Fortran
        programs. This is done for reconciling the sizes in the Python and Fortran routines, when needed.


    (4) To make the python scripts executable from anywhere, do:
    	- Write in the very first line of the python pogram:
          #!/usr/bin/python
          
        - Make your program executable by typing in the terminal:
          chmod +x name.py

        - Move your program to the folder in which the Nemesis programs (bin folder) are stored 

        - Give an alias to your program in your .tcshrc file. For example:
          alias nemesisSO $BIN/nemesisSO.py
