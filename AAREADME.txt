Main Radtrans/Nemesis development repository.

The directories are:

frecipes/	Numerical recipes Fortran programs. A makefile is included.

raddata/	Reference spectral data for use by Radtrans/Nemesis, including
   		CIA tables, solar/stellar reference spectra, gas 
		continuum spectra, lineshape factors, SVP curves and last,
		but definitely not least, a list of the current planets 
		defined, together with their mass and radius, etc.

radtran/	Main Radtrans directory, containing a number of sub-directories.
		Several executables are defined, including Radtrans itself and
		a number of line data and correlated-k processing programs.
		The subdirecoties makefiles/ includes shell scipts for making 
		the libraries in one go and for making all the Radtrans
		executables. The makefiles for each sub-directory are also
		included.

nemesis/	The Nemesis subroutines and programs. The Nemesis programs are
		built upon the Radtrans libraries in radtran/ so these must
		be compiled first.

idl/		Useful IDL programs for reading and plotting Radtrans/Nemesis
		calculations.

manuals/	The Radtrans and Nemesis manuals.

FOVgreg/	Library of routines for defining MCS Field of View. Only needs 
		compiling if you're going to run NemesisMCS or GenerateMCSspx.

Pat Irwin	30/4/12
