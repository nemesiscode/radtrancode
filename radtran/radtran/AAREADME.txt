Radiative transfer calculation routines for performing LBL, band and c-k 
calculations with Nemesis. A number of subroutines are also used by cirsrad,
cirsradg and nemesis.

Executables:

Aground 	Adds ground emission times transmission to atmospheric
                emission.

Ave_table	Reads an absorption coefficient look-up table and
                convolves with a square of user defined width.

Calc_fnktable	Calculates a k-table from line data.

Calc_exotable3	Calculates a k-table from a range of line databases calculated
		at different temperatures by Makedbloop

Calc_exotable4	Calc_exotable is ALMOST same as calc_exotable3. However, to
		use Calc_exotable4, we need to create a blank k-table first 
		and then run the code and fill it up with k-coefficients for 
		each grid point. However,this process is totally manual, not 
		automated, so is probably not a usable code for all users. 

Aband_ktable	Calculates a k-table from band data over a regular grid of
		of wavelengths/wavenumbers

Aband_ktablec	Calculates a k-table from band data at a set of output 
		wavelengths/wavenumbers read in from user.

Conv_spec	Program to smooth Radtrans spectra.

Convert_table 	Program to convert Radtran .kta binary k-coefficient
		files into more easily read ASCII .par format

Direct_ktable	Converts ASCII output of Aband_ktable to direct-access

Direct_ktablec	Converts ASCII output of Aband_ktablec to direct-access

Par_ktable	Program to convert a band_ktable ASCII output file into a
	     	more formatted .par file

Li_spec		Lists standard binary output of Radtrans.

Pl_spec		Plots the binary data output by Radtrans

Radtrans	Top level RT code.

Read_table	Reads in and examines a k-table.

Concat_table	Concatenates two k-tables together.

Cut_table	Extracts specified wavelength range from a k-table



Pat Irwin	25/4/12
