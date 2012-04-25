Subroutines and programs for manipulating line data, band data and 
correlated-k data.

Executables:

Makeband	Creates a band data base file (.ban) using line data.

Mergeband	Merges two .ban files together.

Modifyband	Converts old .ban format files to extended scheme required
		for Kam Sihra's CH4 measurements.

Makedb		Makes Radtran lbl direct access data base file from a
		variety of sources.

Merge		Merges two sequential access line database files together.

Scan		Lists line data to screen - ascii as stored in file

Li_lines	Lists lines from radtran data bases to screen

Pl_lines	Outputs lines to ASCII format file plottable by IDL
		routine pl_lines.pro

Select		Copies a subset of a line data base to a new data base
		and allows modification of the lines.

Cp_lines	As Select, but no line modification allowed.

Summary		Outputs a summery of numbers of lines and strengths
     		for all species in given intervals

Write_xsec	Simple program to write out a dust cross section file

Inv_scan	Modified version of Scan to just list inversion lines.

Makedbloop	Modified from MAKEDB to run over a range of input files
   	        (assumed generated from Selecttemploop) calculated over a range
		of temperatures. Generates a separate line database for each
		temperature.
 
Sortlines	Program to sort a sequential-access linedata file

Pl_elines	Modification of Pl_lines that also prints out the
		lower state energies of the lines.
	
Selecttemp	Revised version of Select to select lines relevant for
               	a given temperature.

Selecttemploop	Revised version of Select to extract appropriate lines for a 
		set of different temperatures.

Pat Irwin	24/4/12
