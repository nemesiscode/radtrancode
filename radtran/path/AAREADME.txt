Subroutines for dealing with paths and a number of executables that manipulate
profiles.

Path		Program for reading a <runname>.pat file and calculating the 
		layer temperatures, amounts etc. and output to the 
		<runname>.drv metafile.

Profile		Program for reading in and manipulating profile files 
		(i.e. <runname>.prf files)

Parah2_profile	Calculates a para-H2 profile for use in Nemesis.

Dust_profile	Generates dust specific concentration profile file
		where specific concentration is the number of dust
		particles per gram of atmosphere.

Nemesisprofile	Program for reading in a <runname.prf> file and extracting 
		a T or vmr profile for use as an apriori profile for Nemesis.

Mod_profile	Program for scaling the opacities in a dust profile.

Interp_prf	Program for interpolating a <runname>.prf file onto a new
		log(pressure) grid.

Write_path	Program for leading the user in setting up a <runname>.pat
		file.

Pat Irwin	24/4/12
