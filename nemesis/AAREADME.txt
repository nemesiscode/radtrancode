Nemesis Radiative Transfer and Retrieval Tools

Codes simulate spectra in either wavenumber or wavelength space. 

Default units for wavelength space
wavelength in microns (10-6 m)
Radiance in W cm-2 sr-1 micron-1
** Units output to .mre file are uW cm-2 sr-1 micron-1

Default units for wavenumber space
units of wavenumbers are cm-1
Radiance in W cm-2 sr-1 (cm-1)-1
** Units output to .mre file are nW cm-2 sr-1 (cm-1)-1

Transmissions are just transmissions.

There are several slightly different versions of Nemesis, the standard version
and others which deal with specific observations geometries:

Nemesis		Standard model for modelling individual observations on a 
		planet

NemesisL	As Nemesis, but optimised to deal with limb-observing 
		geometries. Model uses different method of combining individual
		layers to make the calculations faster.

NemesisMCS	Extension of NemesisL to model MCS observations of Mars. Model
		uses additional FOV data to model observations and also
		pointing data.

Nemesisdisc	Version of Nemesis for specifically modelling disc-averaged
		spectra of planets or secondary transit observations. 
		If a solar reference file is detected, then the code computes 
		the flux ratio for a secondary transit 
                If no solar reference file is detected, then the code computes
		the surface spectral irradiance of the planet in units of
		W cm-2 um-1 or W cm-2 (cm-1)-1

NemesisPT	Version of Nemesis for specifically modelling the primary 
		transit spectra of exoplanets.
		Output units are 100*planet_area/stellar_area

In addition there are some programs for generating synthetic spectra for 
retrieval testing purposes. These are:

Generatespx	General purpose program for generating test spectra

GenerateMCSspx	Version of Generatespx specifically adapted to generate
		a set of synthetic MCS observations

Pat Irwin	25/4/12 

