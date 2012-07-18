Reference files used by different Radtrans/Nemesis subroutines and also a 
convenient place to store reference CIA tables and reference solar/stellar
spectra.

There are a number of CIA tables in here, for which there are two formats: 
a) standard, which lists 9 pairs of CIA coefficients
b) para-H2, which lists H2-H2 and H2-He only, but for variable ortho/para
    ratios.

In addition, CIA tables are pre-calculated with a grid spacing of usually 
1 cm-1 or 10 cm-1. All tables have 25 temperatures and are calculated at 1501
wavenumbers. You just have to know!

Here's a list of the current CIA tables:

a) standard
NewBorysow.tab			50-410K		dnu=10.
Borysowsub.tab			50-410K		dnu=10. (extended to < 1 micron)
GlennBorysow.tab		50-410K		dnu=10.
isotest.tab			70-310K 	dnu=1.
isotest_lowt.tab		57.5-285.5K	dnu=1.
glenniso.tab			70-310K 	dnu=1.
exociatable_fletcher.tab	400-2800K	dnu=10.
exociatable.tab			400-2800K	dnu=10.

b) para-H2
fparah2he.tab			50-410K		dnu=1.0	npara=12 df=0.02
fpara24_dn1.0_40-280K_df0.02.tab 40-280K	dnu=1.0 npara=24 df=0.02
fpara24_dn1.0_40-280K_df0.03.tab 40-280K	dnu=1.0 npara=24 df=0.03
magnus_dn1.0_40-280K_df0.03.tab 40-280K		dnu=1.0 npara=24 df=0.02

Pat Irwin  25/4/12	Original
Pat Irwin  18/7/12	Revised

