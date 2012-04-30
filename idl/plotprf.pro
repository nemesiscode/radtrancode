;****************************************************************************** 
;_TITL:	plotprf.pro
;
;_DESC:	Procedure to plot IDL .prf files written out by profile.f
;
;_ARGS:	NP	INTEGER		Number of points in profile.
;	NCONT	INTEGER		Number of cloud decks.
;
;_HIST:	PGJI	23/9/96		ORIGINAL VERSION.
;	PDP	18mar2002	Major overhaul of the code: updated the code to
;				the modern/current format for .prf files;
;				plotted height and pressure on the same plot
;				using the axis command; provided the user with
;				an interactive dialog box for selecting an
;				input .prf file; added the comments; removed
;				references to obsolete programs like plot_io,
;				plot_oo and plot_oi.
;	PDP	20aug2002	Added the option to plot dust/aerosol profiles.
;				Also capitalised all the IDL commands to make
;				the code easier to read.
;******************************************************************************
;******************************************************************************

answer=' '
PRINT,' '
PRINT,' Select "a" to plot a temperature & compostion profile, else'
READ, Answer, PROMPT= 'select "b" to plot a dust/aerosol profile: '
PRINT,' '
IF (answer EQ 'a') OR (answer EQ 'A') THEN iprf= 0
IF (answer EQ 'b') OR (answer EQ 'B') THEN iprf= 1
IF (iprf NE 0) AND (iprf NE 1) THEN BEGIN
  PRINT, 'Not valid entry. Stopping program.'
  STOP
ENDIF

filename= DIALOG_PICKFILE(/READ, TITLE='Select a profile', FILTER='*.prf')
OPENR,1,filename


;==============================================================================
;
;	Plot a temperature & composition profile ...
;
;==============================================================================
IF iprf EQ 0 THEN BEGIN

;==============================================================================
;=========================== Read the header and gas identification information
  header=' '
  READF,1,header
  WHILE (STRMID(header,0,1) eq '#') DO READF,1,header
  iform= FIX(header)
  READF,1,iplanet,latitude,np,nvmr,molwt
  idgas= INTARR(nvmr) & isogas= INTARR(nvmr)
  temp= np & height= np & press= np
  FOR i=0, nvmr-1 DO BEGIN
    READF,1,a,b
    idgas(i)= a
    isogas(i)= b
;    PRINT,i,idgas(i),isogas(i)
  ENDFOR

;==============================================================================
;==================== Read the VMRs. Did not account for profiles with iform= 1
;======================== where there is no "press" column. Also hard-coded the
;============================= program to expect 11 gases (although that should
;================================== be easy to correct if one used a while loop
;============================================= in combination with a for loop).
  height= FLTARR(np)
  press= FLTARR(np)
  temp= FLTARR(np)
  vmr= FLTARR(nvmr,np)
  READF,1,header
  FOR i=0, np-1 DO BEGIN
    READF,1,a,b,c,d,e,f
    height(i)= a
    press(i)= b
    temp(i)= c
    vmr(0,i)= d
    vmr(1,i)= e
    vmr(2,i)= f
  ENDFOR
  READF,1,header
  FOR i=0, np-1 DO BEGIN
    READF,1,a,b,c,d,e,f
    vmr(3,i)= a
    vmr(4,i)= b
    vmr(5,i)= c
    vmr(6,i)= d
    vmr(7,i)= e
    vmr(8,i)= f
  ENDFOR
  READF,1,header
  FOR i=0, np-1 DO BEGIN
    READF,1,a,b
    vmr(9,i)= a
    vmr(10,i)= b
  ENDFOR

  CLOSE,1

;==============================================================================
;============= Plot the temperature and VMR profiles. The axis command plots on
;=================== the second y-axis. The bits at the bottom were kept should
;===================================== the user desire further/different plots.
  !P.MULTI=[0,2,1]

  hmin= MIN(height)
  hmax= MAX(height)
  pmin= MIN(press)
  pmax= MAX(press)
  PLOT, temp, height, YRANGE= [hmin,hmax], YSTYLE= 1, $
	XTITLE= 'Temperature [K]', YTITLE= 'Altitude [km]'
  AXIS, YAXIS= 1, YRANGE= [pmax,pmin], /YLOG, YSTYLE= 1, $
	YTITLE= 'Pressure [atm]  |  Altitude [km]'
  PLOT, vmr(0,*), height, XRANGE= [1E-12,1E-2], /XLOG, $
	YRANGE= [hmin,hmax], YSTYLE= 1, $
	XTITLE= 'Volume Mixing Ratio', LINESTYLE= 0
  AXIS, YAXIS= 1, YRANGE= [pmax,pmin], /YLOG, YSTYLE= 1, $
	YTITLE= 'Pressure [atm]'
  FOR i= 1, nvmr-1 DO OPLOT, vmr(i,*), height, LINESTYLE= i

  cap_1= 'NH!D3!N'
  cap_2= 'PH!D3!N'
  cap_3= 'C!D2!NH!D2!N'
  cap_4= 'C!D2!NH!D4!N'
  cap_5= 'C!D2!NH!D6!N'
  cap_6= 'C!D4!NH!D2!N'
;  cap_7= 'H!D2!N'		; commented out because with such large
;  cap_8= 'He'			; values, their inclusion ruins the plot.
  cap_9= 'CH!D4!N'
  cap_10= '!E13!NCH!D4!N'
  cap_11= 'CH!D3!ND'
  XYOUTS, 0.9E-4, 15, cap_1, /DATA
  XYOUTS, 9.0E-7, -10, cap_2, /DATA
  XYOUTS, 1.8E-9, 80, cap_3, /DATA
  XYOUTS, 2.2E-10, 130, cap_4, /DATA
  XYOUTS, 9.0E-7, 70, cap_5, /DATA
  XYOUTS, 0.8E-11, 340, cap_6, /DATA
;  XYOUTS, 0.18, -30, cap_7, /DATA
;  XYOUTS, 2.0E-2, -30, cap_8, /DATA
  XYOUTS, 2.5E-4, 300, cap_9, /DATA
  XYOUTS, 2.4E-5, 60, cap_10, /DATA
  XYOUTS, 4.5E-9, 320, cap_11, /DATA


;==============================================================================
;================================================================= Old code ...
;ans= ' '
;read, ans
;
;pmin= min(press)
;pmax= max(press)
;plot, temp, press, yrange= [pmax,pmin], /ylog, $
;	xtitle= 'Temperature [K]', ytitle= 'Pressure [atm]'
;plot, vmr(0,*), press, xrange=[1e-12,1], yrange=[pmax,pmin], /xlog, /ylog, $
;	xtitle= 'Volume Mixing Ratio', ytitle='Pressure [atm]',linestyle= 0
;for i= 1, nvmr-1 do oplot, vmr(i,*), press, linestyle= i
;
;ans= ' '
;read, ans
;
;logpmin= min(alog10(press))
;logpmax= max(alog10(press))
;plot, temp, alog10(press), yrange=[logpmax,logpmin], $
;	xtitle= 'Temperature [K]', ytitle= 'Log!D10!N Pressure [atm]'
;plot, vmr(0,*), alog10(press), xrange= [1e-12,1], /xlog, $
;	yrange= [logpmax,logpmin], linestyle= 0, $
;	xtitle= 'Volume Mixing Ratio', ytitle='Log!D10!N Pressure [atm]'
;for i= 1, nvmr-1 do oplot, vmr(i,*), alog10(press), linestyle= i


;==============================================================================
;
;	Plot a dust/aerosol profile ...
;
;==============================================================================
ENDIF ELSE BEGIN

;==============================================================================
;================================================= Read the header information.
  header=' '
  READF,1,header
  WHILE (STRMID(header,0,1) eq '#') DO READF,1,header
  np= STRMID(header,0,3) & ncont= STRMID(header,4,5)

;==============================================================================
;=== Read the cloud number densities and approximate the corresponding pressure
;========== for plotting (reference pressure [atm] and height [km] set equal to
;============================ the tropopause). Hard-coded the program to expect
;=========================================================== three cloud decks.
  height= FLTARR(np)
  press= FLTARR(np)
  cloud= FLTARR(ncont,np)

  p0= 0.1 & z0= 45 & scaleh= 18		; [atm], [km] and [km] respectively
  FOR i=0, np-1 DO BEGIN
    READF,1,a,b,c,d
    height(i)= a
    press(i)= p0*EXP(-(height(i) - z0)/scaleh)
    cloud(0,i)= b
    cloud(1,i)= c
    cloud(2,i)= d
  ENDFOR

  CLOSE, 1

;==============================================================================
;== Plot the dust/aerosol profile. The axis command plots on the second y-axis.
  !P.MULTI=[0,1,1]

  hmin= MIN(height)
  hmax= MAX(height)
  pmin= MIN(press)
  pmax= MAX(press)
  PLOT, cloud(0,*), height, XRANGE= [1E-12,1E-2], /XLOG, $
	YRANGE= [hmin,hmax], YSTYLE= 1, $
	XTITLE= 'Cloud number density [/cm^3]', YTITLE= 'Altitude [km]', $
	LINESTYLE= 0
  AXIS, YAXIS= 1, YRANGE= [pmax,pmin], /YLOG, YSTYLE= 1, $
	YTITLE= 'Pressure [atm]'
  FOR i= 1, ncont-1 DO OPLOT, cloud(i,*), height, LINESTYLE= i

ENDELSE

!P.MULTI=[0,0,0]

END
