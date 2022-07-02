PRO mytv,image,param2,param3,param4, $
         range=range,ncol=ncol,noscale=noscale,table=table,_extra=_extra
;
;================================================================
; Display an image either on true colour or index colour display.
;
; MYTV, IMAGE, ... , [KEYWORDS]
;
; KEYWORDS:
;
;    TABLE: Set to [NCOL,3] element array with required colour tables,
;           second dimension is in order red, green, blue.
;
;           If not specified, assumes that you have loaded a colour table,
;           e.g. with loadct.
;
;    RANGE: Image values are in range 0 to RANGE-1.
;           If not specified, size is taken from TABLE array.
;           If neither RANGE nor TABLE is specified, 256 is assumed.
;
;    NCOL:  How many colours of the colour table to use.
;           If not specified, whole TABLE array is used.
;           If neither NCOL nor TABLE is specified, !d.table_size is assumed.
;
;    NOSCALE: Don't scale image values at all.
;             Equivalent to setting both RANGE and NCOL to the same value.
;
;    Other keywords and parameters: passed to TV command.
;
; Alan Iwi - 30/01/2001
;
;
; USAGE EXAMPLES:
;
;     (1)   GIF image: colour table passed to MYTV
;
;                 read_gif,'my.gif',image,r,g,b
;                 mytv,image,table=[[r],[g],[b]]
;
;     (2)   GIF image: colour table loaded beforehand
;
;                 read_gif,'my.gif',image,r,g,b
;                 tvlct,r,g,b
;                 mytv,image,/noscale
;
;     (3)   Image in range 0-255, colour table loaded
;
;                 image=replicate(1,200)#indgen(256)
;                 loadct,13
;                 mytv,image
;
;     (4)   Image in range 0-100, colour table loaded
;
;                 image=replicate(1,200)#indgen(100)
;                 loadct,13
;                 mytv,image,range=100
;
;     (5)   Image in range 0-255, but only 10-colour table loaded
;
;                 image=replicate(1,200)#indgen(256)
;                 loadct,13,ncol=10
;                 mytv,image,ncol=10
;
;     (6)   Example with other parameters understood by TV command,
;           e.g. positioning on page
;
;                 image=replicate(1,200)#indgen(256)
;                 window,xsize=400,ysize=256
;                 loadct,13
;                 mytv,image,100,0
;
;================================================================

; open window if not already open
if !d.name eq 'X' and !d.window eq -1 then begin
    s=size(image)
    window,xsize=s[1],ysize=s[2]
endif

;-----------
; work out number of input and output colours, and scale image values if 
; these don't match
;
tabsize=(size(table))[1]

if n_elements(range) eq 0 then begin
    if tabsize ne 0 then range=tabsize else range=256
endif

if n_elements(ncol) eq 0 then begin
    if tabsize ne 0 then ncol=tabsize else ncol=!d.table_size
endif

if range ne ncol and not keyword_set(noscale) then begin
    myimage = byte(image * float(ncol)/(range-1)) < (ncol-1)
endif else begin
    myimage=image
endelse
;-----------

if (!d.n_colors eq 2L^24) then begin
    
    ;; TRUE COLOUR
    
    if tabsize eq 0 then begin        
        tvlct,r,g,b,/get
    endif else begin
        r=table[*,0]
        g=table[*,1]
        b=table[*,2]
    endelse
    
    myimage=[[[r[myimage]]],[[g[myimage]]],[[b[myimage]]]]
    true=3
    
endif else begin
    
    ;; INDEX COLOUR
    
    if tabsize ne 0 then tvlct,table
    true=0
    
endelse

;; call TV, with the right number of parameters
case n_params() of
    1: tv,myimage,true=true,_extra=_extra
    2: tv,myimage,param2,true=true,_extra=_extra
    3: tv,myimage,param2,param3,true=true,_extra=_extra
    4: tv,myimage,param2,param3,param4,true=true,_extra=_extra
endcase

end
