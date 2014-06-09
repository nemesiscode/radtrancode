;+
; NAME:
;	Gauss_Leg_Quadr
;
; PURPOSE:
;	Setup Gaussian-Legendre quadrature for numerical integration.
;	Get roots of Legendre polynomials and corresponding
;	weights for Gaussian quadrature of a function:
;	Approximate integral = total( weights * function( roots ) )
;
; CALLING:
;	Gauss_Leg_Quadr, Npq, xq, wq
;
; INPUTS:
;	Npq = integer, number of points desired for quadrature,
;		forced to be an even number. The quadrature will be
;		exact for polynomials up to degree 2*Npq-1.
;
; KEYWORDS:
;	XRANGE = range of integration, default = (-1,1).
;		Range can be changed later by rescaling points and weights,
;		(see code at end of this routine).
;	EPSILON = accuracy of root convergence, default = 3.d-14.
;
; OUTPUTS:
;	xq = points at which function should be evaluated.
;	wq = corresponding weights:  Integral = total( wq * func( xq ) )
;
; EXAMPLE:
;	IDL> Gauss_Leg_Quadr, 2, xq, wq, XRANGE=[0,3]
;	IDL> print, total( wq * xq^2 )
;	       9.0000000		; = integral of x^2 from 0 to 3.
;	IDL> print, total( wq * xq^3 )
;	       20.250000		; = integral of x^3 from 0 to 3.
;
; PROCEDURE:
;	G. Rybicki algorithm as described in Numerical Recipes section 4.5.
;	Uses Newton's method to find roots of Legendre polynomials,
;	since the Leg.Poly. derivatives are needed anyway to compute weights.
; HISTORY:
;	Written: Frank Varosi NASA/GSFC 1994.
;	F.V.1995: fixed mistake in using xrange.
;-

pro Gauss_Leg_Quadr, Npq, xq, wq, XRANGE=xrange, EPSILON=eps

	if N_elements( Npq ) NE 1 then begin
		print,"syntax:"
		print,"	Gauss_Leg_Quadr, Npoint, xpoints, weights, XRANGE=[ , ]"
		return
	   endif

	Npq = ( 2 * fix( (Npq+1)/2 ) ) > 2
	Nroot = Npq/2
	if N_elements( eps ) NE 1 then eps = 3.d-14
	zr = cos( !DPI * (dindgen( Nroot ) + 0.75)/(Npq+0.5) )
	zr1 = zr

	p1 = make_array( Nroot, /DOUBLE )
	p2 = p1
	p3 = p1
	pp = p1
	w = indgen( Nroot )

	REPEAT BEGIN	;apply Newton's method selectively until all converged.

		p1(w) = 1
		p2(w) = 0

		for j = 1, Npq do begin
			p3(w) = p2(w)
			p2(w) = p1(w)
			p1(w) = ( (2*j-1)*zr(w)*p2(w) - (j-1)*p3(w) ) /j
		  endfor

		pp(w) = Npq * ( zr(w)*p1(w) - p2(w) ) / ( zr(w)^2 - 1 )
		zr1(w) = zr(w)
		zr(w) = zr(w) - p1(w)/pp(w)

		w = where( abs( zr - zr1 ) GT eps, newton )

	 ENDREP UNTIL (newton LE 0)

	xq = [ -zr, rotate( zr, 2 ) ]
	wq = 2/( (1-zr^2) * pp^2 )
	wq = [ wq, rotate( wq, 2 ) ]

	if N_elements( xrange ) EQ 2 then begin
		xL = ( xrange(1) - xrange(0) )/2.0
		wq = xL * wq
		xq = xL * xq + total( xrange )/2.0
	   endif
end

