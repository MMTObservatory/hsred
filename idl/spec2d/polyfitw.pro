; $Id: polyfitw.pro,v 1.6 2000/06/19 17:48:58 chris Exp $
;+
; NAME:
;POLYFITW
;
; PURPOSE:
;Perform a least-square polynomial fit with optional error estimates.
;
; CATEGORY:
;Curve fitting.
;
; CALLING SEQUENCE:
;Result = POLYFITW(X, Y, W, NDegree [, Yfit, Yband, Sigma, Corrm]
;          [, /DOUBLE] [, STATUS=status])
;
; INPUTS:
;    X:The independent variable vector.
;
;    Y:The dependent variable vector.  This vector should be the same
;length as X.
;
;    W:The vector of weights.  This vector should be same length as
;X and Y.
;
;     NDegree:The degree of polynomial to fit.
;
; OUTPUTS:
;POLYFITW returns a vector of coefficients of length NDegree+1.
;
; OPTIONAL OUTPUT PARAMETERS:
; Yfit:The vector of calculated Y's.  Has an error of + or - Yband.
;
;Yband:Error estimate for each point = 1 Sigma.
;
;Sigma:The experimental standard deviation in Y units.
;
;Corrm:Correlation matrix of the coefficients.
;
; KEYWORDS:
;   DOUBLE: Set this keyword to force computations to be done in
;           double-precision arithmetic.
;
;   STATUS = Set this keyword to a named variable to receive the
;status
;            of the operation. Possible status values are:
;             0 = successful completion
;             1 = singular array (invalid inversion)
;             2 = warning that a small pivot element was used
;             3 = undefined Yband error estimate
;
;    Note: if STATUS is not specified then any error messages will be
;output
;          to the screen.
;
; COMMON BLOCKS:
;None.
;
; SIDE EFFECTS:
;None.
;
; MODIFICATION HISTORY:
;Written by: George Lawrence, LASP, University of Colorado,
;December, 1981.
;
;Adapted to VAX IDL by: David Stern, Jan, 1982.
;
;Weights added, April, 1987,  G. Lawrence
;
;   CT, RSI, March 2000: Added STATUS keyword.
;-

FUNCTION POLYFITW,x,y,w,ndegree,yfit,yband,sigma,coorm, $
                  DOUBLE=double,STATUS=status

COMPILE_OPT strictarr

ON_ERROR,2                      ;return to caller if an error occurs

n = N_ELEMENTS(x)
IF (n NE N_ELEMENTS(y)) THEN MESSAGE, $
   'X and Y must have same number of elements.'
m = ndegree + 1; # of elements in coeff vec

IF (N_ELEMENTS(double) EQ 0) THEN $
   double = (SIZE(x,/TNAME) EQ 'DOUBLE') OR (SIZE(y,/TNAME) EQ 'DOUBLE') $
ELSE $
   double = KEYWORD_SET(double)

no_weight = (N_ELEMENTS(w) EQ 1) AND (w[0] EQ 1) ; used for call from POLY_FIT()


; construct work arrays
IF (double) THEN BEGIN
   coorm = DBLARR(m,m) ; least square matrix, weighted matrix
   b = DBLARR(m); will contain sum w*y*x^j
   z = DBLARR(n) + 1; basis vector for constant term
   wy = DOUBLE(w)*DOUBLE(y)
   yfit = DBLARR(n)
   yband = DBLARR(n)
   sigma = !VALUES.D_NAN
ENDIF ELSE BEGIN
   coorm = FLTARR(m,m) ; least square matrix, weighted matrix
   b = FLTARR(m); will contain sum w*y*x^j
   z = FLTARR(n) + 1.; basis vector for constant term
   wy = FLOAT(w)*FLOAT(y)
   yfit = FLTARR(n)
   yband = FLTARR(n)
   sigma = !VALUES.F_NAN
ENDELSE


coorm[0,0] = no_weight ? n : TOTAL(w, DOUBLE=double)
b[0] = TOTAL(wy)


FOR p = 1L,2*ndegree DO BEGIN; power loop
   z = TEMPORARY(z)*x; z is now x^p
   IF p LT m THEN b[p] = TOTAL(wy*z); b is sum w*y*x^j
   sum = no_weight ? TOTAL(z) : TOTAL(w*z, DOUBLE=double)
   FOR j = 0 > (p-ndegree), ndegree < p DO coorm[j,p-j] = sum
ENDFOR ; end of p loop, construction of coorm and b


coorm = INVERT(TEMPORARY(coorm), status)
IF NOT ARG_PRESENT(status) THEN BEGIN
   CASE status OF
      1: MESSAGE, "Singular matrix detected."
      2: MESSAGE,/INFO, "Warning: Invert detected a small pivot element."
      ELSE:
   ENDCASE
ENDIF
IF (status EQ 1) THEN RETURN, double ? !VALUES.D_NAN : !VALUES.F_NAN


c = (TEMPORARY(b) # coorm)  ; construct coefficients


; compute optional output parameters.

; one-sigma error estimates, init
yfit = TEMPORARY(yfit) + c[ndegree]
FOR k = ndegree-1L, 0, -1 DO yfit = c[k] + TEMPORARY(yfit)*x  ; sum basis vectors


; Experimental variance estimate, unbiased
var = (n GT m) ? TOTAL((yfit-y)^2 )/(n-m) : ((double) ? 0d : 0.0)

sigma = SQRT(var)
z = double ? DBLARR(n) + 1 : FLTARR(n) + 1
yband = TEMPORARY(yband) + coorm[0,0]
FOR p=1L,2*ndegree DO BEGIN; compute correlated error estimates on y
   z = TEMPORARY(z)*x; z is now x^p
   sum = 0
   FOR j=0 > (p - ndegree), ndegree  < p DO sum = sum + coorm[j,p-j]
   yband = TEMPORARY(yband) + sum * z  ; add in all the error sources
ENDFOR; end of p loop


yband = TEMPORARY(yband)*var
IF (MIN(yband) LT 0) OR (MIN(FINITE(yband)) EQ 0) THEN BEGIN
   status = 3
   IF NOT ARG_PRESENT(status) THEN MESSAGE, $
      'Undefined (NaN) error estimate encountered.'
ENDIF ELSE yband = SQRT( TEMPORARY(yband) )

RETURN,c
END


