function monotonic, x
;+
; NAME:
;	MONOTONIC
;
; PURPOSE:
;	Determine whether a vector is monotonically increasing or
;	decreasing. 
;
; CALLING SEQUENCE:
;	truefalse = monotonic(x)
;
; INPUTS:
;	x - one-dimensional vector of any datatype
;
; OUTPUTS:
;	Returns 1 if the vector is monotonically increasing or
;	decreasing, and 0 otherwise.
;
; COMMENTS:
;	Does not distinguish between increasing and decreasing
;	vectors.  Could be extended to two-dimensional vectors
;	easily. 
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 July 23, U of A, written
;-

    nx = n_elements(x)
    vec = lindgen(nx)   ; comparison vector
    srt = sort(x)       ; sort the vector       
    rsrt = reverse(srt) ; reverse the sorted vector       
    
    true = 0L
    if total(abs(srt-vec)) eq float(0) then true = 1L else $ ; monotonically increasing
      if total(abs(rsrt-vec)) eq float(0) then true = 1L     ; monotonically decreasing
    
return, true
end
