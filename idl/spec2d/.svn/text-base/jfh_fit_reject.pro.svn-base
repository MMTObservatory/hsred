;;+
; NAME:
;      jfh_fit_reject
;
; PURPOSE:
;      Fit dispersion function with proper outlier rejection
;
; MODIFICATION HISTORY:
;   original code part of XIDL package by J. Hennawi, X. Prochaska et al. 
;   copied and included in hsred distribution for convenience

FUNCTION jfh_fit_reject, x, y, ncoeff, inmask = inmask, outmask = outmask $
                         , ivar = ivar1, maxrej = maxrej $
                         , TOL = TOL, YFIT = YFIT, SIGREJ = SIGREJ $
                         , FUNC = FUNC, XMIN = XMIN, XMAX = XMAX

IF NOT KEYWORD_SET(SIGREJ) THEN SIGREJ = 3.0
IF NOT KEYWORD_SET(TOL) THEN TOL = 1.0e10
nx = n_elements(x)
IF KEYWORD_SET(INMASK) THEN OUTMASK = INMASK $
ELSE outmask = lonarr(nx)+1
use_inds = where(outmask, nuse)
xset = 0
IF NOT KEYWORD_SET(MAXREJ) THEN maxrej = nuse/2
; Now fit a polynomial to the fwhm and reject outliers
FOR j = 0, maxrej-1L DO BEGIN
    IF KEYWORD_SET(ivar1) THEN ivar = ivar1[use_inds]
;    IF FUNC EQ 'GRISM' THEN $
;      xset = grism_solve(x[use_inds], y[use_inds], ivar $
;                         , pix0 = ny/2, pixfit = yfit) $
;    ELSE 
    xy2traceset, x[use_inds], y[use_inds], xset, invvar = ivar $
                 , ncoeff = ncoeff, maxiter = 0, /silent, FUNC = FUNC $
                 , XMIN = XMIN, XMAX = XMAX $
                 , yfit = yfit
    resids = y[use_inds]-yfit
    sig = djsig(resids)
    rej_inds = WHERE(abs(resids) GE sigrej*sig OR abs(resids) GE TOL, nrej)
    IF nrej NE 0 THEN BEGIN
        max_resid = max(abs(resids[rej_inds]), jbad)
        outmask[use_inds[rej_inds[jbad]]] = 0
        use_inds = WHERE(outmask, nkeep)
    ENDIF ELSE BREAK
ENDFOR
RETURN, xset
END
