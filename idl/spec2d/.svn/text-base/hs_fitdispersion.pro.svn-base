;Function to measure the line dispersion function based  on the widths
;of arc lines
function hs_fitdispersion, flux, ivar, linecen

junk=' '
ratios=fltarr(n_elements(flux[0,*]))

for i=0, n_elements(flux[0,*])-1 do begin
    ngd = n_elements(linecen[0,*])
    ny=n_elements(flux[*,0])
    fwhm = fltarr(ngd)
    fwqm = fltarr(ngd)
    FOR j = 0L, ngd-1L DO BEGIN
        fwhm[j] = long_arc_fw(flux[*,i], linecen[i,j] $
                              , 7, 0.5D)/2.35 ;lets do this in sigma space
        fwqm[j] = long_arc_fw(flux[*,i], linecen[i,j] $
                              , 7, 0.25D)/(sqrt(2)*2.35) ;lets do this in sigma space
;        yy=flux[fix(linecen[i,j]-10 > 0):fix(linecen[i,j]+10 < ny-1),i]
;        xx=findgen(n_elements(yy))+fix(linecen[i,j]-10 > 0)
;        result = gaussfit(xx,yy,A) 
;        fwhm[j]=A[2]

    ENDFOR
    inmask = (fwhm GT 1.5) and (fwqm GT 1.5) and (fwqm/fwhm ge 0.8) and (fwqm/fwhm lt 1.25)
    fwhmmask=inmask
    ratios[i]=median(fwhm[where(inmask)]/fwqm[where(inmask)])
    best_sigma=fltarr(n_elements(inmask))
;    for n=0, n_elements(best_sigma)-1 do best_sigma[n]=avg([1.06*fwqm[n],fwhm[n]])
    for n=0, n_elements(best_sigma)-1 do best_sigma[n]=min([1.06*fwqm[n],fwhm[n]])
    ;; Now fit a polynomial to the fwhm and reject outliers
    fwhmset_temp1 = jfh_fit_reject(linecen[i,*], best_sigma, fwcoeff $
                                   , TOL = FWTOL $
                                   , FUNC = 'poly', YFIT = FWFIT $
                                   , outmask = fwhmmask, inmask = inmask $
                                   , XMIN = 0, XMAX = ny-1L)
                                  
    good_lines = WHERE(fwhmmask GT 0, ngline)
;    splot, linecen[i,good_lines], fwhm[good_lines], ps=1
;    soplot, linecen[i,good_lines], fwqm[good_lines], ps=1, color='red'
;    soplot, linecen[i,good_lines], best_sigma[good_lines], ps=1, color='green'
;    soplot, fwhmset_temp1.coeff[0]+linecen[i,*]*fwhmset_temp1.coeff[1]+linecen[i,*]^2.0*fwhmset_temp1.coeff[2], color='red'
;    soplot, fwhmset_temp1.coeff[2]+linecen[i,*]*fwhmset_temp1.coeff[1]+linecen[i,*]^2.0*fwhmset_temp1.coeff[0], color='green'
 ;   read, junk
        fwhm=best_sigma

    IF ngline NE 0 THEN BEGIN
        fwhm_med = djs_median(fwhm[good_lines])
    ENDIF ELSE BEGIN
        fwhm_med = 0.0
    ENDELSE

    if size(fwhmset_temp1,/type) eq 8 then begin
      fwhmset_temp = struct_addtags(fwhmset_temp1 $
                                  , create_struct('MEDIAN', fwhm_med))
    endif else begin
      fwhmset_temp={func:'poly', xmin:0, $
         xmax:(ny-1L), coeff:fltarr(3), median:fwhm_med} 
    endelse

    if i eq 0 then begin
        fwhmset={func:fwhmset_temp.func, xmin:fwhmset_temp.xmin, $
         xmax:fwhmset_temp.xmax, coeff:fltarr(n_elements(fwhmset_temp.coeff),$
         n_elements(flux[0,*])), median:fltarr(n_elements(flux[0,*]))} 
        fwhmset.median[0]=fwhmset_temp.median
        fwhmset.coeff[*,0]=fwhmset_temp.coeff
    endif else begin
        fwhmset.coeff[*,i]=fwhmset_temp.coeff
        fwhmset.median[i]=fwhmset_temp.median
    endelse
endfor
;xy2traceset, temp1, temp, tempset, ncoeff=6
;print, median(ratios)
;stop
return, fwhmset
end
