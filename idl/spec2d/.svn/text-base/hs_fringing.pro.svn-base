; procedure to calculate and apply a fringing correction in the red
pro hs_fringing, rerun

;wdir=DEF_WDIR(LOGFILE)
;instrument=DEF_INST(LOGFILE)
;flat=readfits(wdir+'flat_s.fts',head_flat,/silent)
flat1=mrdfits('calibration/'+rerun+'/traceflat.fits',0,head_flat)
flat2=mrdfits('calibration/'+rerun+'/traceflat.fits',1)
for i=0, 149 do flat1[*,i] = flat1[*,i]/median(flat1[*,i])
for i=0, 149 do flat2[*,i] = flat2[*,i]/median(flat2[*,i])
flat=[[flat1],[flat2]]
head_image=head_flat
;maskname=(instrument eq 'GMOS')? strcompress(sxpar(head_image,'MASKNAME'),/remove_all) : ''
;binc=long(strsplit(sxpar(head_image,'CCDSUM'),' ',/extract))
;if(n_elements(binc) ne 2) then 
binc=[1,1]

;;; get rid of vignetting
a=size(flat)
if a(0) eq 0 then return
nx=a[1]
n_fib=a[2]
x=findgen(nx)

;kk=readfits(wdir+'neon_distortion.fts')
kk=fltarr(n_fib)

;if((instrument eq 'GMOS' and n_fib eq 1500 and $
;        (maskname eq 'IFU' or maskname eq 'IFU-R'))) then $
;    flat[950:*,750:*]=flat[950:*,750:*]*1.06634
;for i=0,n_fib-1 do flat[*,i]=shift_s(flat[*,i],-kk[i])

wy=20
w_glob=nx/10
w_loc=nx/50


;if(instrument eq 'GMOS' and n_fib eq 1500 and $
;        (maskname eq 'IFU' or maskname eq 'IFU-R')) then begin
    flat0=flat
    flat_new=flat
    fringes=flat*0.0
    for j=0,1 do begin
        flat_prof1 = median(flat[*,0:149],dim=2)
        bad_fp1=where(finite(flat_prof1) ne 1,cbfp1,compl=good_fp1)
        if(cbfp1 gt 0) then begin
            ff1=findgen(Nx)
            t1=interpol(flat_prof1[good_fp1],ff1[good_fp1],ff1)
            flat_prof1=t1
        endif
        flat_prof2 = median(flat[*,150:*],dim=2)
        bad_fp2=where(finite(flat_prof2) ne 1,cbfp2,compl=good_fp2)
        if(cbfp2 gt 0) then begin
            ff2=findgen(Nx)
            t2=interpol(flat_prof2[good_fp2],ff2[good_fp2],ff2)
            flat_prof2=t2
        endif
        if(j eq 1) then begin
            lowess,x,flat_prof1,70,ff1
            lowess,x,flat_prof2,70,ff2
            flat_prof1=ff1
            flat_prof2=ff2
        endif
        print,'Creating fringing pattern, step #'+string(j+1,format='(i1)')+'...'
        print,'Progress: 0%',format='(a,$)'
        perc0=0
        for i=0,n_fib-1 do begin
            perc=100.0*(i+1)/double(n_fib)
            if(perc ge perc0+2) then begin
                perc0=fix(perc)
                if(20*(perc0/20) eq perc0 and perc0 lt 100) then $
                    print,string(perc0,format='(i2)')+'%',format='(a,$)' $
                else $
                    print,'.',format='(a,$)'
            endif
            ;if(instrument eq 'GMOS' and n_fib eq 1500) then begin
                flat_prof=(i lt 150)? flat_prof1 : flat_prof2
                frel=flat[*,i]/flat_prof
                gg=where(finite(frel) eq 1,cgg)
                if(j eq 0) then begin
                    kv=robust_poly_fit(x[gg],frel[gg],4)
                    frel=poly(x,kv)
                    flat_new[*,i]=flat[*,i]/frel
                endif else begin
                    ;fringes[*,i]=shift_s(frel-1.0,kk[i])
                    fringes[*,i]=frel-1.0
                endelse
            ;endif
        endfor
        print,'100%'
        flat=flat_new
    endfor
    mwrfits, fringes[*,0:149],'calibration/'+rerun+'/fringes_pattern.fits', /create
    mwrfits, fringes[*,150:*],'calibration/'+rerun+'/fringes_pattern.fits'

;endif



final:
end
