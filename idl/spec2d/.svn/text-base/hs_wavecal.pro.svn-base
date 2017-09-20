;+
; NAME:
;   hs_wavecal
;
; PURPOSE:
;   Determine wavelength calibration from arclines
;
; CALLING SEQUENCE:
;   hs_wavecal, arc, arcivar, xcen, ycen, wset, [wfirst=, $
;    color=color, lampfile=lampfile, fibermask=fibermask, $
;    func=func, aset=aset, ncoeff=ncoeff, lambda=lambda, $
;    thresh=thresh, row=row, nmed=nmed, /gauss, $
;    xdif_tset=xdif_tset, bestcorr=bestcorr , ccdnum=, arcguess=arcguess
;
; INPUTS:
;   arc        - Extracted arc spectra with dimensions [NY,NFIBER]
;   arcivar    - Inverse variance of ARC
;
; OPTIONAL KEYWORDS:
;   color      - 'red' or 'blue'; not required if ANS is set
;   lampfile   - Name of file describing arc lamp lines;
;                default to the file 'lamphgcdne.dat' in $IDLSPEC2D_DIR/etc.
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
;   func       - Name of fitting function; default to 'legendre'
;   aset       - Trace set for initial wavelength solution in row number ROW.
;   ncoeff     - Number of coefficients in fits.  This may be different than
;                the number of coefficients in the initial guess ASET.
;                Default to 5.
;   thresh     - Threshhold counts for significant lines;
;                default to 200 if COLOR='blue' or 500 if COLOR='red'
;   row        - Row to use in initial guess of wavelength solution;
;                default to (NFIBER-30)/2
;   nmed       - Number of rows around ROW to median filter for initial
;                wavelengths solution; default to 5
;   maxdev     - max deviation in log lambda to allow (default 1.0e-5=7 km/s)
;   gauss      - Use gaussian profile fitting for final centroid fit
;
;   parity     - Even or odd?
;
;   ccdnum     - CCD number
;
; OUTPUTS:
;   aset       - (Modified)
;   xcen       - pixel position of lines [nfiber, nlambda]
;   ycen       - fiber number [nfiber, nlambda]
;   wset       - traceset (pix -> lambda)
;
; OPTIONAL OUTPUTS:
;   lampfile   - Modified from input to include full path name of file
;   lambda     - Returns wavelengths of good lamp lines [Angstroms]
;   fibermask  - (Modified)
;   xdif_tset  - Fit residual of lamp lines to fit positions [pixels]
;   bestcorr   - Correlation coefficient with simulated arc spectrum
;   wfirst     - traceset from first iteration on arc fits
;
; COMMENTS:
;   Return from routine after computing BESTCORR if XCEN, YCEN and WSET
;   are not to be returned.
;
; EXAMPLES:
;
; BUGS:
;
;
; PROCEDURES CALLED:
;   arcfit_guess()
;   djs_median
;   djsig()
;   fibermask_bits()
;   trace_crude()
;   trace_fweight()
;   traceset2pix()
;   traceset2xy()
;   xy2traceset
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/lamphgcdne.dat
;
; REVISION HISTORY:
;   15-Oct-1999  Written by S. Burles, D. Finkbeiner, & D. Schlegel, APO.
;   09-Nov-1999  Major modifications by D. Schlegel, Ringberg.
;   20-Jan-2000  Gone back to very simple procedure: replacement (S. Burles)
;   April 2004 - Modified for HS data - 
;                R Cool U of A
;   April 2014 - Enabled and tested 600-line reductions for hectospec,
;                refined 270-line reductions, fixed bug where dcoef
;                was altered in HS_IARCFIT, and then returned back to 
;                here, to be used again (incorrectly)
;-
;------------------------------------------------------------------------------
pro hs_wavecal, arc, arcivar, xcen, ycen, wset, wfirst=wfirst, $
                lampfile=lampfile, fibermask=fibermask, $
                func=func, aset=aset, ncoeff=ncoeff, lambda=lambda, $
                thresh=thresh, $
                row=row, nmed=nmed, xdif_tset=xdif_tset, bestcorr=bestcorr, $
                gauss=gauss, maxdev=maxdev, parity=parity,lamptype=lamptype, $
                ccdnum=ccdnum, arcguess=arcguess, doplot=doplot, $
                chelle=chelle, guessfile=guessfile, center600=center600
   
  if not keyword_set(lamptype) then  lamptype = 'henear'
  
  if (NOT keyword_set(ccdnum)) then $
     splog, "You must specify the ccdnumber"
  if (NOT keyword_set(func)) then func = 'legendre'
  if (NOT keyword_set(ans)) then ans = 0
  if NOT keyword_set(maxdev) then maxdev = 1.0d-5
  
  
  IF keyword_set(guessfile) THEN BEGIN
     guess = mrdfits(guessfile, ccdnum)
     arcguess = guess.coeff
  ENDIF
    
  if NOT keyword_set(arcguess) then begin
     arcguess = fltarr(4,4)
     arcguess[*,0] =  [6390.0557, 2783.7124  ,20.804438  ,-24.199860]
     arcguess[*,1] =  [6431.2490, 2783.2424  ,20.712286,  -24.158333]
     arcguess[*,2] =  [6406.5127, 2812.5986,  20.695004,  -24.229319]
     arcguess[*,3] =  [6445.8882, 2812.6067,  20.457697,  -24.949741]
     if keyword_set(center600) then begin
         case center600 of
            4800: begin
                      arcguess[*,0]=[4778.7974, 1270.4028, 6.6831183, -10.866278]
                      arcguess[*,1]=[4795.7915,1267.7461,7.2083111,-10.929058]
                      arcguess[*,2] =[4790.8833, 1283.6621, 6.8641562, -11.046491]
                      arcguess[*,3] =[4808.3394, 1283.7007, 6.8449359, -11.034074]
                      end
            5300: begin
                      arcguess[*,0]=[5278.7974, 1270.4028, 6.6831183, -10.866278]
                      arcguess[*,1]=[5296.5698, 1270.5659, 6.7214098, -10.880716]
                      arcguess[*,2] =[5290.8833, 1283.6621, 6.8641562, -11.046491]
                      arcguess[*,3] =[5308.3394, 1283.7007, 6.8449359, -11.034074]
                      end
            5800: begin
                      arcguess[*,0] =  [5782.8042,1253.9785,6.0308204,-10.435096]
                      arcguess[*,1] =  [5782.0503,1253.9689,6.7895374,-10.411214]
                      arcguess[*,2] =  [5800.2446,1271.7238 ,5.9599128,-10.674212]
                      arcguess[*,3] =  [5799.4746,1271.7133,6.7269101,-10.670746]
                      end
            6300: begin
                      arcguess[*,0]=[6282.3188,1275.9530,5.8132463,-11.233899]
                      arcguess[*,1]=[6299.6802,1276.1913,5.7900157,-11.192946]
                      arcguess[*,2]=[6295.2065,1289.6493,5.8571243,-11.318105]
                      arcguess[*,3]=[6311.4155,1289.8351,5.6694608,-10.989108]
                      end
            6600: begin
                      arcguess[*,0]=[6582.3188,1275.9530,5.8132463,-11.233899]
                      arcguess[*,1]=[6599.6802,1276.1913,5.7900157,-11.192946]
                      arcguess[*,2]=[6595.2065,1289.6493,5.8571243,-11.318105]
                      arcguess[*,3]=[6611.4155,1289.8351,5.6694608,-10.989108]
                      end
            6800: begin
                      arcguess[*,0] =  [6788.9418, 1290.3749,3.6254,-6.1661]
                      arcguess[*,1] =  [6811.37616,1291.03918,3.5781509,-6.49563]
                      arcguess[*,2] =  [6811.4636,1290.9379,3.4650,-6.3766]
                      arcguess[*,3] =  [6810.9045,1290.96257,4.0164,-6.3843]
                      end
            7300: begin
                      arcguess[*,0]=[7269.2183,1278.7847,4.8519149,-11.052988]
                      arcguess[*,1]=[7286.4189,1278.9600,4.8373795,-11.018242]
                      arcguess[*,2]=[7283.2954,1292.4255,4.9463630,-11.149488]
                      arcguess[*,3]=[7299.3516,1292.4120,4.9743080,-11.166810]
                      end
            7800: begin
                      arcguess[*,0]=[7772.7485,1280.1208,4.2924218,-11.025270]
                      arcguess[*,1]=[7789.7837,1280.3291,4.3078752,-11.023109]
                      arcguess[*,2]=[7787.5361,1293.8228,4.4102864,-11.182652]
                      arcguess[*,3]=[7803.4487,1293.7925,4.4336557,-11.171949]
                      end
            ELSE: begin
                      arcguess[*,0] =  [5283.0146,1254.0004,6.2916389,-10.410913]
                      arcguess[*,1] =  [5327.5391,1209.8246,4.8536382,-9.3512888]
                      arcguess[*,2] =  [5282.4985,1253.8240,6.7775106,-10.168923]
                      arcguess[*,3] =  [5327.0796,1209.6476,5.3118129,-9.1558685]
                      arcguess[0,*]=arcguess[0,*]+center600-5300.
            END
         endcase
     endif
  endif
  t_begin = systime(1)
  
  ndim = size(arc, /n_dim)
  if (ndim NE 2) then $
     message, 'Expecting 2-D arc image'
  dims = size(arc, /dim)
  npix = dims[0]
  nfiber = dims[1]
  if (total(dims NE size(arcivar, /dim))) then $
     message, 'ARC and ARCIVAR must have same dimensions'
  
  if (NOT keyword_set(fibermask)) then fibermask = bytarr(nfiber)
  if (NOT keyword_set(row)) then row = (nfiber-30)/2
  if (NOT keyword_set(nmed)) then nmed = 5
  if (NOT keyword_set(thresh)) then begin
     thresh=000
  endif
  if (NOT keyword_set(ncoeff)) then ncoeff = 5
    
  ;;---------------------------------------------------------------------------
  ;; Read LAMPLIST file for wavelength calibration
  ;;---------------------------------------------------------------------------
  ;; Read this file into a structure
  
  CASE strupcase(lamptype) OF                                              
     'HEAR': lampfile = 'lamphear.dat'                                        
     'HENEAR': if keyword_set(center600) then lampfile = 'lamphenear_'+strtrim(string(fix(center600), format='(I4)'),2)+'.dat' else lampfile = 'lamphenear.dat'                                    
     'THAR' : lampfile = 'lampthar.dat'                                       
     'PENRAY' : lampfile = 'lamppenray_'+strtrim(string(fix(center600), format='(I4)'),2)+'.dat'                                       
     ELSE: BEGIN                                                              
        splog, 'Unsupported comparison arc lamp.'                             
        return                                                                
     endelse                                                                  
  endcaSE       
  
  
  if (keyword_set(lampfile)) then begin
     lampfilename = (findfile(lampfile, count=ct))[0]
     if (ct EQ 0) then BEGIN
        lampfilename = filepath(lampfile, $
                                root_dir=getenv('HSRED_DIR'), $
                                subdirectory='etc')
        lampfilename = (findfile(lampfilename, count=ct))[0]
        IF ct EQ 0 THEN $
           message, 'No LAMPFILE found '+lampfile
     ENDIF
        
  endif else begin
     lampdefault = filepath('lamphenear.dat', $
                            root_dir=getenv("HSRED_DIR"), subdirectory='etc')
     lampfilename = (findfile(lampdefault, count=ct))[0]
     if (NOT keyword_set(lampfilename)) then $
        message, 'No LAMPFILE found '+lampdefault
  endelse
  
  splog, 'Reading lamp file ', lampfilename
  readcol, lampfilename, lampwave, lampinten, lampquality, format='D,F,A'
  lamps = {lambda: 0.0d0, loglam: 0.0d0, intensity: 0.0d0, good: 0.0d0}
  lamps = replicate(lamps, N_elements(lampwave))
  lamps.lambda = lampwave
  lamps.loglam = alog10(lampwave)
  lamps.intensity = lampinten
  lamps.good = strupcase(lampquality) EQ 'GOOD' AND lampinten GT 0

  ;;---------------------------------------------------------------------------
  ;; INITIAL WAVELENGTH SOLUTION
  ;;---------------------------------------------------------------------------
  
  ;; Find the initial wavelength solution if ANS is not passed
  ;; One might want to change nave and nmed for first pass???
     
  if (NOT keyword_set(aset)) then begin
     ;; Extract one spectrum from the NMED spectra around fiber number ROW
     ;; by taking the median value at each wavelength.
     ;; Find the NMED fibers nearest to ROW that are not masked.
     ii = where(fibermask EQ 0, ngfiber)
     if (ngfiber EQ 0) then $
        message, 'No unmasked fibers according to FIBERMASK'
     ii = ii[ sort(abs(ii-row)) ]
     ii = ii[0:(nmed<ngfiber)]  ; At most NGFIBER good fibers
     
     spec = djs_median(arc[*,ii], 2)
     
     dimens = size(arc, /dimen)
     nrows = dimens[0]
     nspec = dimens[1]
     
     aset ={func:'' , xmin:0.0d0, xmax:0.0d0, coeff:dblarr(6,nspec)}
     coeff=dblarr(6,nspec)
     
     for i = 0, 20 -1 do begin
        im = dblarr(nrows, 10)
        imivar = im
        
        for j = 0, 9 do begin
           im[*,j] = arc[*,i]
           imivar[*,j] = arcivar[*,i]/sqrt(10)
        endfor
       
        
        xmid = 0.5*4607
        if ccdnum eq 1 then begin
           if i eq 0 then begin
              acoeff = arcguess[*,0]
           endif else if i/2 eq i/2.  then begin
              acoeff = (asetguesseven.coeff)[0:3]
              acoeff[0] = acoeff[0]
           endif else if i/2 ne i/2. then begin
              if i eq 1 then asetguessodd = asetguesseven
              acoeff = (asetguessodd.coeff)[0:3]
              acoeff[0] = acoeff[0]
              if i eq 1 then acoeff = arcguess[*,1]
           endif
           
        endif else if ccdnum eq 2 then begin
                      
           if i eq 0 then begin
              acoeff = arcguess[*,2]
           endif else  if  i/2 eq i/2.  then begin
              acoeff = (asetguesseven.coeff)[0:3]
              acoeff[0] = acoeff[0]
           endif else if i/2 ne i/2. then begin
              if i eq 1 then asetguessodd = asetguesseven
              acoeff = (asetguessodd.coeff)(0:3)
              acoeff[0] = acoeff[0]
              if i eq 1 then acoeff = arcguess[*,3]
           endif
        endif
        
        IF keyword_set(guessfile) THEN acoeff = arcguess(0:3,i)
   
        ;need to reset dcoef each time, because it gets modified in the calls below
        dcoeff = [100,50,5,5]
        if keyword_set(center600) then  begin
             dcoeff = [50,25,2.5,2.5]
             if center600 eq 4800 then  dcoeff = [2.5,1.25,0.1,0.1]
        endif
        IF keyword_set(chelle) THEN dcoeff=[10, 5, 1, 1]

        hs_iarcfit, im, imivar,  wset, lamptype=lamptype, $
           acoeff=acoeff, dcoeff=dcoeff, xdiff=xdiff1, $
           arclambda = arclambda1, norder_arc = 5, $
           mintol=0.5, doplot=doplot, chelle=chelle, center600=center600
        
        if i/2 eq i/2. then asetguesseven = wset
        if i/2 ne i/2. then asetguessodd = wset
        if i eq 0 then begin
           aset.func = wset.func
           aset.xmin = wset.xmin
           aset.xmax = wset.xmax
        endif   
        coeff[*,i] = (wset.coeff)[*,5]
     endfor
     
     
     ;; Now reverse
     splog, 'STARTING REVERSAL'
     
     
     for kl = 0, 19 do begin
        i = 19-kl
        im = dblarr(nrows, 10)
        imivar = im
        for j = 0, 9 do begin
           im[*,j] = arc[*,i]
           imivar[*,j] = arcivar[*,i]/sqrt(10.)
        endfor
                
        if ccdnum eq 1 then begin
           
           if  i/2 eq i/2.  then begin
              acoeff = (asetguesseven.coeff)[0:3]
              acoeff[0] = acoeff[0]
           endif else if i/2 ne i/2. then begin
              acoeff = (asetguessodd.coeff)[0:3]
              acoeff[0] = acoeff[0]
           endif
           
        endif else if ccdnum eq 2 then begin
           
           if  i/2 eq i/2.  then begin
              acoeff = (asetguesseven.coeff)[0:3]
              acoeff[0] = acoeff[0]
           endif else if i/2 ne i/2. then begin
              acoeff = (asetguessodd.coeff)[0:3]
              acoeff[0] = acoeff[0]
           endif
        endif
        
        IF keyword_set(guessfile) THEN acoeff = arcguess(0:3,i)

        ;need to reset dcoef each time, because it gets modified in the calls below
        dcoeff = [100,50,5,5]
        if keyword_set(center600) then  begin
           dcoeff = [50,25,2.5,2.5]
           if center600 eq 4800 then  dcoeff = [2.5,1.25,0.1,0.1]
        endif
        IF keyword_set(chelle) THEN dcoeff=[10, 5, 1, 1]

        hs_iarcfit, im, imivar,  wset, lamptype=lamptype, $
           acoeff=acoeff, dcoeff=dcoeff, xdiff=xdiff1, $
           arclambda = arclambda1,  norder_arc=5, $
           mintol=0.25, doplot=doplot, chelle=chelle, center600=center600
        
        if i /2 eq i/2. then asetguesseven = wset
        if i/2 ne i/2. then asetguessodd = wset
        if i eq 0 then begin
           aset.func = wset.func
           aset.xmin =  wset.xmin
           aset.xmax = wset.xmax
        endif   
        
        coeff[*,i] = (wset.coeff)(*,5)
     endfor 
     

     for i = 0, nspec-1 do begin
        
        im = dblarr(nrows, 10)
        imivar = im
        for j = 0, 9 do begin
           im(*,j) = arc(*,i)
           imivar(*,j) = arcivar(*,i)/sqrt(10.)
        endfor
        
        
        if ccdnum eq 1 then begin
           
           
           if  i/2 eq i/2.  then begin
              acoeff = (asetguesseven.coeff)[0:3]
              acoeff[0] = acoeff[0]
           endif else if i/2 ne i/2. then begin
              acoeff = (asetguessodd.coeff)[0:3]
              acoeff[0] = acoeff[0]
           endif
           
        endif else if ccdnum eq 2 then begin
           
           if  i/2 eq i/2.  then begin
              acoeff = (asetguesseven.coeff)[0:3]
              acoeff[0] = acoeff[0]
           endif else if i/2 ne i/2. then begin
              acoeff = (asetguessodd.coeff)[0:3]
              acoeff[0] = acoeff[0]              
           endif
        endif
        
        IF keyword_set(guessfile) THEN acoeff = arcguess(0:3,i)

        ;need to reset dcoef each time, because it gets modified in the calls below
        dcoeff = [100,50,5,5]
        if keyword_set(center600) then  begin
           dcoeff = [50,25,2.5,2.5]
           if center600 eq 4800 then  dcoeff = [2.5,1.25,0.1,0.1]
        endif
        IF keyword_set(chelle) THEN dcoeff=[10, 5, 1, 1]

        hs_iarcfit, im, imivar,  wset, lamptype=lamptype, $
           acoeff=acoeff, dcoeff=dcoeff, xdiff=xdiff1, $
           arclambda = arclambda1,  norder_arc=5, $
           mintol=0.25, doplot=doplot, chelle=chelle, center600=center600
        
        if i /2 eq i/2. then asetguesseven = wset
        if i/2 ne i/2. then asetguessodd = wset
        if i eq 0 then begin
           aset.func = wset.func
           aset.xmin =  wset.xmin
           aset.xmax = wset.xmax
        endif   
        
        coeff[*,i] = (wset.coeff)(*,5)
     endfor 
     
     aset.coeff = coeff   
     wset = aset
  endif 

  return                                                        
  
end
;------------------------------------------------------------------------------
