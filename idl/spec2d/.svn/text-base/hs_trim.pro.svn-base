;+
; NAME:
;   hs_trim
;
; PURPOSE:
;   trim an image
;
; CALLING SEQUENCE:
;  hs_trime, image, header=, trim=, immask=, 
;
; INPUTS:
;   image = image to be trimmed
;
; OPTIONAL KEYWORDS:
;   header = image header
;   trim   - trim section
;   immask - image mask
;
; OUTPUTS:
;   
; OPTIONAL OUTPUTS:
;  
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
; BUGS:
;
;;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA
;                adapted from code in Schlegals IDLSPEC2d
;                 package
;-
;------------------------------------------------------------------------------

;;This program carries our simple trimming of the image.  

PRO hs_trim, image, header=header, trim=trim, immask=immask
  
  imsize = size(image,/dimension)
;  if n_elements(trim) eq 0L then xy = [0,imsize[0]-1,0,imsize[1]-1] else $
;      xy = long(trim)
  
  if keyword_set(trim) then trim = long(trim) else begin
     
   trim = long(strsplit(sxpar(header, 'TRIMSEC'), $
                        '[:,]', /EXTRACT))-1l
   xy=trim


   endelse

    xy[0] = xy[0] > 0L
    xy[1] = xy[1] < (imsize[0]-1L)
    xy[2] = xy[2] > 0L
    xy[3] = xy[3] < (imsize[1]-1L)
    
    ncols = xy[1]-xy[0]+1L
    nrows = xy[3]-xy[2]+1L

    origimage = image ; original data
    if keyword_set(immask) then origmask = immask
    
    image = origimage[xy[0]:xy[1],xy[2]:xy[3]] ; trim
    if keyword_set(immask) then immask = origmask[xy[0]: $ 
                                xy[1],xy[2]:xy[3]] ; trim mask

    
  if keyword_set(header) then begin

       sxdelpar, header, 'DATASEC'
       sxdelpar, header, 'TRIMSEC'

       sxaddpar, header, 'NAXIS1', ncols
       sxaddpar, header, 'NAXIS2', nrows

       sxaddpar, header, 'TRIM', im_today()+' Trim section '+$
         '['+strn(xy[0]+1)+':'+strn(xy[1]+1)+','+strn(xy[2]+1)+$
         ':'+strn(xy[3]+1)+'].', before='HISTORY'

       xy = xy+1L
       sxaddpar, header, 'CCDSEC', '['+strn(xy[0])+':'+strn(xy[1])+','+$
         strn(xy[2])+':'+strn(xy[3])+']'

    endif

return
end
 
