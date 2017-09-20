;+
; NAME:
;       ICLEANUP
;
; PURPOSE:
;       Free RD2DSPEC() structure pointer and array memory. 
;
; CALLING SEQUENCE:
;       icleanup, cube, main=main
;
; INPUTS:
;       cube - data cube from RD2DSPEC()
;
; KEYWORD PARAMETERS:
;       main - call RETALL
;
; OUTPUTS:
;       None.
;
; PROCEDURES USED:
;
; EXAMPLE:
;
; OLD MODIFICATION HISTORY:
;       J. Moustakas, 2001 August 23, U of A, written
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 April 9, U of A, checked out ISPEC v1.0.0
;       jm03dec8uofa - added SKY structure entry
;-

pro icleanup, cube, main=main

    ncube = size(cube,/n_elements)
    if ncube eq 0L then return
    
    for i = 0L, ncube-1L do begin

       if (ptr_valid(cube[i].header))[0] then ptr_free, cube[i].header

       if tag_exist(cube[i],'image') then cube[i].image = 0.0
       if tag_exist(cube[i],'sky') then cube[i].sky = 0.0
       if tag_exist(cube[i],'vmap') then cube[i].vmap = 0.0
       if tag_exist(cube[i],'mask') then cube[i].mask = 0.0
       if tag_exist(cube[i],'spec') then cube[i].spec = 0.0
       if tag_exist(cube[i],'sigspec') then cube[i].sigspec = 0.0
       if tag_exist(cube[i],'wave') then cube[i].wave = 0.0

    endfor
    
    if keyword_set(main) then retall

return    
end

