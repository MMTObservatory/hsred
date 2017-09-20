function morph2type, morph, help=help
; jm02may21uofa
; convert a morphology (string) to RC3 numerical type

    if keyword_set(help) then begin

       print, 'The available morphological types are:'
       
       print, ['cE','E0','E0-1','E+','S0-','S0',$
                   'S0+','S0/a','Sa','Sab','Sb','Sbc',$
                   'Sc','Scd','Sd','Sdm','Sm','I0',$
                   'Im','cI','Pec'] 
       return, -99L

    endif
    
    nmorph = n_elements(morph)
    if (nmorph eq 0L) or (size(morph,/type) ne 7L) then begin
       print, 'Syntax - type = morph2type(morph,help=help)'
       return, -99L
    endif

    type = lonarr(nmorph)
    for i = 0L, nmorph-1L do begin
    
       case strlowcase(strcompress(morph[i],/remove)) of
          'ce'  : type[i] = -6
          'e0'  : type[i] = -5
          'e0-1': type[i] = -5
          'e+'  : type[i] = -4
          's0-' : type[i] = -3
          's0'  : type[i] = -2
          's0+' : type[i] = -1
          's0/a': type[i] = 0
          'sa'  : type[i] = 1
          'sab' : type[i] = 2
          'sb'  : type[i] = 3
          'sbc' : type[i] = 4
          'sc'  : type[i] = 5
          'scd' : type[i] = 6
          'sd'  : type[i] = 7
          'sdm' : type[i] = 8
          'sm'  : type[i] = 9
          'i0'  : type[i] = 90
          'im'  : type[i] = 10
          'ci'  : type[i] = 11
          'pec' : type[i] = 99
          else: begin
             print, 'Unrecognized type '+morph[i]+'...classifying as Pec.'
             type[i] = 99
          end
       endcase

    endfor

return, type
end    
