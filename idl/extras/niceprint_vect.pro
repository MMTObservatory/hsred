pro niceprint_vect, input

    on_error, 2
    
    nloop = n_elements(input(*,0))
    nvar = n_elements(input(0,*))
    if nvar eq 0 then nvar =1 
    if nvar eq 1 then nloop = n_elements(input)
    
    
 
    
    for i=0,nloop-1 do begin
       
       output = ''
       if nvar gt 1 then $
          for j = 0, nvar-1 do output = output + $
          string(input(i,j), format='(A20)') + ' ' 

       if nvar eq 1 then $
          output = string(input(i))
       print, output
    
    endfor
 

    return
    
end
