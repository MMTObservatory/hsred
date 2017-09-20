;;This program will read the plugplan and use it to set up the proper directory stucture

PRO hs_makebootesstruct, rerun=rerun, basepath=basepath

   spawn, 'mkdir ' + basepath + rerun
   spawn, 'mkdir ' + basepath + rerun + '/field01'
   spawn, 'mkdir ' + basepath + rerun + '/field02'
   spawn, 'mkdir ' + basepath + rerun + '/field03'
   spawn, 'mkdir ' + basepath + rerun + '/field04'
   spawn, 'mkdir ' + basepath + rerun + '/field05'
   spawn, 'mkdir ' + basepath + rerun + '/field06'
   spawn, 'mkdir ' + basepath + rerun + '/field07'
   spawn, 'mkdir ' + basepath + rerun + '/field08'
   spawn, 'mkdir ' + basepath + rerun + '/field09'
   spawn, 'mkdir ' + basepath + rerun + '/field10'
   spawn, 'mkdir ' + basepath + rerun + '/field11'
   spawn, 'mkdir ' + basepath + rerun + '/field12'
   spawn, 'mkdir ' + basepath + rerun + '/field13'
   spawn, 'mkdir ' + basepath + rerun + '/field14'
   spawn, 'mkdir ' + basepath + rerun + '/field15'
   
   for i = 1, 3 do begin
      for j = 1, 15 do begin
         spawn, 'mkdir ' + basepath + rerun + '/field' + $
            string(j, format='(i2.2') + '/pass' $
            + string(i, format='(i1.1)')
      endfor
   endfor
   
END
