PRO hs_plugout, infile, outfile=outfile
  
  ;;Here infile is the spHect file that you want the information for
  ;;outfile is the text file that will go to default -> same root as 
  ;;infile
  
  
  if not keyword_set(outfile) then begin
     
     junk = strsplit(infile, ' ', length=ct)
     junk = strmid(infile, 0, ct-5)
     outfile = junk + '.dat'
     
  endif
  
  ;;now read
  
  plugmap = mrdfits(infile, 5)
  
  output = string(plugmap.fiberid, format='(i3)') + '  ' + $
     string(plugmap.objtype, format='(A10)') + '  ' + $
     string(plugmap.ra, format='(D10.6)') + '  ' + $
     string(plugmap.dec, format='(D10.6)') + '  ' +$
     string(plugmap.rmag, format='(D8.5)')
  
  
  header = '# Fiberid    Objecttype    RA   DEC   RMAG'
  
  openw,1,outfile
  printf, 1, header
  for i = 0, n_elements(output)-1 do printf, 1, output(I)
  close, 1
  
END


