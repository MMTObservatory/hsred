PRO hs_convolve_models, newres=newres
  
  if not keyword_set(newres) then newres = 6.0
  sdssres = 2.78
  
  ;;Create the needed kernal to get new res from sdssres
  fwhm = sqrt( newres^2 - sdssres^2)
  
  data = gauss1(findgen(100)-50,[0,2.2,1])
  psf = gauss1(findgen(20)-10, [0,1.0, 1])
  

  
  test = convolve(data,psf)
  
  
  
  plot, findgen(100)-50, data
  djs_oplot, findgen(20)-10, psf,color='red'
  djs_oplot, findgen(100)-50, test,color='green'
  
  
  
  
END

