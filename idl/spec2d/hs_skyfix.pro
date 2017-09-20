;;This function will create an array of scales that will correct the 
;;overall scale differences accross the focal plane of HECTOSPEC

FUNCTION hs_skyfix, plugmap
  
  rad = sqrt(plugmap.xfocal^2+plugmap.yfocal^2)
  ratio = 1.11748 + 0.0001866*rad -1.88622e-6*rad^2
                                ; this was from dohecto.cl
  
  return, ratio
end

