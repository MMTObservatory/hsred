FUNCTION strreplace, instring1, matchcode, replacecode
   
   outstring1 = strarr(n_elements(instring1))
   
   FOR io = 0, n_elements(instring1) -1 DO BEGIN
      instring = instring1(io)
      
   
      ;;Split instring into an N element array
      
      junk = strsplit(instring, ' ', length=ct)
      inarray = strarr(ct)
      FOR i = 0, ct(0) -1 DO BEGIN
         inarray(i) = strmid(instring, i, 1)
      ENDFOR
   
      ;;Do the same thing for the matchcode
      
      junk = strsplit(matchcode, ' ', length=ct)
      matcharray = strarr(ct)
      FOR i = 0, ct(0) -1 DO BEGIN
         matcharray(i) = strmid(matchcode, i, 1)
      ENDFOR
      
      ;;Do the same thing for the replacecode
      
      junk = strsplit(replacecode, ' ', length=ct)
      replacearray = strarr(ct)
      FOR i = 0, ct(0) -1 DO BEGIN
         replacearray(i) = strmid(replacecode, i, 1)
      ENDFOR
      
      
      IF n_elements(matcharray) GT n_elements(inarray) THEN $
         message, 'Your match string should be shorter than your input'
      
      rct = n_elements(replcatearray)
      mct = n_elements(matcharray)
      diff = rct-mct
      delvarx, index
      
      FOR i = 0, n_elements(inarray)-n_elements(matcharray) DO BEGIN
         IF inarray(i) EQ matcharray(0) THEN BEGIN
            IF total(inarray(i:i+mct-1) EQ matcharray) EQ mct THEN BEGIN
               IF n_elements(index) EQ 0 THEN index = i ELSE $
                  index = [index, i]
            ENDIF
         ENDIF
      ENDFOR
      
      
      char = ['!','@','#','$','%','^','&','*','(',')','{','}','[',']','|']
      FOR i = 0, n_elements(char) -1 DO BEGIN
         k = where(inarray EQ char(i), ct1)
         k = where(matcharray EQ char(i), ct2)
         k = where(replacearray EQ char(i), ct3)
         ct = ct1+ct2+ct3
         IF ct EQ 0 THEN charindex = i
      endfor
      
      char = char(charindex)
      outarray = inarray
      
      IF n_elements(index) GT 0 THEN BEGIN
         FOR i = 0, n_elements(index) -1 DO BEGIN
            outarray(index(i):index(i)+mct-1) = replicate(char, mct)
         endfor
      endif
      
      outstring = outarray(0)
      FOR i = 1, n_elements(outarray) -1 DO BEGIN
         outstring = outstring + outarray(i)
      ENDFOR
      
      split = strsplit(outstring, char, /extract)
      out = split(0)
      FOR i = 1, n_elements(split) -1 DO BEGIN
         out = out + replacecode + split(i)
      endfor
      
      IF strmid(out, 0, 1) NE strmid(outstring, 0, 1) THEN $
         out = replacecode + out
      junk = strsplit(outstring, ' ', length=ct)
      junk = strsplit(out, ' ', length=ct1)
      
      IF strmid(out, ct1-1, 1) NE strmid(outstring, ct-1, 1) THEN $
         out = out + replacecode
      
      outstring1(io) = out
   ENDFOR
   
   
   return, outstring1
      
      
END
