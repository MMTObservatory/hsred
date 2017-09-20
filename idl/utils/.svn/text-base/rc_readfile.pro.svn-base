;+
; NAME:
;   rc_readfile
;
; PURPOSE:
;   Will read ascii file into arrays
;
; CALLING SEQUENCE:
;   outarray = rc_readfile(file, nskip, columns=columns)
;
; INPUTS:
;   file   - ascii file to be read
;   nskip  - how many lines to skip in the begining
;
; REQUIRED KEYWORDS:;
;
; OPTIONAL INPUTS:
;  columns -  columns of data (assumes all float).  This doesn't work for 
;             strings
;
;
; OUTPUTS:
;  outarray  - array of file
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   readfits()
;   sxpar()
;   sxaddpar
;   writefits
;
; REVISION HISTORY:
;   May 2004 - Written ; Cool UofA
;-
;-----------------------------------------------------------------------------
FUNCTION rc_readfile, fname, nskip, COLUMNS=columns
   
   
   IF NOT KEYWORD_SET(COLUMNS) THEN BEGIN
      
      GET_LUN, FILE
      OPENR, FILE, FNAME
      JUNK = ''
      N = 0L
      IF NSKIP GT 0 THEN BEGIN
         FOR I = 0, NSKIP-1 DO BEGIN
            READF, FILE, JUNK
         ENDFOR
      ENDIF
      WHILE NOT(EOF(FILE)) DO BEGIN
         READF, FILE, JUNK
         N = N+1
      ENDWHILE
      CLOSE, FILE
      
      IF N GT 0 THEN BEGIN
         DATATEMP = STRARR(N)
         
         OPENR, FILE, FNAME
         IF NSKIP GT 0 THEN BEGIN
            FOR I = 0, NSKIP-1 DO BEGIN
               READF, FILE, JUNK
            ENDFOR
         ENDIF
         READF, FILE, DATATEMP
         CLOSE, FILE
         FREE_LUN, FILE
         
         RETURN, DATATEMP
      ENDIF
      
      IF KEYWORD_SET(COLUMNS) THEN BEGIN
         GET_LUN, FILE
         OPENR, FILE, FNAME
         JUNK = ''
         N = 0L
         IF NSKIP GT 0 THEN BEGIN
            FOR I = 0, NSKIP-1 DO BEGIN
               READF, FILE, JUNK
            ENDFOR
         ENDIF
         WHILE NOT(EOF(FILE)) DO BEGIN
            READF, FILE, JUNK
            N = N+1
         ENDWHILE
         
         CLOSE, FILE
         
         DATATEMP = FLTARR(COLUMNS, N)
         
         OPENR, FILE, FNAME
         IF NSKIP GT 0 THEN BEGIN
            FOR I = 0, NSKIP-1 DO BEGIN
               READF, FILE, JUNK
            ENDFOR
         ENDIF
         READF, FILE, DATATEMP
         CLOSE, FILE
         FREE_LUN, FILE
      ENDIF ELSE  DATATEMP = '' 
      
      RETURN, DATATEMP
   ENDIF
   
END
