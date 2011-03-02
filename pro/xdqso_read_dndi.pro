;+
;   NAME:
;      xdqso_read_dndi
;   PURPOSE:
;      read the dndi number count prior files
;   INPUT:
;      filename - name of the file that holds the prior
;   KEYWORDS:
;      correct - correct the nunber counts for SDSS incompleteness
;      attenuate - attenuate the number counts
;      everythingcorrect - de-correct the number counts for SDSS incompleteness
;   OUTPUT:
;      i, dndi - i-mag and dndi number count prior (for xdqso_eval_iprob)
;   HISTORY:
;      2010 - Written - Bovy
;-
PRO XDQSO_READ_DNDI, filename, i, dndi, correct=correct, attenuate=attenuate, $
                     everythingcorrect=everythingcorrect
_CORRECTMEAN= 21.9
_CORRECTWIDTH= 0.2
_ATTENUATECUT= 24.5D
_ATTENUATESLOPE= -0.5D
OPENR, lun, filename, /GET_LUN
;;Read the header (1 line)
hdr= ""
READF, lun, hdr
cmd= "wc -l "+filename
spawn, cmd, nlines
nlines= nlines[0]-1
i= dblarr(nlines)
dndi= dblarr(nlines)
FOR ii=0L, nlines-1 DO BEGIN
    READF, lun, itmp, dnditmp
    i[ii]= itmp
    IF keyword_set(correct) THEN dndi[ii]= dnditmp/(exp((itmp-_CORRECTMEAN)/_CORRECTWIDTH)+1D0) $
    ELSE dndi[ii]= dnditmp
    IF keyword_set(attenuate) AND i[ii] GE _ATTENUATECUT THEN dndi[ii]= dndi[ii]*exp(_ATTENUATESLOPE*(i[ii]-_ATTENUATECUT))
    IF keyword_set(everythingcorrect) AND i[ii] GT 20. THEN BEGIN
        eight= where(i EQ 19.)
        twenty= where(i EQ 20.)
        dndi[ii]= dndi[eight]+(dndi[eight]-dndi[twenty])/(i[eight]-i[twenty])*(i[ii]-i[eight])
        dndi[ii]= dndi[ii]/(exp((itmp-_CORRECTMEAN)/_CORRECTWIDTH)+1D0)
    ENDIF
ENDFOR
; ESS free logical file unit
free_lun, lun
END
