;+
;   NAME:
;      xdqso_calculate_prob
;   PURPOSE:
;      calculate the extreme-deconvolution XDQSO QSO probability
;   USE:
;      out= xdqso_calculate_prob(in,/dereddened)
;   INPUT:
;      in - structure containing PSFFLUX, PSFFLUX_IVAR, EXTINCTION
;   KEYWORDS:
;      dereddened - psfflux, and psfflux_ivar is already dereddened
;      galex - GALEX fluxes are included in psfflux, psfflux_ivar, and
;              extinction; use them
;   OUTPUT:
;      out - structure containing pqso, ... (see catalog description)
;   HISTORY:
;      2010-04-30 - Written - Bovy (NYU)
;      2010-05-29 - Added Galex - Bovy
;      2010-09-26 - Re-written to be stand-alone - Bovy
;-
FUNCTION XDQSO_CALCULATE_PROB, in, dereddened=dereddened
ndata= n_elements(in.psfflux[0])
IF ~keyword_set(dereddened) THEN BEGIN
    flux= xdqso_sdss_deredden(in.psfflux,in.extinction)
    flux_ivar=xdqso_sdss_deredden_error(in.psfflux_ivar,in.extinction)
ENDIF ELSE BEGIN
    flux= in.psfflux
    flux_ivar= in.psfflux_ivar
ENDELSE

;;Read the differential number counts
path= getenv('XDQSODATA')
if strcmp(path,'') then dataDir= '../data/' else dataDir = '$XDQSODATA/'
lumfunc= 'HRH07'
dndi_qsobosszfile= xdqso_dndipath(2.2,3.5,lumfunc)
dndi_qsofile= xdqso_dndipath(3.5,6.,lumfunc)
dndi_qsolowzfile= xdqso_dndipath(0.3,2.2,lumfunc)
dndi_everythingfile= dataDir+'dNdi_everything_coadd_1.4.prt'

xdqso_read_dndi, dndi_qsobosszfile, i_qsobossz, dndi_qsobossz, /correct
xdqso_read_dndi, dndi_qsofile, i_qso, dndi_qso, /correct
xdqso_read_dndi, dndi_qsolowzfile, i_qsolowz, dndi_qsolowz, /correct
xdqso_read_dndi, dndi_everythingfile, i_everything, dndi_everything

xdstruct= {qsolowzlike:0D, qsohizlike:0D, qsomidzlike: 0D, $
           starlike:0D, qsolowznumber:0D, qsohiznumber:0D,qsomidznumber:0D, $
           starnumber:0D, pstar:0D, pqsolowz:0D, pqsomidz:0D, pqsohiz:0D, $
           pqso: 0D}
out= replicate(xdstruct,ndata)

;;Now calculate all of the factors in turn
out.qsomidzlike= xdqso_eval_colorprob(flux,flux_ivar,/qso,/midz)
out.qsomidznumber= xdqso_eval_iprob(flux[3,*],i_qsobossz,dndi_qsobossz)
out.qsohizlike= xdqso_eval_colorprob(flux,flux_ivar,/qso)
out.qsohiznumber= xdqso_eval_iprob(flux[3,*],i_qso,dndi_qso)
out.qsolowzlike= xdqso_eval_colorprob(flux,flux_ivar,/qso,/lowz)
out.qsolowznumber= xdqso_eval_iprob(flux[3,*],i_qsolowz,dndi_qsolowz)
out.starlike= xdqso_eval_colorprob(flux,flux_ivar)
out.starnumber= xdqso_eval_iprob(flux[3,*],i_everything,dndi_everything)

;;Calculate the probabilities
pstar= out.starlike*out.starnumber
nonzero= where(pstar NE 0.)
IF nonzero[0] NE -1 THEN pstar[nonzero]= pstar[nonzero]/$
  (out[nonzero].qsolowzlike*out[nonzero].qsolowznumber+$
   out[nonzero].qsohizlike*out[nonzero].qsohiznumber+$
   out[nonzero].qsomidzlike*out[nonzero].qsomidznumber+$
   out[nonzero].starlike*out[nonzero].starnumber)
out.pstar= pstar
pqsolowz= out.qsolowzlike*out.qsolowznumber
nonzero= where(pqsolowz NE 0.)
IF nonzero[0] NE -1 THEN pqsolowz[nonzero]= pqsolowz[nonzero]/$
  (out[nonzero].qsolowzlike*out[nonzero].qsolowznumber+$
   out[nonzero].qsohizlike*out[nonzero].qsohiznumber+$
   out[nonzero].qsomidzlike*out[nonzero].qsomidznumber+$
   out[nonzero].starlike*out[nonzero].starnumber)
out.pqsolowz= pqsolowz
pqsohiz= out.qsohizlike*out.qsohiznumber
nonzero= where(pqsohiz NE 0.)
IF nonzero[0] NE -1 THEN pqsohiz[nonzero]= pqsohiz[nonzero]/$
  (out[nonzero].qsolowzlike*out[nonzero].qsolowznumber+$
   out[nonzero].qsohizlike*out[nonzero].qsohiznumber+$
   out[nonzero].qsomidzlike*out[nonzero].qsomidznumber+$
   out[nonzero].starlike*out[nonzero].starnumber)
out.pqsohiz= pqsohiz
pqsomidz= out.qsomidzlike*out.qsomidznumber
nonzero= where(pqsomidz NE 0.)
IF nonzero[0] NE -1 THEN pqsomidz[nonzero]= pqsomidz[nonzero]/$
  (out[nonzero].qsolowzlike*out[nonzero].qsolowznumber+$
   out[nonzero].qsohizlike*out[nonzero].qsohiznumber+$
   out[nonzero].qsomidzlike*out[nonzero].qsomidznumber+$
   out[nonzero].starlike*out[nonzero].starnumber)
out.pqsomidz= pqsomidz
out.pqso= out.pqsolowz+out.pqsomidz+out.pqsohiz
RETURN, out
END
