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
;      2014-04-02 - Added WISE, galex, ukidss - DiPompeo (UWyo)
;-
FUNCTION XDQSO_CALCULATE_PROB, in, dereddened=dereddened, galex=galex, ukidss=ukidss, wise=wise
ndata= n_elements(in.psfflux[0])

flux=in.psfflux
flux_ivar=in.psfflux_ivar
IF tag_exist(in,'extinction') THEN extinction=in.extinction ELSE extinction=dblarr(n_elements(flux[*,0]),n_elements(in))

;;MAD set keywords to pass to xdqso_eval_colorprob.pro
IF (keyword_set(galex)) THEN BEGIN
 galex=1
 uvflux=dblarr(2.,n_elements(in))
 uvflux_ivar=dblarr(2.,n_elements(in))
 uvextinction=dblarr(2.,n_elements(in))
 uvflux[0,*]=in.nuv
 uvflux[1,*]=in.fuv
 uvflux_ivar[0,*]=in.nuv_ivar
 uvflux_ivar[1,*]=in.fuv_ivar
 IF ~keyword_set(dereddened) THEN BEGIN
    uvextinction[0,*]= in.extinction[0]/5.155D0*8.18D
    uvextinction[1,*]= in.extinction[0]/5.155D0*8.29D
 ENDIF
ENDIF ELSE BEGIN
 galex=0
 uvflux=0.
 uvflux_ivar=0.
 uvextinction=0.
ENDELSE
IF (keyword_set(ukidss)) THEN BEGIN
 ukidss=1 
 nirflux=dblarr(4.,n_elements(in))
 nirflux_ivar=dblarr(4.,n_elements(in))
 nirextinction=dblarr(4.,n_elements(in))
 _SITONMGY= 1D35/3631
 nirflux[0,*]=in.APERCSIFLUX3_Y*_SITONMGY
 nirflux[1,*]=in.APERCSIFLUX3_J*_SITONMGY
 nirflux[2,*]=in.APERCSIFLUX3_H*_SITONMGY
 nirflux[3,*]=in.APERCSIFLUX3_K*_SITONMGY
 nirflux_ivar[0,*]=1D0/(in.APERCSIFLUX3ERR_Y*_SITONMGY)^2D0
 nirflux_ivar[1,*]=1D0/(in.APERCSIFLUX3ERR_J*_SITONMGY)^2D0
 nirflux_ivar[2,*]=1D0/(in.APERCSIFLUX3ERR_H*_SITONMGY)^2D0
 nirflux_ivar[3,*]=1D0/(in.APERCSIFLUX3ERR_K*_SITONMGY)^2D0
 IF ~keyword_set(dereddened) THEN BEGIN
    nirextinction[0,*]= in.extinction[0]/5.155D0*1.259D0
    nirextinction[1,*]= in.extinction[0]/5.155D0*0.920D0
    nirextinction[2,*]= in.extinction[0]/5.155D0*0.597D0
    nirextinction[3,*]= in.extinction[0]/5.155D0*0.369D0
 ENDIF
ENDIF ELSE BEGIN
 ukidss=0
 nirflux=0.
 nirflux_ivar=0.
 nirextinction=0.
ENDELSE
IF (keyword_set(wise)) THEN BEGIN 
 wise=1 
 wiseflux=dblarr(2.,n_elements(in))
 wiseflux_ivar=dblarr(2.,n_elements(in))
 wiseextinction=dblarr(2.,n_elements(in))
 wiseflux[0,*]=in.w1_nanomaggies
 wiseflux[1,*]=in.w2_nanomaggies
 wiseflux_ivar[0,*]=in.w1_nanomaggies_ivar
 wiseflux_ivar[1,*]=in.w2_nanomaggies_ivar
 IF ~keyword_set(dereddened) THEN BEGIN
    wiseextinction[0,*]= in.extinction[0]/5.155D0*0.184
    wiseextinction[1,*]= in.extinction[0]/5.155D0*0.113
 ENDIF
ENDIF ELSE BEGIN 
 wise=0
 wiseflux=0.
 wiseflux_ivar=0.
 wiseextinction=0.
ENDELSE

IF (keyword_set(galex) OR keyword_set(ukidss) OR keyword_set(wise)) THEN BEGIN
    xdqso_combine_fluxes, flux, flux_ivar, extinction, anirflux=nirflux, $
      bnirflux_ivar=nirflux_ivar, $
      cnirextinction= nirextinction, duvflux=uvflux, $
      euvflux_ivar=uvflux_ivar, fuvextinction=uvextinction, $
      gwiseflux=wiseflux,hwiseflux_ivar=wiseflux_ivar,iwiseextinction=wiseextinction, $
      nir=ukidss,uv=galex, mir=wise, fluxout=outflux, ivarfluxout=outflux_ivar, $
      extinctionout=outextinction
      flux= outflux
      flux_ivar= outflux_ivar
      extinction= outextinction
ENDIF

IF ~keyword_set(dereddened) THEN BEGIN
    flux= xdqso_sdss_deredden(flux,extinction)
    flux_ivar=xdqso_sdss_deredden_error(flux_ivar,extinction)
ENDIF


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
out.qsomidzlike= xdqso_eval_colorprob(flux,flux_ivar,/qso,/midz,galex=galex,ukidss=ukidss,wise=wise)
out.qsomidznumber= xdqso_eval_iprob(flux[3,*],i_qsobossz,dndi_qsobossz)
out.qsohizlike= xdqso_eval_colorprob(flux,flux_ivar,/qso,galex=galex,ukidss=ukidss,wise=wise)
out.qsohiznumber= xdqso_eval_iprob(flux[3,*],i_qso,dndi_qso)
out.qsolowzlike= xdqso_eval_colorprob(flux,flux_ivar,/qso,/lowz,galex=galex,ukidss=ukidss,wise=wise)
out.qsolowznumber= xdqso_eval_iprob(flux[3,*],i_qsolowz,dndi_qsolowz)
out.starlike= xdqso_eval_colorprob(flux,flux_ivar,galex=galex,ukidss=ukidss,wise=wise)
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
