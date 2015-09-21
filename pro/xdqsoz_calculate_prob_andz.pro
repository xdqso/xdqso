;+
;   NAME:
;      xdqsoz_calculate_prob_andz
;   PURPOSE:
;      calculate the extreme-deconvolution probability ratio,
;      marginalizing over an arbitrary redshift range, and calculate z
;      PDF in bins of 0.01
;   USE:
;      out= xdqsoz_calculate_prob(in,zmin,zmax,/deredenned,/galex,/ukidss)
;   INPUT:
;      in - structure containing PSFFLUX, PSFFLUX_IVAR, EXTINCTION
;      zmin, zmax - lower, upper bound of redshift interval
;   OPTIONAL INPUT:
;      dereddened - psfflux, and psfflux_ivar is already dereddened
;   KEYWORDS:
;      galex - GALEX fluxes are included in psfflux, psfflux_ivar, and
;              extinction; use them
;      ukidss - use UKIDSS (like /galex)
;      wise - use WISE (like /galex)
;   OUTPUT:
;      out - structure containing pqso, ...
;   HISTORY:
;      2010-04-30 - Written - Bovy (NYU)
;      2010-05-29 - Added Galex - Bovy
;      2010-10-30 - Added UKIDSS - Bovy
;      2014-03-31 - Added WISE - DiPompeo (UWyo)
;      2014-06-30 - Wrapped in xdqsoz_zpdf.pro, modified to output
;                   structure with z PDF included
;-
FUNCTION XDQSOZ_CALCULATE_PROB_ANDZ, in, z_min, z_max, $
                                dereddened=dereddened, $
                                galex=galex, $
                                ukidss=ukidss, $
                                wise=wise
_DEFAULTZMIN=0.3
_DEFAULTZMAX=5.5

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

;;Now calculate all of the factors in turn
qsolike= xdqsoz_marginalize_colorzprob(z_min,z_max,flux,flux_ivar,galex=galex,ukidss=ukidss,wise=wise,norm=allqsolike)
bossqsonumber= xdqso_eval_iprob(flux[3,*],i_qsobossz,dndi_qsobossz)
qsonumber= xdqso_eval_iprob(flux[3,*],i_qso,dndi_qso)
qsolowznumber= xdqso_eval_iprob(flux[3,*],i_qsolowz,dndi_qsolowz)
qsonumber= bossqsonumber+qsonumber+qsolowznumber
everythinglike= xdqso_eval_colorprob(flux,flux_ivar,galex=galex,ukidss=ukidss,wise=wise)
everythingnumber= xdqso_eval_iprob(flux[3,*],i_everything,dndi_everything)

;;Calculate the probability that a target is a zmin <= z <= zmax QSO
pqso= (qsolike*qsonumber)
nonzero= where(pqso NE 0.)
IF nonzero[0] NE -1 THEN pqso[nonzero]= pqso[nonzero]/(everythinglike[nonzero]*everythingnumber[nonzero]+$
                                                       allqsolike[nonzero]*qsonumber[nonzero])
pqsotot= (allqsolike*qsonumber)
nonzero= where(pqsotot NE 0.)
IF nonzero[0] NE -1 THEN pqsotot[nonzero]= pqsotot[nonzero]/(everythinglike[nonzero]*everythingnumber[nonzero]+$
                                                       allqsolike[nonzero]*qsonumber[nonzero])

;MAD Calculate P(z) for range of z
xdqsoz_zpdf,flux,flux_ivar,galex=galex,ukidss=ukidss,wise=wise,zmean=zmean,zcovar=zcovar,zamp=zamp
zarray=(findgen(z_max*100)/100)+0.01
zarray=zarray[where(zarray GE z_min)]
pz=dblarr(n_elements(zarray),ndata)
FOR i=0L,ndata-1 DO BEGIN
  pz[*,i]=xdqsoz_eval_zpdf(zarray,zmean[*,i],zcovar[*,i],zamp[*,i])
ENDFOR

;;Create output structure
outStruct= {pqso:0.D, $
            pqso_allz:0.D, $
            allqsolike:0D,$
            qsonumber:0.D,$
            everythinglike:0.D,everythingnumber:0.D, $
            z:zarray, zpdf:dblarr(n_elements(zarray))}

out= replicate(outStruct,ndata)
out.pqso= pqso
out.pqso_allz= pqsotot
out.qsonumber= qsonumber
out.allqsolike= allqsolike
out.everythinglike= everythinglike
out.everythingnumber= everythingnumber
out.zpdf=pz

RETURN, out
END
