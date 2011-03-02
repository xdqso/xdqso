;+
;   NAME:
;      xdqsoz_calculate_prob
;   PURPOSE:
;      calculate the extreme-deconvolution probability ratio,
;      marginalizing over an arbitrary redshift range
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
;   OUTPUT:
;      out - structure containing pqso, ...
;   HISTORY:
;      2010-04-30 - Written - Bovy (NYU)
;      2010-05-29 - Added Galex - Bovy
;      2010-10-30 - Added UKIDSS - Bovy
;-
FUNCTION XDQSOZ_CALCULATE_PROB, in, z_min, z_max, $
                                galex=galex, $
                                ukidss=ukidss
_DEFAULTZMIN=0.3
_DEFAULTZMAX=5.5
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

;;Now calculate all of the factors in turn
qsolike= xdqsoz_marginalize_colorzprob(z_min,z_max,flux,flux_ivar,galex=galex,ukidss=ukidss,norm=allqsolike)
bossqsonumber= xdqso_eval_iprob(flux[3,*],i_qsobossz,dndi_qsobossz)
qsonumber= xdqso_eval_iprob(flux[3,*],i_qso,dndi_qso)
qsolowznumber= xdqso_eval_iprob(flux[3,*],i_qsolowz,dndi_qsolowz)
qsonumber= bossqsonumber+qsonumber+qsolowznumber
everythinglike= xdqso_eval_colorprob(flux,flux_ivar,galex=galex,ukidss=ukidss)
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

;;Create output structure
outStruct= {pqso:0.D, $
            pqso_allz:0.D, $
            allqsolike:0D,$
            qsonumber:0.D,$
            everythinglike:0.D,everythingnumber:0.D}
out= replicate(outStruct,ndata)
out.pqso= pqso
out.pqso_allz= pqsotot
out.qsonumber= qsonumber
out.allqsolike= allqsolike
out.everythinglike= everythinglike
out.everythingnumber= everythingnumber
RETURN, out
END
