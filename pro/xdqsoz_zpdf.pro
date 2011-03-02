;+
;   NAME:
;      xdqsoz_zpdf
;   PURPOSE:
;      calculate the photometric redshift pdf for XDQSOZ
;   USE:
;      xdqsoz_zpdf, flux, flux_ivar, /galex, /ukidss, zmean=zmean,
;      zcovar=zcovar, zamp=zamp
;   INPUT:
;      flux - [nfluxes-1] or [nfluxes-1,ndata] array of fluxes
;      flux_ivar - [nfluxes-1] or [nfluxes-1,ndata] array of flux_ivars
;   KEYWORDS:
;      galex - use GALEX fits
;      ukidss - use UKIDSS
;   OUTPUT:
;      zmean - [ngauss,ndata] array of means
;      zcovar - [ngauss,ndata] array of covars
;      zamp - [ngauss,ndata] array of amplitudes
;   HISTORY:
;      2011-01-18 - Written - Bovy (NYU)
;-
PRO XDQSOZ_ZPDF, flux, flux_ivar, galex=galex, ukidss=ukidss, $
                 zmean=zmean, zcovar=zcovar, zamp=zamp
;;check for environment variable
path= getenv('XDQSODATA')
if strcmp(path,'') then _SAVEDIR= '../data/' else _SAVEDIR = '$XDQSODATA/'
savefilename= _SAVEDIR+'xdqsoz_relflux_fits'
IF keyword_set(galex) THEN savefilename+= '_galex'
IF keyword_set(ukidss) THEN savefilename+= '_ukidss'
savefilename+= '.fits'

b= 1.8;;Magnitude softening
_IMIN= 17.7
_IMAX= 22.5
_ISTEP= 0.1
_IWIDTH= 0.2
_NGAUSS= 60
nbins= (_IMAX-_IMIN)/_ISTEP

nfi= n_elements(flux[0,*])
ndim= n_elements(flux)/nfi-1

if nfi EQ 1 THEN scalarOut= 1B ELSE scalarOut= 0B
mi= xdqso_sdss_flux2mags(flux[3,*],b)
zamp= dblarr(_NGAUSS,nfi)
zmean= dblarr(_NGAUSS,nfi)
zcovar= dblarr(_NGAUSS,nfi)
;;Just loop through the solutions bin
FOR ii=0L, nbins-1 DO BEGIN
    indx= where(mi GE (_IMIN+(ii+0.5)*_ISTEP) AND $
                mi LT (_IMIN+(ii+1.5)*_ISTEP))
    IF indx[0] EQ -1 THEN CONTINUE ;;Nothing here
    ;;Prep the data
    if scalarOut THEN xdqso_prep_data, flux, flux_ivar, mags=ydata,var_mags=ycovar, $
      /relfluxes ELSE $
      xdqso_prep_data, flux[*,indx], flux_ivar[*,indx], mags=ydata,var_mags=ycovar, $
      /relfluxes
    ndata= n_elements(ydata[0,*])
    ;;Load solution
    fits= mrdfits(savefilename,ii+1,/silent)
    ;;Marginalize over redshift first
    ngauss= n_elements(fits.amp)
    ndimx= n_elements(fits.xmean)/ngauss
    thisxmean= fits.xmean[1:ndimx-1,*]
    thisxcovar= fits.xcovar[1:ndimx-1,1:ndimx-1,*]
    xdqso_calc_membership_prob, thismember_prob, ydata, ycovar,thisxmean,thisxcovar,$
      fits.xamp, loglike=thisnorm
    zamp[*,indx]= exp(thismember_prob)
    ;;Now get the means and covars
    fluxxmean= fits.xmean[1:ndimx-1,*]
    fluxxcovar= fits.xcovar[1:ndimx-1,1:ndimx-1,*]
    zxmean= fits.xmean[0,*]
    zxcovar= fits.xcovar[0,0,*]
    zfluxxcovar= fits.xcovar[0,1:ndimx-1,*]
    for jj=0L, n_elements(indx)-1 do begin
        for kk=0L, ngauss-1 do begin
            fluxinvxcovar= invert(fluxxcovar[*,*,kk]+ycovar[*,*,jj],/double)
            zmean[kk,indx[jj]]= zxmean[0,kk]+zfluxxcovar[0,*,kk]#(fluxinvxcovar#(ydata[*,jj]-fluxxmean[*,kk]))
            zcovar[kk,indx[jj]]= zxcovar[0,0,kk]-zfluxxcovar[0,*,kk]#(fluxinvxcovar#transpose(zfluxxcovar[0,*,kk]))
        endfor
    endfor
ENDFOR
if scalarOut then begin
    zmean= zmean[*,0]
    zamp= zamp[*,0]
    zcovar= zcovar[*,0]
endif
END
