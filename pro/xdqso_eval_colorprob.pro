;+
;   NAME:
;      xdqso_eval_colorprob
;   PURPOSE:
;      evaluate the relative flux likelihood for various classes
;   INPUT:
;      flux - dereddened flux
;      flux_ivar - dereddened flux_ivar
;   KEYWORDS:
;      qso - class qso (hiz by default)
;      lowz - low-z quasars
;      bossz - mid-z quasars;
;      galex - use GALEX
;      ukidss - use UKIDSS
;   OUTPUT:
;      likelihood
;   HISTORY:
;      2010 - Written - Bovy (NYU)
;-
FUNCTION XDQSO_EVAL_COLORPROB, flux, flux_ivar, qso=qso, lowz=lowz, $
                               midz=midz, galex=galex, ukidss=ukidss
;;check for environment variable
path= getenv('XDQSODATA')
if strcmp(path,'') then _SAVEDIR= '../data/' else _SAVEDIR = '$XDQSODATA/'
IF keyword_set(qso) AND keyword_set(lowz) THEN BEGIN
    savefilename= _SAVEDIR+'xdqso_relflux_fits_qsolowz'
ENDIF ELSE IF keyword_set(qso) AND keyword_set(midz) THEN BEGIN
    savefilename= _SAVEDIR+'xdqso_relflux_fits_qsomidz'
ENDIF ELSE IF keyword_set(qso) THEN BEGIN
    savefilename= _SAVEDIR+'xdqso_relflux_fits_qsohiz'
ENDIF ELSE BEGIN
    savefilename= _SAVEDIR+'xdqso_relflux_fits_star'
ENDELSE
IF keyword_set(galex) THEN savefilename+= '_galex'
IF keyword_set(ukidss) THEN savefilename+= '_ukidss'
savefilename+= '.fits'

b= 1.8;;Magnitude softening
_IMIN= 17.7
_IMAX= 22.5
_ISTEP= 0.1
_IWIDTH= 0.2
_NGAUSS= 20
nbins= (_IMAX-_IMIN)/_ISTEP

nfi= n_elements(flux[0,*])
if nfi EQ 1 THEN scalarOut= 1B ELSE scalarOut= 0B
mi= xdqso_sdss_flux2mags(flux[3,*],b)
out= dblarr(nfi)
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
    ;;Load solution
    fits= mrdfits(savefilename,ii+1,/silent)
    out[indx]= exp(xdqso_calc_loglike(ydata,ycovar,fits.xmean,fits.xcovar,fits.xamp))
ENDFOR
RETURN, out
END
