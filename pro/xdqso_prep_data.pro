;+
;   NAME:
;      xdqso_prep_data
;   PURPOSE:
;      Prepare the data for analysis: deredden and convert to
;      magnitudes or colors
;   INPUT:
;      flux - 
;      flux_ivar - inverse error-squared on the flux
;      extinction - extinction along the los (optional)
;   KEYWORDS:
;      colors - if set, return colors and color errors (u-g, g-r, r-i,
;               i-z)
;      nbyncovar - if set, return the uncertainty variance as an n by
;                  n matrix (with correlations set to 0) instead of an
;                  n-array (ONLY FOR MAGS, SINCE THE COLORS ARE
;                  CORRELATED THEY ALWAYS RETURN AN NbyN MATRIX)
;      fluxes - if set, return the fluxes
;      relfluxes - if set, return relative fluxes to i-band
;   OUTPUT:
;      mags - magnitudes (or colors if /colors is set) [nmags,ndata]
;      var_mags - magnitude errors-squared (or color errors if /colors
;                 is set) [nmags,ndata]
;   HISTORY:
;      2010-03-04 - Written based on Hennawi's code snippets - Bovy (NYU)
;      2010-04-17 - Added relfluxes keyword - Bovy
;      2010-05-05 - Added nir and uv keywords and took them out again
;                   - Bovy
;-
PRO XDQSO_PREP_DATA, flux, flux_ivar, extinction=extinction, $
                     mags=mags, var_mags=var_mags, $
                     colors=colors, nbyncovar=nbyncovar, fluxes=fluxes, $
                     relfluxes=relfluxes
_BIGVAR= 1D5
; softening parameters from EDR paper in units of 1.0e-10
; (Stoughton et al. 2002)
b_u = 1.4
b_g = 0.9
b_r = 1.2
b_i = 1.8
b_z = 7.4
bs= [b_u,b_g,b_r,b_i,b_z]

if keyword_set(extinction) THEN BEGIN
    flux= xdqso_sdss_deredden(flux,extinction)
    flux_ivar= xdqso_sdss_deredden_error(flux_ivar,extinction)
ENDIF

IF keyword_set(fluxes) THEN BEGIN
    nobjs= n_elements(flux[0,*])
    nfluxes= n_elements(flux[*,0])
    mags= dblarr(nfluxes,nobjs)
    if keyword_set(nbyncovar) THEN var_mags= dblarr(nfluxes,nfluxes,nobjs) ELSE $
      var_mags= dblarr(nfluxes,nobjs)
    mags= flux
    IF keyword_set(nbyncovar) THEN BEGIN
        for jj=0L, nfluxes-1 DO var_mags[jj,jj,*]= 1D0/flux_ivar[jj,*]
    ENDIF ELSE BEGIN
        var_mags= 1D0/flux_ivar
    ENDELSE
    RETURN
ENDIF

IF keyword_set(relfluxes) THEN BEGIN
    nobjs= n_elements(flux[0,*])
    nfluxes= n_elements(flux[*,0])
    mags= dblarr(nfluxes-1,nobjs)
    var_mags= dblarr(nfluxes-1,nfluxes-1,nobjs)
    f_i= flux[3,*]
    mags[0,*]= flux[0,*]/f_i
    mags[1,*]= flux[1,*]/f_i
    mags[2,*]= flux[2,*]/f_i
    FOR jj=0L, 2 DO mags[jj,*]= flux[jj,*]/f_i
    FOR jj=4L, nfluxes-1 DO mags[jj-1,*]= flux[jj,*]/f_i

    ;;Diagonal elements
    FOR jj=0L, 2 DO var_mags[jj,jj,*]= 1D0/flux_ivar[jj,*]/f_i^2D0+flux[jj,*]^2D0/flux_ivar[3,*]/f_i^4D0
    FOR jj=4L, nfluxes-1 DO var_mags[jj-1,jj-1,*]= 1D0/flux_ivar[jj,*]/f_i^2D0+flux[jj,*]^2D0/flux_ivar[3,*]/f_i^4D0

    ;;Off-diagonal elements
    FOR jj=0L, 2 DO FOR kk=jj+1, 2L DO BEGIN
        var_mags[jj,kk,*]= flux[jj,*]*flux[kk,*]/flux[3,*]^4D0/flux_ivar[3,*]
        var_mags[kk,jj,*]= flux[jj,*]*flux[kk,*]/flux[3,*]^4D0/flux_ivar[3,*]
    ENDFOR
    FOR jj=4L, nfluxes-1 DO FOR kk=jj+1, nfluxes-1 DO BEGIN
        var_mags[jj-1,kk-1,*]= flux[jj,*]*flux[kk,*]/flux[3,*]^4D0/flux_ivar[3,*]
        var_mags[kk-1,jj-1,*]= flux[jj,*]*flux[kk,*]/flux[3,*]^4D0/flux_ivar[3,*]
    ENDFOR
    FOR jj=0L, 2 DO FOR kk=4L, nfluxes-1 DO BEGIN
        var_mags[jj,kk-1,*]= flux[jj,*]*flux[kk,*]/flux[3,*]^4D0/flux_ivar[3,*]
        var_mags[kk-1,jj,*]= flux[jj,*]*flux[kk,*]/flux[3,*]^4D0/flux_ivar[3,*]
    ENDFOR
    ;;fix infinities
    indx= where(finite(var_mags,/infinity),cnt)
    if cnt gt 0 then var_mags[indx]= _BIGVAR
    RETURN
ENDIF

;;Calculate magnitudes and uncertainties
nobjs= n_elements(flux[0,*])
nfluxes= n_elements(flux[*,0])
mags= dblarr(nfluxes,nobjs)
if keyword_set(nbyncovar) THEN var_mags= dblarr(nfluxes,nfluxes,nobjs) ELSE $
  var_mags= dblarr(nfluxes,nobjs)
FOR ii=0L, 4 DO mags[ii,*]= xdqso_sdss_flux2mags(flux[ii,*],bs[ii])
FOR ii=5L, nfluxes-1 DO mags[ii,*]= xdqso_flux2mags(flux[ii,*])
IF keyword_set(nbyncovar) THEN BEGIN
    FOR ii=0L, 4 DO var_mags[ii,ii,*]= xdqso_sdss_ivar2magerr(flux_ivar[ii,*],flux[ii,*],bs[ii])^2D0
    FOR ii=5L, nfluxes-1 DO var_mags[ii,ii,*]= xdqso_ivar2magerr(flux_ivar[ii,*],flux[ii,*],bs[ii])^2D0
ENDIF ELSE BEGIN
    FOR ii=0L, 4 DO var_mags[ii,*]= xdqso_sdss_ivar2magerr(flux_ivar[ii,*],flux[ii,*],bs[ii])^2D0
    FOR ii=5L, nfluxes-1 DO var_mags[ii,*]= xdqso_ivar2magerr(flux_ivar[ii,*],flux[ii,*])^2D0
ENDELSE
IF keyword_set(colors) THEN BEGIN
    colors= dblarr(nfluxes-1,nobjs)
    var_colors=dblarr(nfluxes-1,nfluxes-1,nobjs)
    FOR ii=0L, nfluxes-2 DO colors[ii,*]= mags[ii,*]- mags[ii+1,*]

    ;;covariance matrix, diagonal part
    FOR ii=0L, nfluxes-2 DO var_colors[ii,ii,*]= var_mags[ii,*] + var_mags[ii+1,*]
    ;;Off-diagonal elements
    FOR ii=0L, nfluxes-3 DO BEGIN
        var_colors[ii,ii+1,*]= -var_mags[ii+1,*]
        var_colors[ii+1,ii,*]= var_colors[ii,ii+1,*]
    ENDFOR
    
    mags= colors
    var_mags= var_colors
ENDIF
END
