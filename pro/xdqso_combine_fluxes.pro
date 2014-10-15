;+
;   NAME:
;      combine_fluxes
;   PURPOSE:
;      combine the sdss, NIR, UV, and MIR fluxes
;   INPUT:
;      flux - sdss flux
;      flux_ivar - sdss flux_ivar
;      extinction - sdss extinction
;      nirflux - UKIDSS flux
;      nirflux_ivar - UKIDSS flux_ivar
;      nirextinction -UKIDSS extinction
;      uvflux - GALEX flux
;      uvflux_ivar - GALEX flux_ivar
;      uvextinction - GALEX extinction
;      wiseflux - WISE MIR flux
;      wiseflux_ivar - WISE MIR flux_ivar
;      wiseextinction - WISE extinction
;   KEYWORDS:
;      nir - NIR included
;      uv - UV included
;      mir - WISE included
;   OUTPUT:
;      outflux, outflux_ivar, outextinction - combination of inputs
;   HISTORY:
;      2010-05-05 - Written - Bovy (NYU)
;      2014-03-19 - Updated for WISE - DiPompeo (UWyo)
;-
PRO XDQSO_COMBINE_FLUXES, flux, flux_ivar, extinction, anirflux=anirflux, $
                    bnirflux_ivar=bnirflux_ivar, $
                    cnirextinction= cnirextinction, duvflux=duvflux, $
                    euvflux_ivar=euvflux_ivar, fuvextinction=fuvextinction, $
                    gwiseflux=gwiseflux, hwiseflux_ivar=hwiseflux_ivar,iwiseextinction=iwiseextinction, $
                    nir=nir,uv=uv, mir=mir, fluxout=fluxout, $
                    ivarfluxout=ivarfluxout, extinctionout=extinctionout
_NSDSS= 5
_NNIR= 4
_NUV= 2
_NWISE=2
IF keyword_set(nir) AND keyword_set(uv) AND keyword_set(mir) THEN BEGIN
    nfluxes= _NSDSS+_NNIR+_NUV+_NWISE
ENDIF ELSE IF keyword_set(nir) AND keyword_set(uv) THEN BEGIN
    nfluxes= _NSDSS+_NNIR+_NUV
ENDIF ELSE IF keyword_set(nir) AND keyword_set(mir) THEN BEGIN
    nfluxes= _NSDSS+_NNIR+_NWISE
ENDIF ELSE IF keyword_set(uv) AND keyword_set(mir) THEN BEGIN
    nfluxes= _NSDSS+_NUV+_NWISE
ENDIF ELSE IF keyword_set(uv) THEN BEGIN
    nfluxes= _NSDSS+_NUV
ENDIF ELSE IF keyword_set(nir) THEN BEGIN
    nfluxes= _NSDSS+_NNIR
ENDIF ELSE IF keyword_set(mir) THEN BEGIN
    nfluxes= _NSDSS+_NWISE
ENDIF ELSE BEGIN
    nfluxes= _NSDSS
ENDELSE
ndata= n_elements(flux[0,*])
outflux= dblarr(nfluxes,ndata)
outflux_ivar= dblarr(nfluxes,ndata)
outextinction= dblarr(nfluxes,ndata)

outflux[0:4,*]= flux[*,*]
outflux_ivar[0:4,*]= flux_ivar[*,*]
outextinction[0:4,*]= extinction[*,*]
IF keyword_set(uv) THEN BEGIN
    outflux[5:6,*]= duvflux[*,*]
    outflux_ivar[5:6,*]= euvflux_ivar[*,*]
    outextinction[5:6,*]= fuvextinction[*,*]
    IF (keyword_set(nir) AND keyword_set(mir)) THEN BEGIN
        outflux[7:10,*]= anirflux[*,*]
        outflux_ivar[7:10,*]= bnirflux_ivar[*,*]
        outextinction[7:10,*]= cnirextinction[*,*]
        outflux[11:12,*]=gwiseflux[*,*]
        outflux_ivar[11:12,*]= hwiseflux_ivar[*,*]
        outextinction[11:12,*]= iwiseextinction[*,*]
    ENDIF ELSE IF keyword_set(nir) THEN BEGIN
        outflux[7:10,*]= anirflux[*,*]
        outflux_ivar[7:10,*]= bnirflux_ivar[*,*]
        outextinction[7:10,*]= cnirextinction[*,*]
    ENDIF ELSE IF keyword_set(mir) THEN BEGIN
        outflux[7:8,*]=gwiseflux[*,*]
        outflux_ivar[7:8,*]= hwiseflux_ivar[*,*]
        outextinction[7:8,*]= iwiseextinction[*,*]
    ENDIF 
ENDIF ELSE IF keyword_set(nir) THEN BEGIN
    outflux[5:8,*]= anirflux[*,*]
    outflux_ivar[5:8,*]= bnirflux_ivar[*,*]
    outextinction[5:8,*]= cnirextinction[*,*]
    IF keyword_set(mir) THEN BEGIN
         outflux[9:10,*]=gwiseflux[*,*]
         outflux_ivar[9:10,*]= hwiseflux_ivar[*,*]
         outextinction[9:10,*]= iwiseextinction[*,*]
    ENDIF
ENDIF ELSE IF keyword_set(mir) THEN BEGIN
    outflux[5:6,*]=gwiseflux[*,*]
    outflux_ivar[5:6,*]= hwiseflux_ivar[*,*]
    outextinction[5:6,*]= iwiseextinction[*,*]
ENDIF
fluxout= outflux
ivarfluxout= outflux_ivar
extinctionout= outextinction
END


