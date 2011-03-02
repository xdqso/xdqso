;+
;   NAME:
;      xdqsoz_eval_zpdf
;   PURPOSE:
;      evaluate the photometric redshift PDF for a given redshift
;      given means, covars, and amps
;   USE:
;      pz= xdqsoz_eval_zpdf(z,zmean,zcovar,zamp)
;   INPUT:
;      z - redshift [nz]
;      zmean, zcovar, zamp - from xdqsoz_zpdf.pro
;   OUTPUT:
;      p(z)
;   HISTORY:
;      2011-01-18 - Written - Bovy (NYU)
;-
FUNCTION XDQSOZ_EVAL_ZPDF, z, zmean, zcovar, zamp
nz= n_elements(z)
if nz eq 1 then scalarOut= 1B else scalarOut= 0B
thisz= alog(z)
jac= 1D0/z
out= jac*exp(calc_loglike(reform(thisz,1,nz),dblarr(1,nz),reform(zmean,1,n_elements(zamp)),$
                          reform(zcovar,1,1,n_elements(zamp)),zamp))
IF scalarOut then return, out[0] else return, out
END

