;+
;   NAME:
;      xdqso_eval_prob
;   PURPOSE:
;      evaluate the number counts at this i-magnitude
;   INPUT:
;      fi - i PSF flux
;      i - i-mag prior (from xdqso_read_dndi)
;      dndi - dndi prior (from xdqso_read_dndi)
;   OUTPUT:
;      number count [deg^{-2}]
;   HISTORY:
;      2010 - Written - Bovy (NYU)
;-
FUNCTION XDQSO_EVAL_IPROB, fi, i, dndi
b= 1.8
nfi= n_elements(fi)
if nfi EQ 1 THEN scalarOut= 1B ELSE scalarOut= 0B
mi= xdqso_sdss_flux2mags(fi,b)
out= INTERPOL(dndi,i,mi,/spline)
IF scalarOut THEN RETURN, out[0] ELSE RETURN, out
END
