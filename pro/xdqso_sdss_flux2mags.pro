;   convert fluxes to asinh mags (f=flux, b=softening)
FUNCTION xdqso_sdss_flux2mags, f, b

a = 1.08574
ln10_min10 = -23.02585
; there is a bug in asinh that makes it crash if you call it with only 
; one object
IF n_elements(f) GT 1 THEN m = -a*(ASINH(5.0*f/b)+ln10_min10 + ALOG(b)) $
ELSE  m = -a*(ASINH(5.0*[f]/b)+ln10_min10 + ALOG(b))    

RETURN, m
END
