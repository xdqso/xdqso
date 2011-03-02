; Converts ivars to magnitude errors, ivar=inverse variance on flux,
; F=flux, b=softening

FUNCTION xdqso_sdss_ivar2magerr, ivar, F, b

dm = 0.0*ivar
zero_ind = WHERE(ivar EQ 0.0, n_zero)
IF n_zero NE 0 THEN ivar[zero_ind] = 1.0 ; prevents problem of div by 0

dF = 1.0/sqrt(ivar)

a = 1.08574
five_ovr_b = 5.0/b
root = 1.0/sqrt((five_ovr_b*F)^2 + 1.0)

dm = abs(a*root*five_ovr_b*dF)

IF n_zero NE 0 THEN dm[zero_ind] = 1.0e8 ; big number indicates infinite error

RETURN, dm
END
