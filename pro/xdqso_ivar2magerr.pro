; Converts ivars to magnitude errors, F=flux
FUNCTION xdqso_ivar2magerr, ivar, F
a= 1.08574
dm = 0.0*ivar
zero_ind = WHERE(ivar EQ 0.0, n_zero)
IF n_zero NE 0 THEN ivar[zero_ind] = 1.0 ; prevents problem of div by 0

dF = 1.0/sqrt(ivar)

a = 1.08574
dm = abs(a*dF/F)

IF n_zero NE 0 THEN dm[zero_ind] = 1.0e8 ; big number indicates infinite error

RETURN, dm
END
