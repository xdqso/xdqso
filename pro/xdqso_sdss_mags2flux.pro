FUNCTION xdqso_sdss_mags2flux, m, b
a = 1.08574
ln10_min10 = -23.02585
f= b/5.0*SINH(-m/a-ln10_min10-ALOG(b))
RETURN, f
END
