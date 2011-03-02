FUNCTION xdqso_sdss_deredden, flux, extinction
exponent = 0.4*extinction
flux_correct = (10.0^exponent)*flux
return, flux_correct
end
