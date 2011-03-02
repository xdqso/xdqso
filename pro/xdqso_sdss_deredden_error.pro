FUNCTION xdqso_sdss_deredden_error, ivar2, extinction
exponent = -0.8*extinction
ivar2_correct = (10.0^exponent)*ivar2
return, ivar2_correct
end
