FUNCTION XDQSO_DNDIPATH, zmin, zmax, lumfunc
;;check for environment variable
path= getenv('XDQSODATA')
if strcmp(path,'') then path= '../data/' else path = '$XDQSODATA/'
path+= 'dNdi_zmin_'+strtrim(string(zmin,format='(F4.2)'),2)+'_zmax_'+strtrim(string(zmax,format='(F4.2)'),2)+'_'+lumfunc+'.prt'
return, path
END
