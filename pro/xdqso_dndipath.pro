;+
;   NAME:
;      xdqso_dndipath
;   PURPOSE:
;      return the path to the file that contains the number count
;      prior for a given redshift interval
;   INPUT:
;      zmin, zmax - minimum and maximum redshift
;      lumfunc - luminosity function
;   OUTPUT:
;      path
;   HISTORY:
;      2010 - Written - Bovy (NYU)
;-
FUNCTION XDQSO_DNDIPATH, zmin, zmax, lumfunc
;;check for environment variable
path= getenv('XDQSODATA')
if strcmp(path,'') then path= path_sep(/parent)+path_sep()+'data' $
else path = '$XDQSODATA'
result= strpos(path,path_sep(),/reverse_search)
if result ne (strlen(path)-1) then path= path+path_sep()
path+= 'dNdi_zmin_'+strtrim(string(zmin,format='(F4.2)'),2)+'_zmax_'+strtrim(string(zmax,format='(F4.2)'),2)+'_'+lumfunc+'.prt'
return, path
END
