;+
;   NAME:
;      xdqsoz_qso_track
;   PURPOSE:
;      calculate the mean quasar locus
;   INPUT:
;      z - redshift or array of redshifts [N]
;   OPTIONAL INPUT:
;      i= dereddened i-band magnitude
;   OUTPUT:
;      mags[5,N]
;   HISTORY:
;      2011-04-01 - Written - Bovy (NYU)
;-
FUNCTION XDQSOZ_QSO_TRACK, z, i=i, galex=galex, ukidss=ukidss
IF ~keyword_set(i) THEN i=19.
thisz= alog(z)
;;check for environment variable
path= getenv('XDQSODATA')
if strcmp(path,'') then _SAVEDIR= '../data/' else _SAVEDIR = '$XDQSODATA/'
savefilename= _SAVEDIR+'xdqsoz_relflux_fits'
IF keyword_set(galex) THEN savefilename+= '_galex'
IF keyword_set(ukidss) THEN savefilename+= '_ukidss'
savefilename+= '.fits'

; softening parameters from EDR paper in units of 1.0e-10
; (Stoughton et al. 2002)
b_u = 1.4
b_g = 0.9
b_r = 1.2
b_i = 1.8
b_z = 7.4
bs= [b_u,b_g,b_r,b_i,b_z]

_IMIN= 17.7
_IMAX= 22.5
_ISTEP= 0.1
_IWIDTH= 0.2
_NGAUSS= 60
nbins= (_IMAX-_IMIN)/_ISTEP

nz= n_elements(thisz)
IF keyword_set(galex) AND keyword_set(ukidss) THEN BEGIN
   ndim= 11
ENDIF ELSE IF keyword_set(galex) THEN BEGIN
   ndim= 7
ENDIF ELSE IF keyword_set(ukidss) THEN BEGIN
   ndim= 9
ENDIF ELSE BEGIN
   ndim= 5
ENDELSE
out= dblarr(ndim-1,nz)
;;restore solution
bin= floor(10.*(double(i)-double(_IMIN)-double(_ISTEP/2.))/long(10.*double(_ISTEP)))
IF i LE 17.75 or i GT 22.45 THEN BEGIN
   print, "i out of range, must be between 17.75 and 22.45"
   print, "Returning ..."
   return, -1
ENDIF
;;Load solution
fits= mrdfits(savefilename,bin+1,/silent)
amp= fits.xamp
xmean= fits.xmean
xcovar= fits.xcovar
ngauss= n_elements(amp)
ndimx= n_elements(xmean)/ngauss
fluxxmean= xmean[1:ndimx-1,*]
fluxxcovar= xcovar[1:ndimx-1,1:ndimx-1,*]
zxmean= xmean[0,*]
zxcovar= xcovar[0,0,*]
zfluxxcovar= xcovar[0,1:ndimx-1,*]
;;First calculate membership probabilities
memz= dblarr(ndimx,nz)
memz[0,*]= thisz
varz= dblarr(ndimx,ndimx,nz)
for jj= 1L, ndimx-1 do varz[jj,jj,*]= 1D5
xdqso_calc_membership_prob, member_prob, memz, varz, xmean, xcovar, amp
for jj=0L, nz-1 do begin
   thisout= dblarr(ndim-1)
   for kk=0L, ngauss-1 do begin
      meanrflux= fluxxmean[*,kk]+zfluxxcovar[0,*,kk]/zxcovar[kk]*(thisz[jj]-zxmean[kk])
      thisout+= exp(member_prob[kk,jj])*meanrflux
   endfor
   out[*,jj]= thisout
ENDFOR
;;process output, i.e., convert to mags and put i in
tmpout= dblarr(ndim,nz)
out*= xdqso_sdss_mags2flux(i,b_i)
tmpout[0,*]= out[0,*];u
tmpout[1,*]= out[1,*];g
tmpout[2,*]= out[2,*];r
tmpout[3,*]= xdqso_sdss_mags2flux(i,b_i);i
for jj= 4L, ndim-1 do tmpout[jj,*]= out[jj-1,*]
mags= dblarr(ndim,nz)
;;asinh for ugriz
FOR ii=0L, 4 DO mags[ii,*]= xdqso_sdss_flux2mags(tmpout[ii,*],bs[ii])
FOR ii=5L, ndim-1 DO mags[ii,*]= xdqso_flux2mags(tmpout[ii,*])
return, mags
END
