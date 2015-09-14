;+
;   NAME:
;      xdqsoz_peaks
;   PURPOSE:
;      calculate the number of peaks of a zpdf as well as the MAP z
;   INPUT:
;      flux - dereddened flux
;      flux_ivar - dereddened flux_ivar
;   OPTIONAL INPUT:
;      nzs - number of points to sample the PDF at
;      peak_threshold - threshold for defining a peak (contiguous
;                       region with p above peak_threshold)
;   KEYWORDS:
;      galex - use GALEX
;      ukidss - use UKIDSS
;      wise - use WISE
;      plot - make QS plot
;   OUTPUT:
;      number of peaks
;   OPTIONAL OUTPUT:
;      peakz - MAP z
;      xdqsoz - structure containing
;               {peakz,peakprob,peakfwhm,otherz,otherprob,otherfwhm}
;               for all peaks
;   HISTORY:
;      2011-01-18 - Written - Bovy (NYU)
;      2014-03-31 - Added WISE - DiPompeo (UWyo)
;-
FUNCTION FWHM_ZPDF, thiszs, zpdf
ii= n_elements(thiszs)
;;find zone
maxp= max(zpdf)
ii= 0L
while zpdf[ii] LT maxp/2. do ii+= 1
zone= thiszs[ii]
ii= n_elements(zpdf)-1
while zpdf[ii] LT maxp/2. do ii-= 1
ztwo= thiszs[ii]
return, (ztwo-zone)
END
FUNCTION XDQSOZ_PEAKS, flux, flux_ivar, galex=galex,ukidss=ukidss, wise=wise, $
                       peakz=peakz, nzs=nzs, plot=plot, $
                       xdqsoz=xdqsoz, peak_threshold=peak_threshold
IF ~keyword_set(nzs) THEN  nzs= 1001
zs= dindgen(nzs)/(nzs-1d0)*5.2+0.3
if ~keyword_set(peak_threshold) then peak_threshold= 1./5.2;;uniform distribution
xdqsoz_zpdf, flux,$
  flux_ivar,$
  galex=galex,ukidss=ukidss, wise=wise,$
  zmean=zmean,zcovar=zcovar,$
  zamp=zamp
zpdf= xdqsoz_eval_zpdf(zs,zmean,zcovar,zamp)
IF keyword_set(plot) THEN BEGIN
    djs_plot, zs, zpdf
    djs_oplot, [zs[0],zs[nzs-1]],[peak_threshold,peak_threshold]
ENDIF
;;Find peaks
peak_indx=  where(zpdf GE peak_threshold,cnt)
peak_indx= peak_indx[sort(peak_indx)]
if cnt eq 0 then begin
    print, "No peaks found, something must be horribly wrong ..."
    xdqsoz= create_struct('peakz',0.,'peakprob',0.,$
                          'peakfwhm',0.,$
                          'otherz',0.,$
                          'otherprob',0.,$
                          'otherfwhm',0.)
    return, 0
endif
npeaks= 0
for jj=0L, cnt-2 do begin
    if (peak_indx[jj+1] ne peak_indx[jj]+1) OR $ ;;end of a peak
      (jj+2) eq cnt then $ ;;end of last peak
      npeaks+= 1
endfor
if arg_present(peakz) or arg_present(xdqsoz) then begin
    maxp= max(zpdf,peakz)
    peakz= zs[peakz]
endif    

;;If xdqsoz output is set, calculate everything about each of the
;;peaks
IF arg_present(xdqsoz) THEN BEGIN
    if npeaks eq 0 then begin
        xdqsoz= create_struct('peakz',0.,'peakprob',0.,$
                              'peakfwhm',0.,$
                              'otherz',0.,$
                              'otherprob',0.,$
                              'otherfwhm',0.)
    ENDIF ELSE BEGIN
        allz= dblarr(npeaks)
        allprob= dblarr(npeaks)
        allfwhm= dblarr(npeaks)
        ;;loop over peaks
        peak_start= peak_indx[0]
        jj= 0L
        for ii=0L, npeaks-1 do begin
            ;;find end
            peak_end= peak_start
            while (peak_indx[jj+1] eq peak_indx[jj]+1) and $ ;;end of a peak
              (jj+2) ne cnt do begin ;;end of last peak
                jj+= 1L
                peak_end+= 1L
            endwhile
            jj+= 1L
            thiszpdf= zpdf[peak_start:peak_end]
            thiszs= zs[peak_start:peak_end]
            ;;find maximum etc.
            thisp= max(thiszpdf,thispeakz)
            allz[ii]= thiszs[thispeakz]
            ;;find fwhm
            allfwhm[ii]= fwhm_zpdf(thiszs,thiszpdf)
            ;;integrate to get the probability
            allprob[ii]= xdqsoz_marginalize_colorzprob(thiszs[0],thiszs[n_elements(thiszs)-1],flux,flux_ivar,galex=galex,ukidss=ukidss,wise=wise,norm=totprob)
            if ii ne (npeaks-1) then peak_start= peak_indx[jj+1]
        endfor
        ;;normalize probabilities
        allprob= allprob/totprob[0]
        ;;find absolute peak
        maxprob= max(allprob,indx)
        peakprob= (allprob[indx])[0]
        peakfwhm= (allfwhm[indx])[0]
        indx= where(allprob ne maxprob,cnt)
        if cnt gt 0 then begin
            otherz= allz[indx]
            otherprob= allprob[indx]
            otherfwhm= allfwhm[indx]
        endif else begin
            otherz= 0.
            otherprob= 0.
            otherfwhm= 0.
        endelse
        ;;create output structure
        xdqsoz= create_struct('peakz',peakz,'peakprob',peakprob,$
                              'peakfwhm',peakfwhm,$
                              'otherz',otherz,$
                              'otherprob',otherprob,$
                              'otherfwhm',otherfwhm)
    ENDELSE
ENDIF
return, npeaks
END
