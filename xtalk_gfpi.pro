pro xtalk_gfpi, sif, sqf, suf, svf, maskv2q=maskv2q, maskv2u=maskv2u, maskqu2v=maskqu2v

tay = size(sif)

print, tay

npos = tay[1]

cnt_QV=0
cnt_UV=0
cnt_QU=0
maxqq=0
maxuu=0
maxvv=0
minqq=0
minuu=0
minvv=0

xi2quv=fltarr(3, tay[1], tay[2])
step = 3.0d-3
xh=findgen(80)*step
xh=xh-max(xh)/2.

print, max(xh)

nmax=10000
xq_V=fltarr(tay[3],nmax)
xu_V=fltarr(tay[3],nmax)
xv_QV=fltarr(tay[3],nmax)
xv_UV=fltarr(tay[3],nmax)
xq_QU=fltarr(tay[3],nmax)
xu_QU=fltarr(tay[3],nmax)
xv_QU=fltarr(tay[3],nmax)

print,'Xtalk I --> Q,U,V'
print,'*****************'


for j=0, tay[1]-1 do begin

   imi = median(reform(sif[j, *, *]), 5)
   imq = median(reform(sqf[j, *, *]), 5)/imi
   imu = median(reform(suf[j, *, *]), 5)/imi
   imv = median(reform(svf[j, *, *]), 5)/imi

;   imi = transpose(imi)
   imq = transpose(imq)
   imu = transpose(imu)
   imv = transpose(imv)

   for i=0,tay[2]-1 do begin

      hq=smooth(histogram(imq(*, i),binsize=step,min=min(xh),max=max(xh)), 3)
      hu=smooth(histogram(imu(*, i),binsize=step,min=min(xh),max=max(xh)), 3)
      hv=smooth(histogram(imv(*, i),binsize=step,min=min(xh),max=max(xh)), 3)
 
      hq = histogram(imq[*,i], binsize=step,min=min(xh),max=max(xh))
      hu = histogram(imu[*,i], binsize=step,min=min(xh),max=max(xh))
      hv = histogram(imv[*,i], binsize=step,min=min(xh),max=max(xh))

      xi2quv(0,j,i)=xh(min(where(hq eq max(hq))))
      xi2quv(1,j,i)=xh(min(where(hu eq max(hu))))
      xi2quv(2,j,i)=xh(min(where(hv eq max(hv))))

   endfor
endfor
print,' '

print,'xtalk I --> Q = ',mean(xi2quv(0,*,*)),' +/- ',stddev(xi2quv(0,*,*))
print,'xtalk I --> U = ',mean(xi2quv(1,*,*)),' +/- ',stddev(xi2quv(1,*,*))
print,'xtalk I --> V = ',mean(xi2quv(2,*,*)),' +/- ',stddev(xi2quv(2,*,*))


xi2q=reform(xi2quv(0,*,*)) < .01 > (-.01)
xi2u=reform(xi2quv(1,*,*)) < .01 > (-.01)
xi2v=reform(xi2quv(2,*,*)) < .01 > (-.01)


xi2q = median(xi2q, 5)
xi2u = median(xi2u, 5)
xi2v = median(xi2v, 5)

maxq=fltarr(tay[1], tay[2])
maxu=fltarr(tay[1], tay[2])
maxv=fltarr(tay[1], tay[2])

print,' '
print,'Checking Q,U,V profiles'	
print,'***********************'
for j=0, tay[1]-1 do begin

   imi = reform(sif[j, *, *])
   imq = reform(sqf[j, *, *])/imi
   imu = reform(suf[j, *, *])/imi
   imv = reform(svf[j, *, *])/imi

   imi = transpose(imi)
   imq = transpose(imq)
   imu = transpose(imu)
   imv = transpose(imv)

   for k=0, tay[2]-1 do begin
      imq(*, k) = imq(*, k)-xi2q(j,k)
      imu(*, k) = imu(*, k)-xi2u(j,k)
      imv(*, k) = imv(*, k)-xi2v(j,k)
   endfor

   for k=0, tay[2]-1 do begin
      maxq(j,k)=max(abs(imq(*, k)))
      maxu(j,k)=max(abs(imu(*, k)))
      maxv(j,k)=max(abs(imv(*, k)))
   endfor
 
endfor

print,' '

nmin=0.01*float(tay[2])*float(tay[1])
cota=2.0
print, 'cota = ', cota
print,' '
maskV2Q=intarr(tay[1], tay[2])
zV2Q=where(maxv gt cota*maxq)
if(zV2Q(0) ne -1 and n_elements(zV2Q) ge nmin) then begin
   maskV2Q(zV2Q)=1
   print,'V   --> Q ',n_elements(zV2Q),' profiles available'
endif else begin
   print,'not enough V --> Q profiles available' 
endelse


maskV2U=intarr(tay[1],tay[2])
zV2U=where(maxv gt cota*maxu)
if(zV2U(0) ne -1 and n_elements(zV2U) ge nmin) then begin
   maskV2U(zV2U)=1
   print,'V   --> U ',n_elements(zV2U),' profiles available'
endif else begin
   print,'not enough V --> U profiles available' 
endelse

maskQU2V=intarr(tay[1],tay[2])
zQU2V=where(maxq gt cota*maxv and maxu gt cota*maxv)
if(zQU2V(0) ne -1 and n_elements(zQU2V) ge nmin) then begin
   maskQU2V(zQU2V)=1
   print,'Q,U --> V ',n_elements(zQU2V),' profiles available'
endif else begin
   print,'not enough Q,U --> V profiles available' 
endelse
      
print,' '

np = tay[3]
if(zV2Q(0) ne -1 and n_elements(zV2Q) ge nmin) then xxV2Q=fltarr(np*n_elements(zV2Q))
if(zV2U(0) ne -1 and n_elements(zV2U) ge nmin) then xxV2U=fltarr(np*n_elements(zV2U))
if(zQU2V(0) ne -1 and n_elements(zQU2V) ge nmin) then xxQU2V=fltarr(np*n_elements(zQU2V),2)

if(zV2Q(0) ne -1 and n_elements(zV2Q) ge nmin) then yV2Q=fltarr(np*n_elements(zV2Q))
if(zV2U(0) ne -1 and n_elements(zV2U) ge nmin) then yV2U=fltarr(np*n_elements(zV2U))
if(zQU2V(0) ne -1 and n_elements(zQU2V) ge nmin) then yQU2V=fltarr(np*n_elements(zQU2V))


print,'Xtalk V --> Q,U and Xtalk Q,U --> V'
print,'***********************************'

j1_V2Q=0
j1_V2U=0
j1_QU2V=0

for j=0,tay[1]-1 do begin
   cntQ2V=0
   cntU2V=0
   imi = reform(sif[j, *, *])
   imq = reform(sqf[j, *, *])/imi
   imu = reform(suf[j, *, *])/imi
   imv = reform(svf[j, *, *])/imi

   imi = transpose(imi)
   imq = transpose(imq)
   imu = transpose(imu)
   imv = transpose(imv)
   
   for k=0,tay[2]-1 do begin
      imq(*, k)=imq(*, k)-xi2q(j,k)
      imu(*, k)=imu(*, k)-xi2u(j,k)
      imv(*, k)=imv(*, k)-xi2v(j,k)
   endfor

   z=where(maskV2Q(j,*) eq 1)
   if(z(0) ne -1) then begin
      j2_V2Q=j1_V2Q+n_elements(z)*np-1
      xxV2Q(j1_V2Q:j2_V2Q)=imv(0:np-1,z)
      yV2Q(j1_V2Q:j2_V2Q)=imq(0:np-1,z)
      j1_V2Q=j2_V2Q+1
   endif   

   z=where(maskV2U(j,*) eq 1)
   if(z(0) ne -1) then begin
      j2_V2U=j1_V2U+n_elements(z)*np-1
      xxV2U(j1_V2U:j2_V2U)=imv(0:np-1,z)
      yV2U(j1_V2U:j2_V2U)=imu(0:np-1,z)
      j1_V2U=j2_V2U+1
   endif   

   z=where(maskQU2V(j,*) eq 1)
   if(z(0) ne -1) then begin
      j2_QU2V=j1_QU2V+n_elements(z)*np-1
      xxQU2V(j1_QU2V:j2_QU2V,0)=imq(0:np-1,z)
      xxQU2V(j1_QU2V:j2_QU2V,1)=imu(0:np-1,z)
      yQU2V(j1_QU2V:j2_QU2V)=imv(0:np-1,z)
      j1_QU2V=j2_QU2V+1
   endif   
endfor

ccV2Q=fltarr(2)
ccV2U=fltarr(2)
ccQU2V=fltarr(2,2)
if(zV2Q(0) ne -1 and n_elements(zV2Q) ge nmin) then ccV2Q=lstsqfit(xxV2Q,yV2Q)
if(zV2U(0) ne -1 and n_elements(zV2U) ge nmin) then ccV2U=lstsqfit(xxV2U,yV2U)
if(zQU2V(0) ne -1 and n_elements(zQU2V) ge nmin) then ccQU2V=lstsqfit(xxQU2V,yQU2V)

xxV2Q=0.
xxV2U=0.
xxQU2V=0.
yV2Q=0.
yV2U=0.
yQU2V=0.


print,' '
if(zV2Q(0) ne -1 and n_elements(zV2Q) ge nmin) then print,'Xtalk V --> Q = ',ccV2Q(0),' +/-',ccV2Q(1)
if(zV2U(0) ne -1 and n_elements(zV2U) ge nmin) then print,'Xtalk V --> U = ',ccV2U(0),' +/-',ccV2U(1)
if(zQU2V(0) ne -1 and n_elements(zQU2V) ge nmin) then begin
   print,'Xtalk Q --> V = ',ccQU2V(0,0),' +/-',ccQU2V(0,1)
   print,'Xtalk U --> V = ',ccQU2V(1,0),' +/-',ccQU2V(1,1)
endif
print,' '



print,'Correcting for Xtalk and saving data'
print,'************************************'
for j=0, tay[1]-1 do begin
   imi = reform(sif[j, *, *])
   imq = reform(sqf[j, *, *])
   imu = reform(suf[j, *, *])
   imv = reform(svf[j, *, *])

   imi = transpose(imi)
   imq = transpose(imq)
   imu = transpose(imu)
   imv = transpose(imv)
  
   for k=0, tay[2]-1 do begin

     imq(*, k)=imq(*, k) - xi2q(j,k) *imi(*, k) - ccV2Q(0)*imv(*, k)
     imu(*, k)=imu(*, k) - xi2u(j,k) *imi(*, k) - ccV2U(0)*imv(*, k)
     imv(*, k)=imv(*, k) - xi2v(j,k) *imi(*, k) - ccQU2V(0,0)*imq(*, k)-  ccQU2V(1,0)*imu(*, k)

   endfor
  sqf[j, *, *] = transpose(imq)
  suf[j, *, *] = transpose(imu)
  svf[j, *, *] = transpose(imv)

endfor

print,' '
print,'Xtalk correction finished !'	
print,'****************************'

end      	    

