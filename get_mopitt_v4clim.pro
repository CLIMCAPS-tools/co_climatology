pro get_mopitt_v4clim,datestr,alat,coprior

; This is the IDL version of get_mopitt_v4clim.F

; read 'mopitt_v4clim_nhsh_20180605.dat' and make it into 1 x 1 deg grid

; mopitt_clim_co has dimensions [nmon,nlon,nlat,npres]
; where nmon=1, nlon=128, nlat=64, npres=100
read_mopitt_v4clim_binary,mopitt_clim_co

; mopitt_clim_co has place for more profiles but is populated 
; with only one profile per hemisphere 
; get date
yy=fix(strmid(datestr,0,4))
mm=fix(strmid(datestr,4,2))
dd=fix(strmid(datestr,6,2))

codim=size(mopitt_clim_co,/dimensions)

mopitt_clim_nmon=codim(0); 12
mopitt_clim_nlat=codim(1); 2
mopitt_clim_npres=codim(2); 100

; Calculate which climatology profile to use based on input [ilat,ilon]
; -------------------------------------------
; This is necessary if mopitt_clim_nlat > 2
;dlat = 180. / FLOAT(mopitt_clim_nlat - 1)
;ilat = FIX((alat + 90.) / dlat + 0.51 )

;if (ilat ge mopitt_clim_nlat) then ilat=mopitt_clim_nlat-1
; Here, ilat is either 0 or 1
ilat=0
; ****** Do this when nlat > 2
; Linear interpolation weighting for latitude
; -------------------------------------------
; wgt_lat = (-90.0+ dlat*float(ilat)-alat)/dlat
; wgt_lat = (1.0 - wgt_lat)

; ****** Do this when nlat = 2
; Two profile nh/sh from J. Warner UMBC 
; interpolates between +/- 15 deg. 
; -------------------------------------------
lat_transit = 30.
lat_boundary = 0.
if (alat lt -lat_transit/2.0) then wgt_lat = 0.
if (alat gt lat_transit/2.0) then wgt_lat = 1.
if (alat ge -lat_transit/2.0) and (alat le lat_transit/2.0) then $
	wgt_lat = ABS(alat -(lat_boundary-lat_transit/2.))/lat_transit

; There is only one Longitudinal value
ilon = 1

; ==================================
; interpolate the co profile in time 
; ----------------------------------     
days=[0.,31.,59.,90.,120.,151.,181.,212.,$
	243.,273.,304.,334]

; Calculate the decimal julian day for the current observation
; - Test if this is a leap year
i = yy/4
i = i*4 
rdays = 365.

day = float(dd) + days(fix(mm)-1)
if (i eq yy) then begin ; Leap year
	rdays = 366.
	if mm gt 2 then day = day + 1.0
endif

time = day / rdays

; The mopitt climatology is monthly averaged with the date
; used as the center of each month
; Find the index of the month to interpolate
; --------------------------------------------
int_time1=-1
int_time2=-1
co_time=fltarr(12); mopitt_clim_maxtime=12 in clim.com
for imonth=0,mopitt_clim_nmon-2 do begin
		co_time(imonth)=(0.5*(days(imonth+1)-days(imonth))+days(imonth))/rdays
		if co_time(imonth) lt time then int_time1 = imonth
endfor

co_time(mopitt_clim_nmon-1)=(0.5*(rdays - days(mopitt_clim_nmon-1))+ $
	days(mopitt_clim_nmon-1))/rdays

if co_time(mopitt_clim_nmon-1) le time then int_time1=mopitt_clim_nmon-1

int_time2 = int_time1 + 1

if int_time1 eq mopitt_clim_nmon - 1 then begin
	int_time1 = mopitt_clim_nmon - 1 
	int_time2 = 0 
endif

if int_time1 eq -1 then begin
	int_time1 = mopitt_clim_nmon - 1
	int_time2 = 0
endif

if int_time1 gt int_time2 then begin
	co_time(int_time1) = co_time(int_time1)-1.0
	if mm eq 12 then time = time - 1.0
endif

coprior=fltarr(mopitt_clim_npres)
dtm = (co_time(int_time2) - co_time(int_time1))

for ip=0,mopitt_clim_npres-1 do begin
; read the two CO profiles to combine with wgt_lat in final answer
	dco1=(mopitt_clim_co(int_time2,0,ip)- $
		mopitt_clim_co(int_time1,0,ip))
	dco2=(mopitt_clim_co(int_time2,1,ip)- $
		mopitt_clim_co(int_time1,1,ip))

	coprior(ip) = (1.0-wgt_lat)*(mopitt_clim_co(int_time1,0,ip)+ $
		dco1 / dtm * (time - co_time(int_time1))) + $
		wgt_lat * (mopitt_clim_co(int_time1,1,ip) + $
		dco2 / dtm * (time - co_time(int_time1)))

; convert from parts to parts per billion
	coprior(ip)=coprior(ip)*1.e9
endfor
end
