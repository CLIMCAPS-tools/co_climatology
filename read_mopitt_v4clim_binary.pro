pro read_mopitt_v4clim_binary,mopitt_clim_co


fn='~/climcaps_v2p00/static_files/all_ret/climatology/mopitt_v4clim_nhsh_20180605.dat'

;mopitt_clim_maxtime = 12
;mopitt_clim_maxlon  = 128
;mopitt_clim_maxlat  = 64
;mopitt_clim_maxpres = 100

b100 = fltarr(100)
b1 = 0.
openr,lun,fn,/get_lun,/swap_endian

readu,lun,b100 ; header
mopitt_clim_version = FIX(b100[0]); version 2.0
mopitt_clim_npres   = FIX(b100[1]); 100 pressure layers
mopitt_clim_nlat    = FIX(b100[2]); 2 latitudes, 1xNH, 1xSH
mopitt_clim_nlon    = FIX(b100[3]); 1 longitude
mopitt_clim_nmonth  = FIX(b100[4]); 12 months

readu,lun,b100  ; spam

buffer=fltarr(mopitt_clim_npres)

mopitt_clim_pres = fltarr(100)
mopitt_clim_co = fltarr(mopitt_clim_nmonth, $
          mopitt_clim_nlon, mopitt_clim_nlat, $
          mopitt_clim_npres)

readu,lun,b100
mopitt_clim_pres=b100

for imonth=0,mopitt_clim_nmonth-1 do begin
    for ilon = 0,mopitt_clim_nlon -1 do begin
       for ipres = 0,mopitt_clim_npres-1 do begin
          readu,lun,buffer
          for ilat = 0, mopitt_clim_nlat -1 do begin
            mopitt_clim_co[imonth,ilon,ilat,ipres] = buffer[ilat]
          endfor
       endfor         
    endfor
endfor

; get rid of all single dimensions
mopitt_clim_co=reform(mopitt_clim_co); [12 x 2 x 100]
; Two CO profiles for each month of the year
; - one profile representing northern hemisphere CO
; - one profile representing southern hemisphere CO
free_lun,lun
end
