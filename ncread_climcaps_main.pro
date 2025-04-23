pro ncread_climcaps_main,fn,outgroup

; Nadia Smith
; nadias@stcnet.com
; Aug 2019
; Routine to read NASA L2 CLIMCAPS netCDF file
; main branch of variables, not the sub-groups
;
; input:
; - fn: filename with full path
; output: 
; - outmain: structure containing all the variables in 
;   the main branch
; ----------------------------------------

; open ncdf file
fid=ncdf_open(fn)

if fid lt 0 then begin
	print,'no such file'
	stop
endif ; end file check

; get list of variables in file
varids=ncdf_varidsinq(fid)
nvar=n_elements(varids)
print,'number of variables',nvar
for i=0,nvar-1 do begin
	a=ncdf_varinq(fid,varids(i))
	varname=a.name
	vartype=a.datatype
	ncdf_varget,fid,varids(i),outval
;	print,'reading '+varname
	if i eq 0 then begin 
		outgroup=create_struct(varname,outval)
	endif else begin
		outgroup=create_struct(varname,outval,outgroup)
	endelse
	
endfor

end
