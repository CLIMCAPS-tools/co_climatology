pro ncread_climcaps_grps,fn,gn,outgroup

; Nadia Smith
; nadias@stcnet.com
; Aug 2019
; Routine to read NASA L2 CLIMCAPS format
; hierarchical netcdf with groups of variables
;
; input:
; - fn: filename with full path
; - gn: name of netcdf subgroup, 'mol_lay', 'ave_kern' or 'aux'
; output: 
; - outgroup: structure containing all the variables in group "gn"
; ----------------------------------------

; open ncdf file
fid=ncdf_open(fn)

if fid lt 0 then begin
	print,'no such file'
	stop
endif ; end file check

pid=ncdf_groupsinq(fid)
ngrp=n_elements(pid)

for i=0,ngrp-1 do begin
	gname=ncdf_groupname(pid(i))
	if strcmp(gname,gn) eq 1 then gid=pid(i)
endfor

fullgname=ncdf_fullgroupname(gid)
; get list of variables in group
varids=ncdf_varidsinq(gid)
nvar=n_elements(varids)
for i=0,nvar-1 do begin
	a=ncdf_varinq(gid,varids(i))
	varname=a.name
	vartype=a.datatype
	ncdf_varget,gid,varids(i),outval
	;print,'reading '+varname
	if i eq 0 then begin 
		outgroup=create_struct(varname,outval)
	endif else begin
		outgroup=create_struct(varname,outval,outgroup)
	endelse
	
endfor

end
