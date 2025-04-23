pro prepare_clim_for_assimilation, fn,nfor,nscan,climco

; -----------------------------------------
; Routine to demonstrate how to prepare CLIMCAPS
; Level 2 retrievals for ingest into a data assimilation system
;
; There are four main steps:
; (1) Read retrieved profiles [molec/cm2] in 100 pressure layers
; (2) Adjust bottom layer according to surface pressure at [nfor,nscan]
; (3) Extract and expand AK matrices to 100 pressure layers
; (4) Generate a-priori climatology profile. 
;
; Step 4 is specifically for CO. The a-priori profiles for 
; temperature, water vapor and ozone are available in the Level 2 file
;
; Nadia Smith, nadias@stcnet.com
; 12 April 2023
; -----------------------------------------
 
; -----------------------------------------
; INPUT: 
; "fn" is a string that defines the filename and directory
; [nfor,nscan] stipulates the position in the 2-D Level 2 file
; 		where nfor is the field-of-regard index ranging [0..29]
;     and nscan is the scanline index ranging [0..44]
; Each CLIMCAPS Level 2 file has 30 FORs in each scanline
; and 45 scanlines in total. 
; 
; OUTPUT:
; climco is a data structure with the following elements
; 		climco.ret is the 100-layer profile retrieval [molec/cm2]
;		climco.prior is the 100-layer climatological a-priori 
; 				that represents CO at the [nfor,nscan] lat/lon in a dry 
; 				atmosphere [molec/cm2]
; 		climco.ak_matrix is the 100-layer averaging kernel matrix
;				derived from the 2D coarse-layer AK matrix
; 		climco.pressure is the 100-layer pressure array 
; 		climco.qc is the quality flag to apply, 0=success, ≥1=fail
; --------------------------------------------

; define constants
totfor = 30
totscan = 45
ret_nlev=100
fillval=9.96921e+36
nan=!Values.F_NAN
navog=6.02214199e+23; NIST Avogadro's number (molecules/mole)
mw_wv=18.0151; gm/mole water
mw_co=28.0104; gm/mole carbon monoxide
mw_d=28.9644; gm/mole dry air
g_std=980.664; acceleration of gravity, cm/s^2
eps=mw_wv/mw_d
co_eps=mw_co/mw_d;
ppb_scale=1.0e9
cm2_scale=1.0e4

; initialize the output structure
climco = create_struct('ret',fltarr(ret_nlev)); [molec/cm2]
climco = create_struct('prior',fltarr(ret_nlev),climco); [molec/cm2]
climco = create_struct('ak_matrix',fltarr(ret_nlev,ret_nlev),climco); [100 x 100]
climco = create_struct('pressure',fltarr(ret_nlev),climco); [100 x 100]
climco = create_struct('qc',[1],climco)

; -------------------------------------------
; STEP 1: read the Level 2 file and extract all relevant values
; --------------------------------------------
; Read main branc of CLIMCAPS netcdf file 
ncread_climcaps_main,fn,clim

; Read AK group
ncread_climcaps_grps,fn,'ave_kern',akgrp

; Read Aux group
ncread_climcaps_grps,fn,'aux',auxgrp

; Read Mol_lay group where all the retrieved gases are stored
ncread_climcaps_grps,fn,'mol_lay',molgrp

; Extract CO relevant fields at [nfor,nscan]
surf_pres=auxgrp.prior_surf_pres(nfor,nscan)/100.; surface pressure [hPa]
htop = akgrp.co_func_htop
hbot = akgrp.co_func_hbot
ak_coarse = reform(akgrp.co_ave_kern(*,*,nfor,nscan)); [9 x 9]
ak_pidx = akgrp.co_func_indxs;
ak_nlev = n_elements(ak_pidx)
; index where AKs intersect Earth surface
ak_nsurf = akgrp.co_func_last_indx(nfor,nscan)

; clim.air_pres define the pressure level boundaries of the pressure 
; 		layers (clim.air_pres_lay) on which the trace gases are retrieved
; All CLIMCAPS gas retrievals (everything except air_temp and co2_vmr) should be associated
; 		with clim.air_pres_lay
ret_pres = clim.air_pres/100.; Pressure levels [hPa]
ret_pres_lay = clim.air_pres_l/ay/100.; Pressure layers [hPa]
climco.pressure = ret_pres_lay

; index where ret_pres_lay intersect the Earth surface
nsurf = clim.air_pres_nsurf(nfor,nscan)

; Retrieval profiles
; Level 2 file report gases in [molec/m2]
ret_prof = reform(molgrp.co_mol_lay(*,nfor,nscan))
ret_prof(where(ret_prof eq fillval)) = nan
ret_prof = ret_prof/cm2_scale; convert to [molec/cm2]

; Read h2o_vap retrieval to convert prior profile to molec/cm2
wv_prof = reform(molgrp.h2o_vap_mol_lay(*,nfor,nscan))
wv_prof(where(wv_prof eq fillval)) = nan
wv_prof = wv_prof/cm2_scale; convert to [molec/cm2]

plat = clim.lat(nfor,nscan)
; use profile if qc=0, otherwise, reject
qc = auxgrp.ispare_2(nfor,nscan)
climco.qc = qc

; extract datestr from filename
pos = strpos(fn,'CRIMSS')
ndum=strlen('CRIMSS')+1
datestr = strmid(fn,pos+ndum,8); YYYYMMDD

; -------------------------------------------
; (2) Adjust bottom layer according to surface pressure at [nfor,nscan]
; -------------------------------------------
; Retrieval
L = nsurf-1 ; adjust for IDL indexing starting at 0
blmult = (surf_pres - ret_pres(L-1))/(ret_pres(L)-ret_pres(L-1))
ret_prof(L) = ret_prof(L)*blmult
climco.ret = ret_prof

; Averaging Kernel
; surf_pres     surface pressure [hPa]
; ret_pres      Standard RTA pressure level grid
; ret_nlev    number of levels (RTA grid)
; ak_nlev    number of coarse layers to be adjusted for topography
; ak_pres     coarse layers to be adjusted for topography 
; OUTPUT
; ak_nlev_scene    number of coarse level boundaries at scene psurf
; ak_pidx_scene    coarse level boundaries at scene psurf
changlev, surf_pres, ret_pres, ret_nlev, ak_nlev, ak_pidx, $
	ak_nlev_scene, ak_pidx_scene

; --------------------------------------------
; (3) Expand AK matrix to 100 pressure layers
; --------------------------------------------
; Now compute the transformation matrix and its pseudo inverse 
; Equation 11 of Maddy & Barnet, IEEE TGRS, 2008 (MB08) 
; func_matrix is the [nj x nL] matrix of trapezoids
; func_inv is the inverse of func_matrix [nL x nj]
; Notation in Maddy & Barnet corresponds to the code here as follows:
; nL = ret_nlev_scene
; nj = ak_nlev_scene
; F = func_matrix
; F+ = func_inv

calc_finv_mp,ak_nlev_scene,ak_pidx_scene,ret_nlev, htop, hbot,$
	ret_pres,func_matrix,func_inv

help,func_matrix
help,func_inv
; Effective 100-layer AK matrix: F*AK*F+
ak_expand = func_matrix#ak_coarse#func_inv; [nL x nL]
climco.ak_matrix(0:nsurf-1,0:nsurf-1) = ak_expand

; --------------------------------------------
; (4) Generate a-priori climatology profile
; --------------------------------------------
; initialize output vector
co_prior=fltarr(ret_nlev)*nan

; get CO climatology profile in units [ppb]
get_mopitt_v4clim,datestr,plat,co_prior_vmr

; adjust bottom layer that intersects Earth surface
co_prior_vmr(nsurf-1) = co_prior_vmr(nsurf-1) * blmult

; convert CO mixing ratio to dry column density
navog=6.02214199e+23; NIST Avogadro's number (molecules/mole)
mw_wv=18.0151; gm/mole water
mw_co=28.0104; gm/mole carbon monoxide
mw_d=28.9644; gm/mole dry air
g_std=980.664; acceleration of gravity, cm/s^2
eps=mw_wv/mw_d
co_eps=mw_co/mw_d;

; cd_d() is the layer column density for dry air
; ----------------------------------------------
; cd_t(L) is the molecules/cm^2 of air in that layer
;         cd_d(L) is the molecules/cm^2 of dry air in that layer
;         svp(L)  is the saturation vapor pressure of this layer
cd_t=fltarr(ret_nlev)*nan
cd_d=fltarr(ret_nlev)*nan

; define first layers before looping over remaining layers
cd_t(0)=1000.*ret_pres(0)*navog/(mw_d*g_std)
cd_d(0)=cd_t(0) - eps*wv_prof(0)

for j=1,nsurf-1 do begin
   cd_t(j)=1000.*(ret_pres(j)-ret_pres(j-1))*navog/(mw_d*g_std)
   cd_d(j)=cd_t(j)-eps*wv_prof(j)
endfor

; adjust bottom layer to surface pressure by multiplying with blmult
cd_d(nsurf-1)=blmult*cd_d(nsurf-1)

for j=0,nsurf-1 do begin
	co_prior(j)=co_prior_vmr(j)*(cd_d(j)/ppb_scale)
endfor

climco.prior = co_prior; [molec/cm2]

end
