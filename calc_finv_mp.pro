PRO calc_finv_mp, num_func, func_indx, ret_nlev, htop, hbot, air_pres, $
		f_matrix, f_inv
; ---------------------------------------------
; PURPOSE: calculates a scene-dependent transformation matrix (F_matrix) and 
;   			its inverse the using Moore-Penrose pseudoinverse technique  
;  	The transformation matrix consists of the retrieval variable
;     	trapezoid state functions [e.g., 21 for water vapor, 9 for ozone]
;     	on the standard 100 pressure level retrieval grid, e.g., dimensions [9 x 100]
;  	This matrix is scene-dependent because we use only those functions and 
;     	pressure levels that are above Earth surface at a target scene.

; INPUT: 
;  Name          Description 
; ---------    ----------------------
; num_func      number of state functions above Earth surface
; func_indx     trapezoid state function hinge-points as reported in [ave_kern/*_func_indxs] 
; ret_nlev      number of retrieval pressure levels (air_pres) above Earth surface  
; htop       	 value in [ave_kern/*_func_hbot]
; hbot       	 value in [ave_kern/*_func_htop]
; air_pres      standard 100 level pressure grid [air_pres]/100 in hPa units

; OUTPUT: 
; Name			 Description
; --------     ---------------------
; f_matrix     [L,j], where L=retrieval levels, j=trapezoid levels, j<L
; 					It is the transformation matrix where each retrieval state function is on the standard 
;					   retrieval pressure grid
; f_inv        [j,L] and it is the pseudoinverse matrix of f_matrix
;
; DEPENDENCIES:
;  slb2fin.pro
;  compile this module if necessary on IDL command line: .run slb2fin
 
; --------------------------------------------------
; Original code by Chris D. Barnet and Eric S. Maddy
; chrisdbarnet@gmail.com
;
; Modified by Nadia Smith
; nadias@stcnet.com
; Science and Technology Corporation
; 6 July 2020
; ------------------------------------------------
; ------------------------------------------------
; Step 1: calculate the transformation matrix: f_matrix
; ------------------------------------------------
	
	ndim = num_func
	f_matrix = FLTARR(ret_nlev,ndim)
	
   FOR ifunc = 0, ndim - 1 DO BEGIN
;  Call slb2fin for one state function at a time setting the corresponding slbval = 1.0
		slbval = FLTARR(ndim)
   	slbval(ifunc) = 1.0
      fine = FLTARR(ret_nlev)
   	SLB2FIN, num_func, func_indx, slbval, htop, hbot, $
                   1100., air_pres, fine 
      f_matrix(*,ifunc) = fine; 
   ENDFOR  

; Subset f_matrix to remove trailing zeros

s=size(f_matrix)
nL=s(1); ret_nlev, max=100
nj=s(2); ak_nlev, max=30

for i=nL-1,0,-1 do if (n_elements(where(f_matrix(i,*) eq 0.0)) lt nj) then break
rpos=i
for i=nj-1,0,-1 do if (n_elements(where(f_matrix(*,i) eq 0.0)) lt nL) then break
cpos=i
f_matrix=f_matrix(0:rpos,0:cpos); [nL x nj] 

; ------------------------------------------------
; Step 2: calculate the inverse of f_matrix using 
;         the Moore-Penrose pseudoinverse method
; ------------------------------------------------
; This matrix, fftr, should have no zeros on diagonal, otherwise inversion fails
	fftr = MATRIX_MULTIPLY(DOUBLE(f_matrix),DOUBLE(f_matrix),/ATRANSPOSE); [nj x nj]
	status=1L
   finv1 = LA_INVERT(fftr,STATUS=status,/double)
	print,status; status should be zero to indicate successful inversion
   f_inv = MATRIX_MULTIPLY(finv1,DOUBLE(f_matrix),/BTRANSPOSE); [nj x nL]
 
END 
