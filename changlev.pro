;                                                                                         
; NAME: calc_finv                                                                 
; PURPOSE: adjust the coarse lev grid accounting for topography
;
; INPUT:
;                                        
;  Name          Description
; ------    ----------------------
; Psurf     surface pressure
; Pstd      Standard RTA pressure grid
; numlev    number of levels (RTA grid)
; nlev_0    number of coarse layers to be adjusted for topography
; lev_0     coarse layers to be adjusted for topography 

; OUTPUT
; nlev_1    number of coarse levels
; lev_1     coarse level boundaries 

pro changlev, Psurf, Pstd, numlev, nlev_0, lev_0, nlev_1, lev_1

thinlev = 50.0

numlevold = numlev
numlevnew = numlev

lsurface, numlevnew, Pstd, Psurf, 100., 1100.
numlev = numlevnew


nlev = nlev_0
lev  = lev_0 

for L = 0, nlev_0-1 do begin

  if ( lev(L) gt numlevnew ) then begin
    if (lev(l-1) eq numlevnew) then begin
      nlev = nlev - 1
    endif else begin
      if ( Pstd(numlevnew) - Pstd(lev(l-1)) lt thinlev ) then begin
        nlev = nlev - 1
        lev(L-1) = numlevnew
      endif
    endelse
    lev(L) = numlevnew
  endif

endfor

nlev_1 = nlev
lev_1  = lev(0:nlev-1)

end
