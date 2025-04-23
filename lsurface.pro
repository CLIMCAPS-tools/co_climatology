; input:
;   numlev = # of levels (e.g., 100)
;   pres = Pobs(0:numlev-1)
;   Psurf = surface pressure
;   Plow = lower limit of pressure for numlev
;
; output:
;   numlev = bottom level index in Pobs(1) system.
;     NOTE: for IDL numlev-1 is the bottom level


pro lsurface, numlev, pres, Psurf, Plow, Phigh
      phigh = 1100.

      if ( psurf gt plow and psurf le phigh ) then begin

         for L = numlev-1, 0, -1 do begin
            if ( psurf ge pres(L-1)+5.0 ) then begin
               lsurface = L
               goto, finish
            endif
         endfor
      endif

      if(psurf le plow) then begin
        lsurface = 0
        for L = 0, numlev-1 do begin
          if(psurf lt pres(L)) then lsurface = L
        endfor
      endif else begin
        lsurface = numlev-1
        for L = numlev-1,0,-1 do begin
          if(psurf gt pres(L)+5.0) then lsurface = L
        endfor
      endelse

  f100 = "('lsurface: ',f7.2,' <= (psurf=',f7.2,') <= ',f7.2,' pres(',i3,')=',f7.2)"
     print, plow, psurf, phigh, lsurface, pres(lsurface), format=f100

finish: 
   numlev=  lsurface + 1

end
