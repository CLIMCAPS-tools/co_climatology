; same function as the FORTRAN version.  Builds trapezoids.
; numslab = number of levels in Lslab
; Lslab   = index values for trapezoid hinge points, starting at 1
; tmpslab = amplitude of trapezoids
; usehalftop = 0 (.false.) for a trapezoid, 1 (.TRUE.) is a wedge
; Psurf = pressure at surface
; Pres  = pressure values at index points
; tmpfine = amplitude on fine levels
; tmpfile(L) = tmpslab_i*F_i(L)


pro slb2fin, numslab, Lslab, tmpslab, usehalftop, usehalfbot,$
    Psurf, Pres, tmpfine
   
     tslab = fltarr(numslab)

;    ------------------------------------------------------
;    Interpolate temperatures ( linear in log of pressure )
;    ------------------------------------------------------

      if(usehalftop gt 0) then begin
        tslab(0) = 0.5 * tmpslab(0)
      endif else begin
        tslab(0) = tmpslab(0)
      endelse

      for n = 1, numslab-2 do begin
        tslab(n) = 0.5 * ( tmpslab(n) + tmpslab(n-1) )
      endfor

      if ( usehalfbot gt 0) then begin
        tslab(numslab-1) = 0.5 * tmpslab(numslab-2)
      endif else begin
        tslab(numslab-1) = tmpslab(numslab-2)
      endelse

      ldn          = lslab(0) - 1
      ztmpdn       = tslab(0)
      zprlndn      = alog (pres(ldn) )

      for n = 0, numslab-2 do begin
        lup        = ldn
        ztmpup     = ztmpdn
        zprlnup    = zprlndn
        ldn        = lslab(n+1) - 1
        ztmpdn     = tslab(n+1)
        zprlndn    = alog(Pres(ldn))
        if (pres(ldn) gt  psurf ) then zprlndn = alog(psurf)
        zslope     = ( ztmpdn - ztmpup ) / ( zprlndn - zprlnup )

        for L = lup, ldn-1 do begin
          tmpfine(L) = ztmpup + zslope * ( alog(Pres(L)) - zprlnup )
        endfor

        lup = L
      endfor

      tmpfine(ldn) = ztmpdn

end

