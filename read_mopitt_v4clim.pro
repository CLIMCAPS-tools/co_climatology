pro read_mopitt_v4clim,array

; CLIMCAPS V2.1 CO a priori estimate

; Define constants
maxhm=2; maximum two hemispheres, 1=North, 2=South
maxmon=12; maximum number of months in a year
maxplev=100; maximum number of pressure layers
totrec=float(maxhm*maxmon*100)
totparam=4; total number of parameters in a line

; Read text file
fn='mopitt_v4clim_nhsh_20180605.txt'

openr, lun, fn, /get_lun

header=''
; read header
readf,lun,header

array=[totparam,totrec]
hm=long(0); hemisphere
m=long(0); month
p=double(0.0); pressure
co=double(0.0); co mixing ratio
help,array
count=0
while not eof(lun) do begin
	readf,lun,hm,m,p,co,format='(I12,I12,D13,G20.13E4)'
	line=[hm,m,p,co]
	print,count,line
	array=[line,count]
	count=count+1
endwhile
free_lun,lun
end
