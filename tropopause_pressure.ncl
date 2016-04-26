load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; Compute the pressure of the tropopause using the WMO criterion 
; coded in trop_wmo.ncl, for all lats and lons in a WACCM/CESM 
; history file  
;
; File requirements: 
; + vertical levels must be in hybrid coordinates and go from the top of the
;	atmosphere to the surface  
; + currently I've designed this to read 'h1' history files, but 
;	just modify the filename_out definition below to deal with something 
;	else  
;
; Lisa Neef  
; GEOMAR Helmholtz Centre for Ocean Research Kiel  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


begin 


;=================================================;
; Read all file names in directory
;=================================================;
diri = "./" ; input directory
fili = systemfunc("cd "+diri+" ; ls *nc")
nfili = dimsizes(fili) 

;=================================================;
; Loop over all available files 
;=================================================;
do nf=0,nfili-1
	filename_in=fili(nf)
	fin = addfile(diri+filename_in, "r")
	print("reading file "+filename_in)

	;=================================================;
	; Open input file and read in array dimensions 
	;=================================================;
	; read in netcdf  
	fin = addfile(filename_in,"r")  

	; get array dimensions for a 2d array (from surface pressure) 
	vdim = new(2,integer)
	vdim = dimsizes(fin->PS(0,:,:))
	time = dimsizes(fin->PS(:,0,0))

	;=================================================;
	; Initialize output file 
	;=================================================;
	filename_out = str_sub_str(filename_in, "h1", "tropopause")
	system( "rm " + filename_out)
	fout = addfile(filename_out,"c")  

	;===================================================================
	; create coordinate variables for output file 
	;===================================================================
	dimNames = (/ "time","lat","lon" /)
	dimSizes = dimsizes(fin->PS(:,:,:))
	dimUnlim = (/ False,False,False /)
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

	;===================================================================
	; predefine the the dimensionality of the variables to be written out
	;===================================================================

	filevardef(fout, "time" ,typeof(fin->time),getvardims(fin->time))
	filevardef(fout, "lat" ,typeof(fin->lat),getvardims(fin->lat))
	filevardef(fout, "lon" ,typeof(fin->lon),getvardims(fin->lon))
	filevardef(fout, "ptrop" ,typeof(fin->PS),getvardims(fin->PS))

	;===================================================================
	; Copy attributes associated with each variable to the file
	; All attributes associated with each variable will be copied.
	;====================================================================

	filevarattdef(fout,"time" ,fin->time)                    ; copy time attributes
	filevarattdef(fout,"lat"  ,fin->lat)                     ; copy lat attributes
	filevarattdef(fout,"lon"  ,fin->lon)                     ; copy lon attributes

	;===================================================================
	; output only the data values since the dimensionality and such have
	; been predefined. The "(/", "/)" syntax tells NCL to only output the
	; data values to the predefined locations on the file.
	;====================================================================

	fout->time   = (/fin->time/)     
	fout->lat    = (/fin->lat/)
	fout->lon    = (/fin->lon/) 

	;====================================================================
	; Compute pressures from hybrid level definition  
	;====================================================================

	hyam = fin->hyam ; read from a file the mid-layer coef
	hybm = fin->hybm ; read from a file
	ps   = fin->PS   ; surface pressure 
	p0   = fin->P0	 ; reference pressure 
	pi   = pres_hybrid_ccm(ps,p0,hyam,hybm)  

	if (fin->P0@units.eq."hPa".or.fin->P0@units.eq."millibars") then
		punit=0
	else
		punit=1
	end if 

	;====================================================================
	; loop over time 
	do tt=0,time-1
		; loop over all lats and lons 
		do ii=0,vdim(0)-1
			do jj=0,vdim(1)-1
				; vector of pressure at the given point 
				p = pi(tt,:,ii,jj)
					   
				; TEMPERATURE (C)   
				t = fin->T(tt,:,ii,jj)

				; compute the pressure of the tropopause  
				fout->ptrop(tt,ii,jj) = trop_wmo(p, t, punit, False)     
			end do
		end do 
	end do 
	;====================================================================
	
	print("writing Netcdf file "+filename_out)

end do 

end  
