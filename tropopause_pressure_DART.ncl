load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; Compute the pressure of the tropopause using the WMO criterion 
; coded in trop_wmo.ncl, for all lats and lons in DART output file, 
; i.e. files that have names like 'Posterior_Diag_XXXX.nc', etc. 
; This script replaces them with names like 'Posterior_tropopause_XXXX.nc'  
;
;
; File requirements: 
; + vertical levels must be in hybrid coordinates and go from the top of the
;	atmosphere to the surface  
; + the files must have time as the first dimension and 'copy' as the 
;	second. 
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
	vdim = dimsizes(fin->PS(0,0,:,:))
	time = dimsizes(fin->PS(:,0,0,0))
	copy = dimsizes(fin->PS(0,:,0,0))
	metadatalength = dimsizes(fin->CopyMetaData(0,:))

	;=================================================;
	; Initialize output file 
	;=================================================;
	filename_out = str_sub_str(filename_in, "Diag", "tropopause")
	system( "rm " + filename_out)
	fout = addfile(filename_out,"c")  

	;===================================================================
	; create coordinate variables for output file 
	;===================================================================
	dimNames = (/ "time","copy","lat","lon","metadatalength" /)
	dimSizes = (/ time,copy,vdim(0),vdim(1),metadatalength /)
	dimUnlim = (/ False,False,False,False,False /)
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

	;===================================================================
	; predefine the the dimensionality of the variables to be written out
	;===================================================================

	filevardef(fout, "time" ,typeof(fin->time),getvardims(fin->time))
	filevardef(fout, "copy" ,typeof(fin->copy),getvardims(fin->copy))
	filevardef(fout, "CopyMetaData" ,typeof(fin->CopyMetaData),getvardims(fin->CopyMetaData))
	filevardef(fout, "lat" ,typeof(fin->lat),getvardims(fin->lat))
	filevardef(fout, "lon" ,typeof(fin->lon),getvardims(fin->lon))
	filevardef(fout, "ptrop" ,typeof(fin->PS),getvardims(fin->PS))

	;===================================================================
	; Copy attributes associated with each variable to the file
	; All attributes associated with each variable will be copied.
	;====================================================================

	filevarattdef(fout,"time" ,fin->time)                    ; copy time attributes
	filevarattdef(fout,"copy" ,fin->copy)                    ; copy copy attributes
	filevarattdef(fout,"CopyMetaData" ,fin->CopyMetaData)    ; copy CopymetaData attributes
	filevarattdef(fout,"lat"  ,fin->lat)                     ; copy lat attributes
	filevarattdef(fout,"lon"  ,fin->lon)                     ; copy lon attributes

	;===================================================================
	; Define attributes for tropopause pressure. 
	;====================================================================
	fout->ptrop@units 	= fin->PS@units
	fout->ptrop@long_name	= "Lapse-rate tropopause pressure"
	fout->ptrop@_FillValue	= 9.96921e+36

	;===================================================================
	; output only the data values since the dimensionality and such have
	; been predefined. The "(/", "/)" syntax tells NCL to only output the
	; data values to the predefined locations on the file.
	;====================================================================

	fout->time   = (/fin->time/)     
	fout->copy   = (/fin->copy/)     
	fout->CopyMetaData   = (/fin->CopyMetaData/)     
	fout->lat    = (/fin->lat/)
	fout->lon    = (/fin->lon/) 


	;====================================================================
	; Read in information about vertical levels common to each DART copy 
	;====================================================================
	hyam = fin->hyam ; read from a file the mid-layer coef
	hybm = fin->hybm ; read from a file
	if (fin->P0@units.eq."hPa".or.fin->P0@units.eq."millibars") then
		punit=0
	else
		punit=1
	end if 

	;====================================================================
	; loop over time 
	do tt=0,time-1
		; loop over copy
		;do cc=0,copy-1
		do cc=0,1

			;====================================================================
			; Compute pressures from hybrid level definition  
			;====================================================================
			ps   = fin->PS(:,cc,:,:)   ; surface pressure 
			p0   = fin->P0(0)	; reference pressure 
			pi   = pres_hybrid_ccm(ps,p0,hyam,hybm)  

			; loop over all lats and lons 
			do ii=0,vdim(0)-1
				do jj=0,vdim(1)-1
					; vector of pressure at the given point 
					p = pi(tt,:,ii,jj)
						   
					; temperature 
					t = fin->T(tt,cc,ii,jj,:)

					; compute the pressure of the tropopause  
					fout->ptrop(tt,cc,ii,jj) = trop_wmo(p, t, punit, False)     

				end do	; lon loop 
			end do 	; lat loop 
		end do 	; copy loop 
	end do 	; time loop 
	;====================================================================
	
	print("writing Netcdf file "+filename_out)

end do 

end  
