load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; Compute the pressure of the tropopause using the WMO criterion 
; coded in trop_wmo.ncl, for all lats and lons in ERA-Interim 
; data on constant pressure levels 
; 
; How to use:  
; 	+ create a file called file_list that is basically a list of all the files 
; 	  that you want to compute tropopause pressure for
;	+ run this program
; File requirements:  
; 	+ data must be on constant pressure levels that go from the top of the 
;		atmo to the surface  
;	+ file must contain a temperature array named t, with shape time-lev-lat-lon 
;		or time-lev-lon-lat and temps given in Kelvin 
;
; Lisa Neef  
; GEOMAR Helmholtz Centre for Ocean Research Kiel  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


begin 

;=================================================;
; Read the files in the  list 
;=================================================;
fili = asciiread("file_list",-1,"string")

;=================================================;
; Loop over all available files 
;=================================================;
nfili = dimsizes(fili)
do nf=0,nfili-1
	filename_in=fili(nf)
	fin = addfile(filename_in, "r")
	print("reading file "+filename_in)

	;=================================================;
	; Open input file and read in array dimensions 
	;=================================================;
	; read in netcdf  
	fin = addfile(filename_in,"r")  

	; get array dimensions for a 2d array 
	vdim = new(2,integer)
	vdim = dimsizes(fin->t(0,0,:,:))
	time = dimsizes(fin->t(:,0,0,0))

	;=================================================;
	; Initialize output file 
	;=================================================;
	filename_out = str_sub_str(filename_in, "t", "ptrop")
	system( "rm " + filename_out)
	fout = addfile(filename_out,"c")  

	;===================================================================
	; create coordinate variables for output file 
	;===================================================================
	dimNames = (/ "time","latitude","longitude" /)
	dimSizes = dimsizes(fin->t(:,0,:,:))
	dimUnlim = (/ False,False,False /)
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

	;===================================================================
	; predefine the the dimensionality of the variables to be written out
	;===================================================================

	filevardef(fout, "time" ,typeof(fin->time),getvardims(fin->time))
	filevardef(fout, "latitude" ,typeof(fin->latitude),getvardims(fin->latitude))
	filevardef(fout, "longitude" ,typeof(fin->longitude),getvardims(fin->longitude))
	filevardef(fout, "ptrop" ,"double",getvardims(fin->t(:,0,:,:)))

	;===================================================================
	; Copy attributes associated with each variable to the file
	; All attributes associated with each variable will be copied.
	;====================================================================

	filevarattdef(fout,"time" ,fin->time)                    ; copy time attributes
	filevarattdef(fout,"latitude"  ,fin->latitude)                     ; copy lat attributes
	filevarattdef(fout,"longitude"  ,fin->longitude)                     ; copy lon attributes

	;===================================================================
	; output only the data values since the dimensionality and such have
	; been predefined. The "(/", "/)" syntax tells NCL to only output the
	; data values to the predefined locations on the file.
	;====================================================================

	fout->time   = (/fin->time/)     
	fout->latitude    = (/fin->latitude/)
	fout->longitude    = (/fin->longitude/) 


	;====================================================================
	;Load pressure levels from input file 
	;====================================================================

	p = fin->level	
	if (p@units.eq."hPa".or.p@units.eq."millibars") then
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
					   
				; TEMPERATURE (K)   
				tflt = short2flt(fin->t)
				t = tflt(tt,:,ii,jj)

				; compute the pressure of the tropopause  
				fout->ptrop(tt,ii,jj) = trop_wmo(p, t, punit, False)     
			end do
		end do 
	end do 
	;====================================================================
	
	print("writing Netcdf file "+filename_out)

end do 

end  
