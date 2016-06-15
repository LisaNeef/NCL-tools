load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; This script computes the following TIL diagnostics for 
; ERA-Interim data on hybrid model levels.  
; this script writes two output files for each input file:  
;	1. Files with names like 'Posterior_TIL_XXXX.nc':
;		theta (potential temperature)
;		pressure  
;		brunt (Brunt-Vaisala frequency in s^-1, i.e. N)  
;		STIL - strength of the TIL, 
;			i.e. the maximum buoyancy frequency 
;			within a few hundred meters above the tropopause
;
; Looping over netcdf files:  
;  This script reads in a list of files given in file_list.
;
; File requirements: 
; + vertical levels must be in hybrid coordinates and go from the top of the
;	atmosphere to the surface  
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
	vdim = new(3,integer)
	vdim = dimsizes(fin->T(0,:,:,:))
	time = dimsizes(fin->T(:,0,0,0))

	;================================================================;
	; Initialize output files
	; Here the string 'TUV' is replaced with 'TIL' to create new filenames 
	; so here choose whatever string you want to replace in input filenames
	;================================================================;
	strs = str_split(filename_in,"/")
	filename_in_stripped = strs(dimsizes(strs)-1)
	filename_out = str_sub_str(filename_in_stripped, "TUV", "TIL")
	system( "rm " + filename_out)
	fout = addfile(filename_out,"c")  

	;========================================================================;
	; create coordinate variables for output files
	; TODO: do we need to define the vertical coordinat differently for the
	;  	TP-based one?  
	;========================================================================;
	dimNames = (/ "time","lat","lon","lev","nhym" /)
	dimSizes = (/ time,vdim(0),vdim(1),vdim(2),vdim(2) /)
	dimUnlim = (/ False,False,False,False,False /)
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

	;===================================================================
	; predefine the dimensionality of the variables to be written out
	;===================================================================

	filevardef(fout, "time" ,typeof(fin->time),getvardims(fin->time))
	filevardef(fout, "lat" ,typeof(fin->lat),getvardims(fin->lat))
	filevardef(fout, "lon" ,typeof(fin->lon),getvardims(fin->lon))
	filevardef(fout, "lev" ,typeof(fin->lev),getvardims(fin->lev))
	filevardef(fout, "hyam" ,typeof(fin->hyam),getvardims(fin->hyam))
	filevardef(fout, "hybm" ,typeof(fin->hybm),getvardims(fin->hybm))
	filevardef(fout, "theta" ,typeof(fin->T),getvardims(fin->T))
	filevardef(fout, "P" ,typeof(fin->T),getvardims(fin->T))
	filevardef(fout, "brunt" ,typeof(fin->T),getvardims(fin->T))
	filevardef(fout, "ptrop" ,typeof(fin->PS),getvardims(fin->PS))

	;===================================================================
	; Copy attributes associated with existing variables to output files
	;====================================================================

	filevarattdef(fout,"time" ,fin->time)                    ; copy time attributes
	filevarattdef(fout,"lat"  ,fin->lat)                     ; copy lat attributes
	filevarattdef(fout,"lon"  ,fin->lon)                     ; copy lon attributes
	filevarattdef(fout,"lev"  ,fin->lev)                     ; copy lev attributes
	filevarattdef(fout,"hyam"  ,fin->hyam)                     ; copy hyam attributes
	filevarattdef(fout,"hybm"  ,fin->hybm)                     ; copy hybm attributes
	if (fin->PS@units.eq."hPa".or.fin->PS@units.eq."millibars") then
		punit=0
	else
		punit=1
	end if 


	;===================================================================
	; Define attributes for new variable arrays 
	;====================================================================
	fout->theta@units 		= "K"
	fout->theta@long_name	= "Potential Temperature"
	fout->theta@_FillValue	= 9.96921e+36

	fout->P@units 		= fin->PS@units
	fout->P@long_name	= "Pressure"
	fout->P@_FillValue	= 9.96921e+36

	fout->brunt@units 		= "s^-2"
	fout->brunt@long_name	= "Brunt-Vaisala Frequency"
	fout->brunt@_FillValue	= 9.96921e+36

	fout->ptrop@units 	= fin->PS@units
	fout->ptrop@long_name	= "Lapse-rate tropopause pressure"
	fout->ptrop@_FillValue	= 9.96921e+36

	;===================================================================
	; output only the data values since the dimensionality and such have
	; been predefined. The "(/", "/)" syntax tells NCL to only output the
	; data values to the predefined locations on the file.
	;====================================================================

	fout->time   		= (/fin->time/)     
	fout->lat    		= (/fin->lat/)
	fout->lon    		= (/fin->lon/) 
	fout->lev    		= (/fin->lev/) 
	fout->hyam    		= (/fin->hyam/) 
	fout->hybm    		= (/fin->hybm/) 

	;====================================================================
	; Read in information about vertical levels common to each DART copy 
	; TODO: is this still necessary?  
	;====================================================================
	hyam = fin->hyam ; read from a file the mid-layer coef
	hybm = fin->hybm ; read from a file

	;====================================================================
	; loop over time and compute 
	;  + potential temperature  
	;  + brunt-vaisala frequency  
	;====================================================================

	; loop over time 
	do tt=0,time-1

		; loop over all lats and lons 
		do ii=0,vdim(0)-1
			do jj=0,vdim(1)-1
				; compute pressure at the given point 
				fout->P(tt,ii,jj,:)=hyam+hybm*fin->PS(tt,ii,jj)

				; vector of pressure at the given point 
				p = fout->P(tt,ii,jj,:)
					   
				; temperature 
				t = fin->T(tt,ii,jj,:)
	
				; potential temperature  
				; function pot_temp only available for NCL v6.2 and later
				;fout->theta(tt,ii,jj,:) = pot_temp(p, t, -1, False)  
				fout->theta(tt,ii,jj,:) = t*(p0/p)^(0.286)
				
				; geometric height in m
				z = 7000.0*log(p0/p)

				; calculate Brunt-Vaisala frequency and potential temp
				; NCL version 6.4 (not released yet) will have this function:
				; fout->brunt = brunt_vaisala_atm(fout->theta(tt,ii,jj,:), z, 0, 1)  

				; N^2 --> for now we have to do it manually
				theta=fout->theta(tt,ii,jj,:)
				dthetadz  = center_finite_diff(theta,z,False,0)
				g=9.80665
				fout->brunt(tt,ii,jj,:) = (g/theta)*dthetadz

				; also compute the tropopause pressure
				fout->ptrop(tt,ii,jj) = trop_wmo(p, t, punit, False)

			end do	; lon loop 
		end do 	; lat loop 
	end do 	; time loop 
	;====================================================================
	

	print("writing Netcdf file "+filename_out)

end do 

end  
