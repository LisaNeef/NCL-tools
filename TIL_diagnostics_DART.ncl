load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; This script computes the following TIL diagnostics for 
; DART output files, 
; i.e. files that have names like 'Posterior_Diag_XXXX.nc'.
; this script writes two output files for each input file:  
;	1. Files with names like 'Posterior_TIL_XXXX.nc':
;		theta (potential temperature)
;		pressure  
;		brunt (Brunt-Vaisala frequency in s^-1, i.e. N)  
;		STIL - strength of the TIL, 
;			i.e. the maximum buoyancy frequency 
;			within a few hundred meters above the tropopause
;	2. Files with names like 'Posterior_TIL_TPbased_XXXX.nc'
;		brunt (N relative to the tropopause)  
;
;
; Looping over netcdf files:  
;  This script reads in a list of files given in file_list.
;
; File requirements: 
; + vertical levels must be in hybrid coordinates and go from the top of the
;	atmosphere to the surface  
; + the files must have time as the first dimension and 'copy' as the 
;	second. 
;
; A note on comparing DART copies:  
; There are a few nonlinear computatons here: 
; 	- The computation `trop_wmo` selects a local minimum. 
;	- The N2max criterion (obvi) selects a local maximum  
; The present script applies these computations to individual ensembe members as well as the ensemble mean, 
; but applying them to the ensemble mean (i.e. a min of a mean) will yield a different answer 
; than taking the mean of the WMO criterion applied to each ensemble member. 
; When visualizing output, it makes more sense to look at these diagnostics for individual ensemble 
; members and then compute the mean of that. 



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

	; get array dimensions for a 2d array (from surface pressure) 
	vdim = new(3,integer)
	vdim = dimsizes(fin->T(0,0,:,:,:))
	time = dimsizes(fin->T(:,0,0,0,0))
	copy = dimsizes(fin->T(0,:,0,0,0))

	; also get information on the "copy" array  
	metadatalength = dimsizes(fin->CopyMetaData(0,:))
	CopyMetaData = fin->CopyMetaData

	; identify the copy indices where computing TIL diagnostics makes sense,
	; i.e. ensemble members and the mean, but not the ensemble spread 
	; or inflation factors 
	; This saves time and avoids contaminating the ptrop array with 
	; ridiculous values 
	ensemble_member_copies = str_match_ind_ic(CopyMetaData,"ensemble member")
	ensemble_mean_copies = str_match_ind_ic(CopyMetaData,"ensemble mean")
	selected_copies = array_append_record (ensemble_member_copies, ensemble_mean_copies, 0)
	nC = dimsizes(selected_copies) 

	;================================================================;
	; Initialize output files
	;================================================================;
	strs = str_split(filename_in,"/")
	filename_in_stripped = strs(dimsizes(strs)-1)
	filename_out = str_sub_str(filename_in_stripped, "Diag", "TIL")
	system( "rm " + filename_out)
	fout = addfile(filename_out,"c")  

	;========================================================================;
	; create coordinate variables for output files
	; TODO: do we need to define the vertical coordinat differently for the
	;  	TP-based one?  
	;========================================================================;
	dimNames = (/ "time","copy","lat","lon","lev","metadatalength" /)
	dimSizes = (/ time,copy,vdim(0),vdim(1),vdim(2),metadatalength /)
	dimUnlim = (/ False,False,False,False,False,False /)
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

	;===================================================================
	; predefine the dimensionality of the variables to be written out
	;===================================================================

	filevardef(fout, "time" ,typeof(fin->time),getvardims(fin->time))
	filevardef(fout, "copy" ,typeof(fin->copy),getvardims(fin->copy))
	filevardef(fout, "CopyMetaData" ,typeof(fin->CopyMetaData),getvardims(fin->CopyMetaData))
	filevardef(fout, "lat" ,typeof(fin->lat),getvardims(fin->lat))
	filevardef(fout, "lon" ,typeof(fin->lon),getvardims(fin->lon))
	filevardef(fout, "lev" ,typeof(fin->lev),getvardims(fin->lev))
	filevardef(fout, "hyam" ,typeof(fin->hyam),getvardims(fin->hyam))
	filevardef(fout, "hybm" ,typeof(fin->hybm),getvardims(fin->hybm))
	filevardef(fout, "P0" ,typeof(fin->P0),getvardims(fin->P0))
	filevardef(fout, "theta" ,typeof(fin->T),getvardims(fin->T))
	filevardef(fout, "P" ,typeof(fin->T),getvardims(fin->T))
	filevardef(fout, "brunt" ,typeof(fin->T),getvardims(fin->T))

	;===================================================================
	; Copy attributes associated with existing variables to output files
	;====================================================================

	filevarattdef(fout,"time" ,fin->time)                    ; copy time attributes
	filevarattdef(fout,"copy" ,fin->copy)                    ; copy copy attributes
	filevarattdef(fout,"CopyMetaData" ,fin->CopyMetaData)    ; copy CopymetaData attributes
	filevarattdef(fout,"lat"  ,fin->lat)                     ; copy lat attributes
	filevarattdef(fout,"lon"  ,fin->lon)                     ; copy lon attributes
	filevarattdef(fout,"lev"  ,fin->lev)                     ; copy lev attributes
	filevarattdef(fout,"hyam"  ,fin->hyam)                     ; copy hyam attributes
	filevarattdef(fout,"hybm"  ,fin->hybm)                     ; copy hybm attributes
	filevarattdef(fout,"P0"  ,fin->P0)                     ; copy P0 attributes


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

	;===================================================================
	; output only the data values since the dimensionality and such have
	; been predefined. The "(/", "/)" syntax tells NCL to only output the
	; data values to the predefined locations on the file.
	;====================================================================

	fout->time   		= (/fin->time/)     
	fout->copy   		= (/fin->copy/)     
	fout->CopyMetaData   	= (/fin->CopyMetaData/)     
	fout->lat    		= (/fin->lat/)
	fout->lon    		= (/fin->lon/) 
	fout->lev    		= (/fin->lev/) 
	fout->hyam    		= (/fin->hyam/) 
	fout->hybm    		= (/fin->hybm/) 
	fout->P0    		= (/fin->P0/) 

	;====================================================================
	; Read in information about vertical levels common to each DART copy 
	; TODO: is this still necessary?  
	;====================================================================
	hyam = fin->hyam ; read from a file the mid-layer coef
	hybm = fin->hybm ; read from a file

	;====================================================================
	; loop over copy and compute 3d pressure array 
	;====================================================================
	do cc=0,copy-1
		; Compute pressures from hybrid level definition  
		ps   = fin->PS(:,cc,:,:)   ; surface pressure 
		p0   = tofloat(fin->P0(0))	; reference pressure 
		pi = pres_hybrid_ccm(ps,p0,hyam,hybm)  
		; the DART output arrays are shaped like time-copy-lat-lon-lev 
		; whereas the above function returns it like tie-lev-lat-lon,
		; so have to loop over lat and lon and sort everything 
		do ii=0,vdim(0)-1
			do jj=0,vdim(1)-1
				fout->P(:,cc,ii,jj,:) = (/pi(:,:,ii,jj)/)
			end do
		end do
	end do

	;====================================================================
	; loop over selected copies and time, and compute 
	;  + potential temperature  
	;  + brunt-vaisala frequency  
	;====================================================================
	; loop over copy
	do icopy=0,nC-1
		cc = selected_copies(icopy)

		; loop over time 
		do tt=0,time-1

			; loop over all lats and lons 
			do ii=0,vdim(0)-1
				do jj=0,vdim(1)-1
					; vector of pressure at the given point 
					p = fout->P(tt,cc,ii,jj,:)
						   
					; temperature 
					t = fin->T(tt,cc,ii,jj,:)
		
					; potential temperature  
					; function pot_temp only available for NCL v6.2 and later
					;fout->theta(tt,cc,ii,jj,:) = pot_temp(p, t, -1, False)  
					fout->theta(tt,cc,ii,jj,:) = t*(p0/p)^(0.286)
					
					; geometric height in m
					z = 7000.0*log(p0/p)

					; calculate Brunt-Vaisala frequency and potential temp
					; NCL version 6.4 (not released yet) will this function:
					; fout->brunt = brunt_vaisala_atm(fout->theta(tt,cc,ii,jj,:), z, 0, 1)  

					; for now we have to do it manually
					theta=fout->theta(tt,cc,ii,jj,:)
					dthetadz  = center_finite_diff(theta,z,False,0)
					g=9.80665
					fout->brunt(tt,cc,ii,jj,:) = (g/theta)*dthetadz

				end do	; lon loop 
			end do 	; lat loop 
		end do 	; copy loop 
	end do 	; time loop 
	;====================================================================
	

	print("writing Netcdf file "+filename_out)

end do 

end  
