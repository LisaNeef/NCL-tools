load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; Compute the pressure of the tropopause using the WMO criterion 
; coded in trop_wmo.ncl, for all lats and lons in DART output file, 
; i.e. files that have names like 'Posterior_Diag_XXXX.nc', etc. 
; This script replaces them with names like 'Posterior_tropopause_XXXX.nc'  
;
; Looping over netcdf files:  
;  There are two options for looping: either go over all DART diagnostic 
;  files in the working directory, or read in a list of files given in file_list.
;  The latter is good if you have a lot of DART output sitting in a different 
;  directory and don't want to copy it all over. 
;  Uncomment the desired option below. 
;
;
; File requirements: 
; + vertical levels must be in hybrid coordinates and go from the top of the
;	atmosphere to the surface  
; + the files must have time as the first dimension and 'copy' as the 
;	second. 
;
; A note on comparing DART copies:  
; The computation `trop_wmo` is nonlinear -- it selects a local minimum. 
; The present script applies this computation to individual ensembe members as well as the ensemble mean, 
; but applying the WMO criterion to the ensemble mean (i.e. a min of a mean) will yield a different answer 
; than taking the mean of the WMO criterion applied to each ensemble member. 
; When visualizing output, it makes more sense to look at the tropopause pressure of individual ensemble 
; members and then compute the mean of that. 



; Lisa Neef  
; GEOMAR Helmholtz Centre for Ocean Research Kiel  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


begin 


;=================================================;
; Read all file names in directory
;=================================================;
;diri = "./" ; input directory
;fili = systemfunc("cd "+diri+" ; ls *Diag*nc")

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
	;fin = addfile(diri+filename_in, "r")
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

	; identify the copy indices where computing tropopause pressure makese sense,
	; i.e. ensemble members and the mean, but not the ensemble spread 
	; or inflation factors 
	; This saves time and avoids contaminating the ptrop array with 
	; ridiculous values 
	ensemble_member_copies = str_match_ind_ic(CopyMetaData,"ensemble member")
	ensemble_mean_copies = str_match_ind_ic(CopyMetaData,"ensemble mean")
	selected_copies = array_append_record (ensemble_member_copies, ensemble_mean_copies, 0)
	nC = dimsizes(selected_copies) 

	;=================================================;
	; Initialize output file 
	;=================================================;
	strs = str_split(filename_in,"/")
	filename_in_stripped = strs(dimsizes(strs)-1)
	filename_out = str_sub_str(filename_in_stripped, "Diag", "tropopause")
	system( "rm " + filename_out)
	fout = addfile(filename_out,"c")  

	;===================================================================
	; create coordinate variables for output file 
	;===================================================================
	dimNames = (/ "time","copy","lat","lon","lev","metadatalength" /)
	dimSizes = (/ time,copy,vdim(0),vdim(1),vdim(2),metadatalength /)
	dimUnlim = (/ False,False,False,False,False,False /)
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

	;===================================================================
	; predefine the the dimensionality of the variables to be written out
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
	filevardef(fout, "ptrop" ,typeof(fin->PS),getvardims(fin->PS))
	filevardef(fout, "P" ,typeof(fin->T),getvardims(fin->T))

	;===================================================================
	; Copy attributes associated with each variable to the file
	; All attributes associated with each variable will be copied.
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
	; Define attributes for tropopause pressure and the 3D pressure array. 
	;====================================================================
	fout->ptrop@units 	= fin->PS@units
	fout->ptrop@long_name	= "Lapse-rate tropopause pressure"
	fout->ptrop@_FillValue	= 9.96921e+36

	fout->P@units 		= fin->PS@units
	fout->P@long_name	= "Pressure"
	fout->P@_FillValue	= 9.96921e+36

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
	fout->lev    = (/fin->lev/) 
	fout->hyam    = (/fin->hyam/) 
	fout->hybm    = (/fin->hybm/) 
	fout->P0    = (/fin->P0/) 


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
	; loop over copy and compute 3d pressure array 
	;====================================================================
	do cc=0,copy-1
		; Compute pressures from hybrid level definition  
		ps   = fin->PS(:,cc,:,:)   ; surface pressure 
		p0   = fin->P0(0)	; reference pressure 
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
	; loop over selected copies and time, and compute tropopause pressure
	;====================================================================
	; loop over copy
	do icopy=0,nC-1
		cc = selected_copies(icopy)
		print("Computing tropopause height for copy "+CopyMetaData(cc,:))

		; loop over time 
		do tt=0,time-1

			; loop over all lats and lons 
			do ii=0,vdim(0)-1
				do jj=0,vdim(1)-1
					; vector of pressure at the given point 
					p = fout->P(tt,cc,ii,jj,:)
						   
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
