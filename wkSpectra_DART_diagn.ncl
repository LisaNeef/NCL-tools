;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
; wkSpectra_DART_diagn.ncl
; Lisa Neef (lneef@geomar.de)
; October 2016  
;------------------------------------------------------------------------------
; This program runs through a list of DART output data (files with names like 
; modelname_Posterior_Diagn.nc) and runs them through the wkSpectra tool, 
; which creats Wheeler and Kiladis-type wavenumber frequency spectra. 
; The overall structure of this program is based on the example wkSpaceTime_5.ncl 
; by Carl Schreck   
;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"

;<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



begin 
	dq         = str_get_dq()
	date_cmd   = "date +" + dq + "%a %b %d %H:%M:%S %Z %Y" + dq
	wcStrt     = systemfunc(date_cmd)

	;-------------------------------------------------
	;		spectrum tool setup
	;-------------------------------------------------

	diro 	= systemfunc("pwd ")+"/../Plots/"	; directory where the output plots go  
	varName = "brunt"			; name of the variable to be processed 
	diagnostic = "Prior"			; which DART diagnostic? 

	; set up latitude boundaries and select a level
	latN	= 60				; northern latitude boundary 
	latS	= 45				; southern latitude boundary 
	lev_in	= 300				; vertical level  

	; time settings  
	spd = 1					; samples per day in the input data  

	; define windows of which to apply the FFT  
	nDayWin = 90		; length of the temporal window in days  
	nDaySkip = -30				; number of days to skip between temporal windows 

	;-------------------------------------------------
	;		reading in data 
	;-------------------------------------------------

	print( "Reading in the available run data... " )

	; find and concatenate the requested file types  
	inDir = "/data/c1/lneef/DART-WACCM/nechpc-waccm-dart-gpsro-ncep-global-02/dart/hist/"
	files = systemfunc("ls *"+diagnostic+"*.20*0000.nc")
	f_all = addfiles(files,"r")
	ListSetType(f_all, "cat")


	;-------------------------------------------------
	;	retrieve grid and time info 
	; 	from a single file 
	;-------------------------------------------------

	; read the time units from the first file in the list 
	inFile=addfile(inDir+files(0),"r")
	timeUnits = inFile->time@units
	filtStart = ut_inv_calendar( 2009, 10, 01, 00, 0, 0, timeUnits, 0 )
	filtEnd = ut_inv_calendar( 2010, 02, 28, 00, 0, 0, timeUnits, 0 )
	filtStart@units = timeUnits
	filtEnd@units = timeUnits


	; identify the copy indices where computing the spectra makes sense,
	; i.e. ensemble members and the mean, but not the ensemble spread 
	; or inflation factors 

	CopyMetaData = inFile->CopyMetaData
	ensemble_member_copies = str_match_ind_ic(CopyMetaData,"ensemble member")
	ensemble_mean_copies = str_match_ind_ic(CopyMetaData,"ensemble mean")
	selected_copies = array_append_record (ensemble_member_copies, ensemble_mean_copies, 0)
	nC = dimsizes(selected_copies) 
	

	;-------------------------------------------------
	;		loop over copies and plot spectra 
	;-------------------------------------------------
	do icopy=0,nC-1

		cc = selected_copies(icopy)
		print( "Creating and plotting spectra for "+CopyMetaData(cc,:))  

		copyname = str_squeeze(CopyMetaData(cc,:))
		if (latN .lt. 0) then
			latNname = tostring(latN)+"S"
		else 
			latNname = tostring(latN)+"N"
		end if 
		if (latS .lt. 0) then
			latSname = tostring(latS)+"S"
		else 
			latSname = tostring(latS)+"N"
		end if 

		caseName = copyname+" "+diagnostic+" "+tostring(lev_in)+"hPa "+latSname+"-"+latNname
		

		;-------------------------------------------------
		;	reading in multi-file data 
		;-------------------------------------------------

		print( "Reading in requested variable.... " )

		inData   = f_all[:]->brunt({filtStart:filtEnd},cc,{latS:latN},:,{lev_in})
		time = f_all[:]->time({filtStart:filtEnd})


		;-------------------------------------------------
		;	replace FillValues with time mean 
		;-------------------------------------------------

		; FFT doesn't like FillValues, so we need to replace them. 
		; an easy option is to fill in the time-mean 

		; print an alert about missing data 
		if (any(ismissing(inData))) then
			N = num(ismissing(inData))
			print("*******************")
			print("We have "+tostring(N)+" missing values -- fill these in with the time mean")

			meanData = dim_avg_n_Wrap( inData, 0 )
			meanDataEx = conform_dims(dimsizes(inData),meanData,(/1,2/))

			Data1d = ndtooned(inData)
			meanData1d = ndtooned(meanDataEx)
			Data1d(ind(ismissing(Data1d))) = meanData1d(ind(ismissing(Data1d)))

			Data = onedtond(Data1d,dimsizes(inData))
			copy_VarCoords(inData,Data)
			copy_VarAtts(inData,Data)

			; print another alert if some data are still missing  
			if (any(ismissing(Data))) then
				N = num(ismissing(Data))
				print("*******************")
				print("We still have "+tostring(N)+" missing values!")
			end if 
		else
			Data = inData
		end if 

		;-------------------------------------------------
		;	plot the wavenumber-frequency spectrum  
		;-------------------------------------------------


		; some figure settings
		opt=True
		opt@pltTitle=caseName  
		opt@pltType="eps"
		;opt@Fig_1=fspan(-3.3,0.9,15)
		;opt@Fig_2=fspan(-1.2,1.2,16)
		;opt@Fig_3=fspan(0.2,2.8,15)

		wkSpaceTime (Data, diro, caseName, varName, latN, spd, nDayWin, nDaySkip, opt)  


	end do ; loop over copies  

end  
