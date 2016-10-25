;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
; filter_DART_diagn.ncl
; Lisa Neef (lneef@geomar.de)
; October 2016  
;------------------------------------------------------------------------------
; This program runs through a list of DART output data (files with names like 
; modelname_Posterior_Diagn.nc) and runs them through the kf_filter code 
; for space-time filter. 
; The overall structure of this program is based on the example wkSpaceTime_5.ncl 
; by Carl Schreck   
;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/kf_filter.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

;<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



begin 
	dq         = str_get_dq()
	date_cmd   = "date +" + dq + "%a %b %d %H:%M:%S %Z %Y" + dq
	wcStrt     = systemfunc(date_cmd)

	;-------------------------------------------------
	;		filter setup
	;-------------------------------------------------

	; set the wave parameters
	tMin = 4.0
	tMax = 25.0
	kMin = 1
	kMax = 10
	hMin = 120
	hMax = 1200
	waveName = "dummy"  

	; set number of time steps per day  
	obsPerDay = 1  

	; set number of days to buffer on each side of the data  
	bufferDays = 60


	;-------------------------------------------------
	;		reading in data 
	;-------------------------------------------------

	print( "Reading in the available run data... " )

	; find and concatenate the requested file types  
	inDir = "/data/c1/lneef/DART-WACCM/nechpc-waccm-dart-gpsro-ncep-global-02/dart/hist/"
	files = systemfunc("ls "+"*Prior_TIL.20*0000.nc")
	print(files)
	f_all = addfiles(files,"r")
	ListSetType(f_all, "cat")


	;-------------------------------------------------
	;		read in variable values 
	;-------------------------------------------------


	; get the time units and set the first and last times to cover  
	; (the dates to which the buffer is then added when we run the filter)  

	; read the time units from the first file in the list 
	inFile=addfile(inDir+files(0),"r")
	timeUnits = inFile->time@units
	startDate = ut_inv_calendar( 2010, 01, 01, 00, 0, 0, timeUnits, 0 )
	endDate   = ut_inv_calendar( 2010, 01, 31, 00, 0, 0, timeUnits, 0 )
	
	; calculate number of time units to include for the buffer 
	if( isStrSubset( timeUnits, "days" ) ) then
		buffer = bufferDays * 1
	else
		if( isStrSubset( timeUnits, "hours" ) ) then
			buffer = bufferDays * 24
		else
			print( "DANGER! DANGER! Unexpected time units" )
			return
		end if
	end if

	filtStart = startDate - buffer
	filtStart@units = timeUnits
	filtEnd   = endDate + buffer
	filtEnd@units = timeUnits

	; NOTE: DART output data have shape time x copy x lat x lon x lev  
	ilat =67 
	ilev = 54
	icopy = 0

	print( "Reading in requested variable.... " )

	inData   = f_all[:]->brunt({filtStart:filtEnd},icopy,ilat,:,ilev)

	; TODO: decide if we want to average of lats of run filter over every lat   
	; average over selected latitude band 
	;origData = dim_avg_n_Wrap( inData, 1 )
	origData = inData


	;-------------------------------------------------
	;		filter
	;-------------------------------------------------


	print( "Filtering...")  

	filtData = kf_filter( origData(time|:,lon|:), obsPerDay \
		, tMin, tMax, kMin, kMax, hMin, hMax, waveName )


	printVarSummary( filtData )
	print( min(origData) + " " + max(origData) )
	print( min(filtData) + " " + max(filtData) )

	wallClockElapseTime(wcStrt, "Total Processing Time", 0)

end  
