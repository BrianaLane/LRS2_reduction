#########################################
# Configuration file for LRS2 reduction # 
#########################################

#reduction_wrapper_lrs2.py should be run inside of a folder containing the date folder containing raw data

#Choose to reduce either LRS2-Red, LRS2-Blue by specifying R or B for LRS2_spec
#Reduction is run on only one spectograph unit at a time because mastertrace and arcs will be different
#LRS2-Blue contains the UV and orange channels (B)
#LRS2-Red contains the red and far-red channels (R)
LRS2_spec    	= 'R' 		#choose R or B only 

redux_dir       = "shela_CL2278A_test"	 	#name of the folder reduction is run 
date_folder     = "20160909"				#date folder containing raw data 

zro_folder      = ["lrs20000007",]		#folder containing raw bias images 
Hg_folder       = ["lrs20000002",]		#folder containing raw Hg arc lamp images 
Cd_folder       = ["lrs20000003",]		#folder containing raw Cd arc lamp images (Only needed for LRS2-B)
FeAr_folder     = ["lrs20007999",]		#folder containing raw FeAr arc lamp images (Only needed for LRS2-R)
flt_folder      = ["lrs20009920",]		#folder containing raw flat images (LRS2-R use Qth) (LRS2-B use LDLS)
sci_folder      = ["lrs20000014",]	#List of folders containing the science observations for that night (Must all be taken in the same night)

configdir       = "/Users/Briana/Documents/Grad_School/LRS2/lrs2_config"	#path to lrs2_config folder

basic           = False		#run basic reduction (normalize, build mastertrace and masterarc)
run_deformer    = False		#run deformer to map spectral traces and build wavelength solution for fiber extraction
subsky          = False		#run sky subtraction on sci images - Need to have run deformer, only runs on non-extracted spectra
fiberextract    = True	 	#extract spectra and save into fits file - Need to have run deformer
makecube 		= False		#builds data cube out of fiber extracted image - Need to have run fiberextract

CLEAN_AFTER_DONE = True 	#If true it will delete intermediate reduction files for the calibration data

########################
# sky subtraction opts #
########################
window_size 	= 250 			#[integer] Size (in image lines) of moving window for sky median
sky_kappa 		=[3.5,3.5]		#[floatarray] Lower and upper kappa for final sky kappa-sigma clipping.
smoothing 		= 5.0 			#[float] Smoothing factor for approximating spline.
sn_thresh		= 5 			#[float] Minimum signal to noise to flag fiber as continuum and ignore it during sky generation.

######################
# fiber extract opts #
######################
wl_resample 	 = True 		#[True/False] If True it will resample in wavelength, Does not resample in wavelength if False 

##################
# make cube opts #
##################
sky_sampling 	= 0.3 			#[float] Regridded sample size on sky in arcsec.
max_distance	= 5.0 			#[float] Samples further away will haver weight=0 [arcsec]. 
cube_sigma		= 0.75			#[float] Gaussian sigma of interpolation kernel [arcsec].
diffAtmRef		= True 			#[True/False] Differential atmospheric refraction correction applied if True.

