#########################################
# Configuration file for LRS2 reduction # 
#########################################

#reduction_wrapper_lrs2.py should be run inside of a folder containing the date folder containing raw data

#Choose to reduce either LRS2-Red, LRS2-Blue by specifying R or B for LRS2_spec
#Reduction is run on only one spectograph unit at a time because mastertrace and arcs will be different
#LRS2-Blue contains the UV (370-470nm) and orange (460-700nm) channels (B)
#LRS2-Red contains the red (650-842nm) and far-red (818-1050nm) channels (R)
LRS2_spec    	= 'R' 		#choose R or B only 

redux_dir       = "SN_IPTFxx_W"	 	#name of the folder reduction is run - folder created by script
date_folder     = "20160506"				#date folder containing raw data 

zro_folder      = ["lrs20000012",]		#[need]   folders containing raw bias images 
drk_folder		= []					#[opt]    folders containing darks for science frames
Hg_folder       = ["lrs20009990",]		#[need]   folders containing raw Hg arc lamp images 
Cd_folder       = []					#[only-B] folders containing raw Cd arc lamp images 
FeAr_folder     = ["lrs20007999",]		#[only-R] folders containing raw FeAr arc lamp images
flt_folder      = ["lrs20000009",]		#[need]   folders containing raw flat images (LRS2-R use Qth) (LRS2-B use LDLS)
sci_folder      = ["lrs20000014",]		#[need]   folders containing the science observations for that night 

configdir       = "/Users/Briana/Documents/Grad_School/LRS2/LRS2_reduction/lrs2_config"	#path to lrs2_config folder

basic           = True		#run basic reduction (normalize, build mastertrace and masterarc)
run_deformer    = True		#run deformer to map spectral traces and build wavelength solution for fiber extraction
subsky          = True		#run sky subtraction on sci images - Need to have run deformer, only runs on non-extracted spectra
fiberextract    = True	 	#extract spectra and save into fits file - Need to have run deformer
makecube 		= False		#builds data cube out of fiber extracted image - Need to have run fiberextract
collapseCube 	= False		#collapse data cube to make an image of a wavelength range of the users choice

CLEAN_AFTER_DONE = True 	#If true it will delete intermediate reduction files for the calibration data

##############
# basic opts #
##############
dividePixFlt 	= False			#[True/False] If True images will be divided by pixel flats (default: False)
rmCosmics	 	= True  		#[True/False] If True the program LAcosmics is used to eliminate cosmic rays (default: True)

########################
# sky subtraction opts #
########################
window_size 	= 200 			#[integer] Size (in image lines) of moving window for sky median
sky_kappa 		=[3.5,3.5]		#[floatarray] Lower and upper kappa for final sky kappa-sigma clipping.
smoothing 		= 2.0 			#[float] Smoothing factor for approximating spline.
sn_thresh		= 15 			#[float] Minimum signal to noise to flag fiber as continuum and ignore it during sky generation.

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
diffAtmRef		= True  		#[True/False] Differential atmospheric refraction correction applied if True.

######################
# collapse cube opts #
######################
col_wave_range 	= [6500,6700]	#[floatarray] Choose the wavelength range (IN ANGSTROMS) that you would like to build and image of 
									#you can choose [0,0] if you would like the collapse the entire data cube 
