#########################################
# Configuration file for LRS2 reduction # 
#########################################

#This is the config file read by reduction_wrapper_lrs2.py and therefore needs to 

#***************************************#
# define paths to data and config files #
#***************************************#

#Choose to reduce either LRS2-Red, LRS2-Blue by specifying R or B for LRS2_spec
#Reduction is run on only one spectograph unit at a time because mastertrace and arcs will be different
#LRS2-Blue contains the UV (370-470nm) and orange (460-700nm) channels (B)
#LRS2-Red contains the red (650-842nm) and far-red (818-1050nm) channels (R)
LRS2_spec    	= 'B' 					#choose R or B only 

#path and name of the folder where reduction is run - NOTE: this folder is created by the script so the path has to exist but the folder does not
redux_dir       = "/home/04195/bindahl/shela_z4_731"	

#path to date folder containing raw data - This points to the folder where the data is stored on Maverick 
date_folder     = "/work/03946/hetdex/maverick/20160731"	

#Choose science objects from this night to reduce (These names can be found in the header keyword OBJECT)
#If left as empty list all science objects from that night will be reduced. 
#Take off the '_R' or '_B' at the end of the object name. The script will choose the correct one based on LRS2_spec you choose
sci_objects 	= ["SHELA_z4gal_0503678",]	

#Make sure to change this so it points to the lrs2_config folder in your LRS2_reduction directory
configdir       = "/home/04195/bindahl/LRS2_reduction/lrs2_config"	#path to your lrs2_config folder

#If you want to use a sky frame for sky subtraction set this to True
#This is only recommended for very bright object or extened sources that fill most of the field
#If False fibers in the science frame are used to build the sky model
#Script automatically finds sky frame from night observations and scales exposure
use_sky_frames 		= False 
#If you would like to force it to use a certain frame or the automatic routine does not work
#enter a list of exposure folders (full path) here, if empty will automatically search for frames
#ex. ["/work/03946/hetdex/maverick/20160731/lrs2/lrs20000005/exp01", "/work/03946/hetdex/maverick/20160731lrs2/lrs2000020/exp02"]
user_skyframes 		= []

#*******************************#
# choose reduction steps to run #
#*******************************#

basic           = True		#run basic reduction (overscan + bias subtract, ccd combine, build mastertrace + masterarc)
run_deformer    = True		#run deformer to map spectral traces and build wavelength solution for fiber extraction
subsky          = True		#run sky subtraction on sci images - Need to have run deformer, only runs on non-extracted spectra
fiberextract    = True 	 	#extract spectra and save into fits file - Need to have run deformer
makecube 		= True		#builds data cube out of fiber extracted image - Need to have run fiberextract
collapseCube 	= True		#collapse data cube to make an image of a wavelength range of the users choice

CLEAN_AFTER_DONE = True 	#If true it will delete intermediate reduction files for the calibration data

#*************************************#
# choose options for reductions steps #
#*************************************#

#------------#
# basic opts #
#------------#
subDarks 		= False			#[True/False] If True darks will be subtracted from science images (default: False)
dividePixFlt 	= False			#[True/False] If True images will be divided by pixel flats (default: False)
rmCosmics	 	= True  		#[True/False] If True the program L.A.Cosmic is used to eliminate cosmic rays (default: True)
#rmCosmics variables
sigclip 		= 5.0 			#[float] detection limit for cosmic rays (sigma). Increase if you detect cosmics where there are none (default: 5.0)
sigfrac 		= 0.3 			#[float] fractional detection limit for neighbouring pixels (default: 0.3)
objlim  		= 5.0 			#[float] contrast limit between CR and underlying object. Increase if normal stars are detected as cosmics (default: 5.0)

#----------------------#
# sky subtraction opts #
#----------------------#
window_size 	= 200 			#[integer] Size (in image lines) of moving window for sky median
sky_kappa 		=[3.5,3.5]		#[floatarray] Lower and upper kappa for final sky kappa-sigma clipping.
smoothing 		= 2.0 			#[float] Smoothing factor for approximating spline.
sn_thresh		= 15 			#[float] Minimum signal to noise to flag fiber as continuum and ignore it during sky generation.
sky_scaling 	= 1 			#[float] Scale the sky by this factor before subtracting. (** only for when sky frames used **)

#--------------------#
# fiber extract opts #
#--------------------#
wl_resample 	 = True 		#[True/False] If True it will resample in wavelength, Does not resample in wavelength if False (default: True)

#----------------#
# make cube opts #
#----------------#
sky_sampling 	= 0.3 			#[float] Regridded sample size on sky in arcsec. (default: 0.3)
max_distance	= 5.0 			#[float] Samples further away will haver weight=0 [arcsec]. (default: 5.0)
cube_sigma		= 0.75			#[float] Gaussian sigma of interpolation kernel [arcsec]. (default: 0.75)
diffAtmRef		= True  		#[True/False] Differential atmospheric refraction correction applied if True.

#--------------------#
# collapse cube opts #
#--------------------#
col_wave_range 	= [0,0]			#[floatarray] Choose the wavelength range (IN ANGSTROMS) that you would like to build and image of 
									#you can choose [0,0] if you would like the collapse the entire data cube 
