==================================
# Instructions for LRS2 Reduction 
==================================

- This file contains instructions for running the CURE based reduction for LRS2 using the python script reduction_wrapper_lrs2.py
- For more information about CURE refer to the the CURE page on the HETDEX wiki: 
	https://luna.mpe.mpg.de/wikihetdex/index.php/Cure_-_the_data_analysis_system_for_HETDEX 

- All LRS2 data can be found on TACC's supercomputer Maverick inside the following directory:
	/work/03946/hetdex/maverick/

- It is recommended you run your reductions on TACC because CURE and python builds already exit. 
	- You will also not have to download and store data onto your computer.
	- However instructions are also provided for setting up reduction on your own computer if you wish to do so. 
	- NOTE: whether you run reduction on TACC or you own computer you will need a TACC account to access your data.

=================================================
# Understanding your data and the data stucture 
=================================================

------------------------------
1) Data structure on Maverick:
------------------------------

- On Maverick inside of the data folder (/work/03946/hetdex/maverick/) you will see all of the date folders
- Each date folder contains all of the data taken that night. 
- Inside each data folder are folders for the instruments used that night on the HET. Your data will be in 'lrs2'
- Inside the lrs2 folder are folders for all of the observations taken that night (cals and science targets) with lrs2
- Inside each observation folder there are folders for each exposure taken for that object (ex exp01)
- Inside each exposure folder there is a folder called 'lrs2' and that contains the 4 or 8 fits files for that exposure 

- An exposure with one LRS2 unit gives you 4 fits files. There are two detectors in each unit and each has 2 amplifiers.
- If you observed with both units (red and blue) you will have 8 fit files per exposure. 

-------------------------------
1) Information from file names:
-------------------------------

- Here is an example filename: 20161202T074404.6_056LL_cmp.fits
- The part before the first underscore (20161202T074404.6) is the date and time the image was taken (the 'T' separates the date and time)

		(format: yyyymmddThhmmss.s)

- The 3 digit number after the first underscore (056) tells you the IFU slot ID for the unit. This tells you which LRS2 unit this is from

		LRS2-Blue: 056
		LRS2-Red : 066

- The next lettter will either be 'L' or 'R' for left or right. This tells you if this is the left or right side of the unit or which channel.

		LRS2-B left side : UV channel      (056L)
		LRS2-B right side: orange channel  (056R)
		LRS2-R left side : red channel     (066L)
		LRS2-R right side: far-red channel (066R)

- The letter after that will either be 'L' or 'U' for lower or upper. This tells with amplifier of that detector it is. 
   * NOTE: basic reduction will orient U and L properly and combine them.
   * So after reduction instead of LU, LL, RU, RL you will just have L and R images. 
- The part after the second underscore tells you the image type. There are 5 image types:

		zro - bias
		flt - flats (taken with either the LDLS or Qth lamps)
		cmp - arc lamps or comps (taken with either Cd, Hg, FeAr, or Kr)
		drk - darks
		sci - science frames (twilight flats are also sci)

================================================
# Running LRS2 reduction on TACC - RECOMMENDED 
================================================

-------------------------------------------------------
1) Getting a TACC account and joining the HETDEX group:
-------------------------------------------------------

- If you don't have an account on TACC you need to set one up.

		go to:  https://portal.tacc.utexas.edu/
		Then, click: Create a TACC account

- Send Karl (gebhardt@astro.as.utexas.edu) your username so he can add you to the HETDEX group.

- All data for LRS2 and VIRUS is stored on the supercomputer Maverick. Everything will be found and run on here. 

		>>> ssh username@maverick.tacc.utexas.edu

- After ssh'ing into Maverick, you will be at your home page.  

---------------------------------------------
2) Setting up TACC account to run reductions:
---------------------------------------------

- CURE and all python packages needed for reduction are already on TACC.
- You must first set up these paths in you "~/.bashrc" file in your account on Maverick.

- Open your "~/.bashrc" file and add the following lines at the end: 
	(NOTE: The last 3 are optional but are useful for running VIRUS reduction scripts or viewing fits files)

		umask 022
		module load  intel/15.0.3 cxx11
		export WORKON_HOME=$HOME/.virtualenvs
		export PROJECT_HOME=$HOME/Devel
		export PATH="/home/00115/gebhardt/anaconda2/bin:/work/03946/hetdex/maverick/bin:$PATH"
		export CURELRS2="/home/04195/bindahl/curelrs2/cure/bin"
		export CUREBIN="/home/03946/hetdex/curehome/cure/bin"
		alias ds9=/work/02426/ngaffney/saods9/bin/ds9
		alias QFitsView='/home/01821/papovich/qfitsview/QFitsView'

- Save the file and exit the text editor 

- Source your "~/.bashrc" file to make the changes 

		>>> source ~/.bashrc

- Change permissions of your home directory (cd) for accesibility

		>>> cd
		>>> chmod a+rx ../username

-----------------------------------------------
3) Obtaining reduction script and config files: 
-----------------------------------------------

- Inside your work directory (cdw) make a copy of the LRS2_reduction folder with the following commands: 

		>>> cdw 
		>>> cp -r /home/04195/bindahl/LRS2_reduction ./ 

- You should now have a folder in your work directory called LRS2_reduction. This folder contains the following files and folder:

	1. __reduction_wrapper_lrs2.py__ - This is the folder that runs the reductions. You should never have to edit this file. 
	2. __lrs2_config.py__ 			 - This is the config file where the user defines the data and opts for their reduction
	3. __cosmics.py__ 				 - This is the script that runs L.A.comsic in the reduction (http://obswww.unige.ch/~tewes/cosmics_dot_py/)
	4. __lrs2_config__ 				 - This is a folder that contains all of the configurations files needed for LRS2 reduction 
	  * __lines_files__	     - These files defined the pixel and wavelength to find the arc lines for building the wavelength solution
	  * __mapping_files__	 - These files contain the mapping of the fibers onto the field for building data cubes
	  * __pixel_flats__	     - These files are the pixels flats for each CCD that can be optionally divided during reduction
	  * __longExpCals__	     - These are 1800sec FeAr exposures used for pinning down the wavelength solution for the far-red channel
	  * __short_OrgFlts__	 - These replace longer exposure flats in early LRS2 data that were saturating the orange channel

--------------------------
4) Running LRS2 reduction:
--------------------------

- open your lrs2_config.py file in a text editor 
- edit fields and paths according to the comments to define your data to be reduced:
	* choose the steps of reduction to run, and opts for those steps
- save the changes to lrs2_config.py 
- run reduction_wrapper_lrs2.py (if run outside of LRS2_reduction you must give the path to reduction_wrapper_lrs2.py)

		>>> python reduction_wrapper_lrs2.py 

- When the reduction script runs it will automatically save a copy of the config file used in the redux directory
- You can also choose to run the reduction script with a different config file 
	* For example you could make your own config file copies with settings you want to save 
	* or you can use the saved config files in one of your redux directories. 
- To run the reduction script with a different config file:
	* First you must copy the config file to your LRS2_reduction directory if it is not already there
	* Run the script with the following command

		>>> python reduction_wrapper_lrs2.py -config "path/name_of_config.py"

----------------------------
5) Upadating LRS2 reduction:
----------------------------

- I recommend periodically making sure your version of the LRS2 reduction software is up to date
- Changes to the script can be tracked on my LRS2_reduction github: https://github.com/BrianaLane/LRS2_reduction
- If you copy the entire LRS2_reduction directory to your work folder you can pull changes made to your version
- To update your LRS2_reduction version in your work directory:

		>>> cd LRS2_reduction
		>>> git pull

- Edits you may have made to lrs2_config.py may be incompatable with the new version you are trying to pull
- You can disregard the changes to pull in the update with the following command inside your LRS2_reduction directory: 

		>>> git stash
		>>> git pull

===============================================
# Running LRS2 reduction on your own computer 
===============================================

-------------------------------------------
1) Obtaining and configuring CURE for LRS2:
-------------------------------------------

- If you don't already have CURE installed you have to do that first
	- Instructions for installing CURE can be found in the CURE cookbook on HETDEX wiki (link at top of page)

- Once CURE is installed you have to edit specconf.h to set up CURE for LRS2 data 
	- Open specconf.h in a text editor (file found inside cure/libcure/)
	- Edit line 34 to define LRS2 instead of VIRUS_HET (#define LRS2). Save file and exit
	- You must recompile CURE to make the change 
		- cd into your cure directory and run the following

				>>> make clean
				>>> make install 

- Add this to your "~/.bashrc" file to set the path to your CURE bin:

		export CURELRS2="/path_to_folder_containing_CURE/cure/bin"

------------------------------------
2) Python build and packages needed:
------------------------------------

- You need to have Python 2.7 or later installed 
- You need the following Python packages: numpy, scipy, pyfits, glob

-------------------------------------
3) Obtaining LRS2 data from Maverick: 
-------------------------------------

- You must scp your data off of Maverick from /work/03946/hetdex/maverick/ onto your computer
- It is important that you maintain the same folder structure which is: 

		date_folder/lrs2/lrs000####/exp##/lrs2/*.fits 

- When running the reduction you will be defining the path to the date_folder in which to find your data 

---------------------------------------------------------------
4) Obtaining reduction script/files and running LRS2 reduction: 
---------------------------------------------------------------

- Refer to sections 3 and 4 under 'Running LRS2 reduction on TACC' above
- The one difference is you will scp the LRS2_reduction directory onto your computer instead of copying it into your TACC directory

=========================================
# Understanding reduction data products 
=========================================

-------------------------------------------
1) Files + folders in your redux directory:
-------------------------------------------

- You will see a ton of files that are all of the master files for the arcs, flats, biases, and darks. 
- There are also mastertrace files with .pmod, .fmod, .dist, ect. These are the data products from deformer needed to run later reduction steps

- If you choose CLEAN_AFTER_DONE to be True then the folders cmp, flt, zro, and drk will be empty 
	* otherwise they are filled with the intermediate files from reducing calibration data to build master files 

- The folder sci contains all of your reduced science images, sky subtracted files, fiber extracted files, data cubes, and collapsed cubes

----------------------
2) Your science files: 
----------------------

- Refer to 'Understanding your data and the data stucture' to understand the filename convention 

- You will see a bunch of versions of your science files with different prefixs appended to there name
- After each step of the reduction a letter is added to the filename as follows:

		pses 	- These are files that have been run through basic reduction 
		S       - These are files that have also been sky subtracted 
		Fe 	    - These are files that have been fiber extracted WITHOUT wavelength resampling 
		FeR 	- These are files that have been fiber extracted WITH wavelength resampling
		Cu 	    - These files are data cubes 
		Col 	- These files are collapsed data cubes 
		e. 	    - These are the error files for all of these frames

- As an example for the file: CuFeRSpses20160731T094822.4_056_sci_R.fits
	- The pses shows that it has been through basic reduction (pses), sky subtracted (S), 
		fiber extracted with wavelength resampling (FeR), and then build into a data cube (Cu)
	- The error file for this frame would be: e.CuFeRSpses20160731T094822.4_056_sci_R.fits

-----------------------------------
2) Finding the Wavelength Solution: 
-----------------------------------

-If you have run fiber extraction with wavelength resampling then you 
