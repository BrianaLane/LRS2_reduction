# -*- coding: utf-8 -*-
"""
Created on Friday April  8 02:10:02 2016

Reduces LRS2 data for either the blue or red channels
This script reads user parameters from lrs2_config.py 

@author: brianaindahl, gregz
"""
from __future__ import print_function

import numpy as np
from scipy import interpolate
import pyfits
import glob
import six
from datetime import datetime
import shutil
import sys
import os
import copy
import os.path as op
from os import environ
import re
import string
import cosmics 
import argparse as ap
import importlib
#from lrs2_config import * 

#blu 370-470nm
#org 460-700nm
#red 650-842nm
#frd 818-1050nm

#############################
# Define config file to use #
#############################

def parse_args(argv=None):
    """Parse the command line arguments
    Parameters
    ----------
    argv : list of string
        arguments to parse; if ``None``, ``sys.argv`` is used
    Returns
    -------
    Namespace
        parsed arguments
    """
    description = "config_file"
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.ArgumentDefaultsHelpFormatter)
                        
    parser.add_argument("-config", nargs='?', type=str, help='''config file. "path_to_config/lrs2_config.py"''', default="lrs2_config.py")

    args = parser.parse_args(args=argv) 
           
    return args

args = parse_args()

config_file_name = args.config
config_arg       = config_file_name.split(".")[0]

config = importlib.import_module(config_arg, package=None)

#################################
# Defining which unit to reduce #
#################################

#specifying LRS2 unit to reduce 
if config.LRS2_spec == 'B':
    redux_dir   = config.redux_dir+'_B'
    print ('#########################################')
    print ('## RUNNING DATA REDUCTION ON LRS2-BLUE ##')
    print ('#########################################')
elif config.LRS2_spec == 'R':
    redux_dir   = config.redux_dir+'_R'
    print ('########################################')
    print ('## RUNNING DATA REDUCTION ON LRS2-RED ##')
    print ('########################################')
else:
    sys.exit('You need to choose either R or B for config.LRS2_spec')

##################################################
# Setting CUREBIN and check LRS2 defined in CURE #
##################################################

#setting CUREBIN
CURELRS2 = None
if not CURELRS2:
    CURELRS2 = environ.get('CURELRS2')
if not CURELRS2:
    sys.exit("Please set CURELRS2 as environment variable")

#checking that LRS2 is defined in specconf.h 
cureversion = os.popen(op.join(CURELRS2, 'cureversion')).readlines()
spec_define = cureversion[4].split(' ')[1]
instrument = spec_define.rstrip('\n')

if instrument == 'LRS2':
    print ('CURE is set for LRS2 reduction')
else: 
    print ('You need to update specconf.h in CURE to define LRS2')
    sys.exit('Right now specconf.h defines '+instrument)

###################################
# Defining which functions to run #
###################################

#if basic reduction is run need to specify specific routines to run 
# divide pixel flat and masterdark are not being used now
if config.basic:
    rmcosmics       = config.rmCosmics 
    fix_chan        = True
    dividepf        = config.dividePixFlt
    normalize       = False
    masterdark      = config.subDarks
    subtractdark    = config.subDarks
    masterarc       = True  
    mastertrace     = True 
    sort_sci        = True
else:
    rmcosmics       = False
    fix_chan        = False
    dividepf        = False
    normalize       = False  
    masterdark      = False
    masterarc       = False  
    mastertrace     = False
    sort_sci        = False

# This makes sure that the redux folder is only overwritten if the user chooses to run basic reduction
# If you user only wants to run deformer, skysubtract, fiberextract, or mkcube it used the data in redux 
if config.basic:
    all_copy = True
    RESTART_FROM_SCRATCH = True
else:
    all_copy = False
    RESTART_FROM_SCRATCH = False

#sets the initial action base to start with. CURE adds bases to data when a change is made 
initial_base     = ''

########################
# Specifying CURE opts #
########################

#specify opts for CURE routines used in this script
meanfitsopts    = "--new_error -k 3"
headfitsopts    = "-m -k EXPTIME -v 1 -V"
darkopts        = "--maxmem 1024 -s -t -m -k 2.8" 
arcopts         = "--maxmem 1024 -s -t -m -k 2.8"
traceopts       = "--maxmem 1024 -s -t -m -k 2.8"
deformeropts    = "-p 7 -n 4 -C 10 --debug --dump_psf_data"
subskyopts      = "-J --output-both -w "+str(config.window_size)+" -k "+str(config.sky_kappa[0])+','+str(config.sky_kappa[1])+" -m "+str(config.smoothing)+" -T "+str(config.sn_thresh)
fibextractopts  = "-P"
cubeopts        = "-a "+str(config.sky_sampling)+" -k "+str(config.max_distance)+" -s "+str(config.cube_sigma)

#########################
# Defining data folders #
######################### 

zro_dir  = "zro"
flt_dir  = "flt"
sci_dir  = "sci"
cmp_dir  = "cmp"
drk_dir  = "drk"

# Dictionary for looping through directory names
DIR_DICT     = {  0:zro_dir, 1:drk_dir, 2:cmp_dir, 3:flt_dir, 4:sci_dir } 

##########################
# Building Spec Libaries #
##########################

#Detector Amps and spectrograph side lists 
SPEC = ["LL","LU","RL","RU"]
SPECBIG = ["L","R"]

#############################
# Define config directories #
#############################

#specifies directories for lines and mapping/cen files for LRS2 in the LRS2 config directory 
linesdir    = config.configdir + '/lines_files'
mapdir      = config.configdir + '/mapping_files/'
pixflatdir  = config.configdir + '/pixel_flats/'

##################
# Define Classes #
##################

#class to build VIRUS frames for all data 
class VirusFrame:
    def __init__ ( self, filename = None):
        '''
        Initializing a VirusFrame for a given filename.
        This includes reading the header and setting reduction keywords.
        From this, the file with have attributes that can be tracked.
        '''
        
        if filename is None:
            print ( "No filename was given for VirusFrame to be initialized." )
            return None
        else:
            ######## OPEN AND KEEP FILES IN SOME SMART WAY ########
            self.filename               = filename
            self.origloc                = op.dirname  ( self.filename )
            self.basename, temp1, temp2 = op.basename ( self.filename ).split('_')
            self.type                   = temp2.split ( '.fits', 1 )[0]
            self.ifuslot                = temp1[0:3]
            self.time                   = self.basename.split('T')[1]
            self.hr                     = int(self.basename.split('T')[1][0:2])
            self.min                    = int(self.basename.split('T')[1][2:4])
            self.sec                    = float(self.basename.split('T')[1][4:8])
            self.date                   = self.basename.split('T')[0]
            self.year                   = int(self.basename.split('T')[0][0:4])
            self.month                  = int(self.basename.split('T')[0][4:6])
            self.day                    = int(self.basename.split('T')[0][6:8])
            self.clean                  = config.CLEAN_AFTER_DONE
            self.trimsec                = "2:2065,1:1032" 
            self.biassec                = "2066:2128,1:1032" 

            ###### READ IN THE HEADER AND NECESSARY KEYWORDS ######
            self.actionbase = {}
            for amp in SPEC:  
                self.actionbase[amp] = ''   

            self.actionbase["L"] = initial_base  
            self.actionbase["R"] = initial_base  

            rootname          = op.join (self.origloc, self.basename + '_' + self.ifuslot + 'LL_' + self.type + '.fits' )
            hdulist           = pyfits.open ( rootname )       
            self.specid      = str(hdulist[0].header['SPECID']) 
            self.orggain     = hdulist[0].header['GAIN']
            self.orgrdnoise  = hdulist[0].header['RDNOISE']
            self.exptime     = hdulist[0].header['EXPTIME']
            self.ccdpos        = hdulist[0].header['CCDPOS']

            #for cal data: object will have R or B appended for the liquid light guide used 
            #for sci data: object will have R or B appended for the unit that was pointed at the sci target
            #this is not true for older data and is sometimes inconsistent so if no tag cal_side = None
            if len(hdulist[0].header['OBJECT'].split('_')) == 1:
                self.object   = hdulist[0].header['OBJECT']
                self.cal_side = None
            elif len(hdulist[0].header['OBJECT'].split('_')) > 1:
                #sci objects and flats can have multiple '_' so have to check if the last one is a tag R or B
                end_term = hdulist[0].header['OBJECT'].split('_')[-1]
                if (end_term == 'R') or (end_term == 'B'):
                    self.object   = hdulist[0].header['OBJECT'][0:-2]
                    self.cal_side = end_term
                else:
                    self.object   = hdulist[0].header['OBJECT']
                    self.cal_side = None

                # if len(hdulist[0].header['OBJECT'].split('_')) == 4 and self.type == 'flt':
                #     self.length   = hdulist[0].header['OBJECT'].split('_')[2]
                # else:
                #     self.length   = None
                

    def addbase(self, action, amp = None, side = None):
        if self.clean:
            if amp is not None:
                filename   = op.join ( self.origloc,        self.actionbase[amp] + self.basename + '_' + self.ifuslot + amp + '_' + self.type + '.fits')
                filename_e = op.join ( self.origloc, 'e.' + self.actionbase[amp] + self.basename + '_' + self.ifuslot + amp + '_' + self.type + '.fits')            
                self.actionbase[amp] = action + self.actionbase[amp]
                if op.exists  ( filename ):
                    os.remove ( filename )
                if op.exists  ( filename_e ):
                    os.remove ( filename_e )
                    
            if side is not None:
                filename   = op.join ( self.origloc,        self.actionbase[side] + self.basename + '_' + self.ifuslot + '_' + self.type + '_' + side + '.fits')
                filename_e = op.join ( self.origloc, 'e.' + self.actionbase[side] + self.basename + '_' + self.ifuslot + '_' + self.type + '_' + side + '.fits')            
                self.actionbase[side] = action + self.actionbase[side]
                if op.exists  ( filename ):
                    os.remove ( filename )
                if op.exists  ( filename_e ):
                    os.remove ( filename_e )                           
        else:
            if amp is not None:
                self.actionbase[amp] = action + self.actionbase[amp]
                
            if side is not None:
                self.actionbase[side] = action + self.actionbase[side]

#Class for building dither file for CURE's mkcube
class ditherinfo(object):
    @classmethod
    def writeHeader(cls, f):
        """Write the header to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = []
        s.append("# basename is the base file name of the data files.")
        s.append("# modelbase is the base file name of the dist, fmod, pmode files corresponding to the data files")
        s.append("# $Id:$")
        s.append("#")
        s.append("# basename modelbase ditherx dithery seeing norm airmass")
        s.append("#")
        f.write('\n'.join(s) + "\n")

    @classmethod
    def writeDither(cls, f, filename, basename, ditherx, dithery, seeing, norm, airmass):
        """Write something to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = ("%s %s %4.2f %4.2f %4.2f %4.2f %4.2f" %
             (filename, basename, ditherx, dithery, seeing, norm, airmass))
        f.write(s)
        f.write("\n")
        f.flush()

#####################################################
# Define For Running Data Reduction Steps Functions #
#####################################################

def run_cure_command(command, suppress_output=0, debug=1):
    '''
       Run a cure command
       
       Parameter suppress_output is used to quiet down. (default)
       Return the value of os.system(command)
    '''
    # store info in log?

    if debug:
        print('Running: \'%s\'' % command)
    if not suppress_output:
        return os.system(op.join(CURELRS2, command) +' 1>>output.log  2>> error.log')
    else:
        return os.system(op.join(CURELRS2, command))
        
        
def mkerrorframe(frames, amp):
    
    filenames = [op.join ( f.origloc, f.actionbase[amp] + f.basename + '_' + f.ifuslot + amp + '_' + f.type + '.fits') for f in frames]

    command = 'mkerrorframe' 
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
    
    run_cure_command( command, 0 )
        
    return command
    
        
def subtractoverscan(biassec, frames, amp):
    
    filenames = [op.join ( f.origloc, f.actionbase[amp] + f.basename + '_' + f.ifuslot + amp + '_' + f.type + '.fits') for f in frames]

    command = 'subtractfits -s -a -k 2.8 -t -o %s -z' % (biassec)

    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    [f.addbase ('s',amp=amp) for f in frames] 

    return command
    
    
def subtractbias(frames, masterbiasname, amp):
        
    filenames = [op.join ( f.origloc, f.actionbase[amp] + f.basename + '_' + f.ifuslot + amp + '_' + f.type + '.fits') for f in frames]

    command = 'subtractfits -f %s' % ( masterbiasname ) 
        
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    [f.addbase ('s',amp=amp) for f in frames] 

    return command

def subtractdark(frames, masterdarkname, drk_scale, amp):
        
    filenames = [op.join ( f.origloc, f.actionbase[amp] + f.basename + '_' + f.ifuslot + amp + '_' + f.type + '.fits') for f in frames]
        
    #scale the masterdark to the exposure time of the science frames
    command_1 = 'multiplyfits -p d%s_ -c %s' % ( drk_scale,drk_scale ) 
    command_1 = command_1 + ' ' + masterdarkname    
    run_cure_command( command_1, 0 )

    #subtract the scaled masterdark from the science frames 
    command_2 = 'subtractfits -p '' -f %s' % ( masterdarkname ) 
    command_2 = command_2 + ' ' + filenames
    run_cure_command( command_2, 0 )

    return command_2
    

def meanbiasfits(amp, specid, dest_dir, basename , opts, frames):
    
    filenames = [op.join ( f.origloc, f.actionbase[amp] + f.basename + '_' + f.ifuslot + amp + '_' + f.type + '.fits') for f in frames]

    mastername = op.join ( dest_dir , basename + '_' + specid + '_' + amp + '.fits')
    
    command = 'meanfits %s -o %s' % (opts,mastername)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    [f.addbase ('master',amp=amp) for f in frames]

    return command
    
    
def meandarkfits(side, specid, dest_dir, basename , opts, frames):
    
    filenames = [op.join ( f.origloc, f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in frames]

    mastername = op.join ( dest_dir , basename + '_' + specid + '_' + side + '.fits')
    
    command = 'meanfits %s -o %s' % (opts,mastername)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    [f.addbase ('master',side=side) for f in frames]

    return command
    

def extractfits(trimsec, frames, amp):

    filenames = [op.join ( f.origloc, f.actionbase[amp] + f.basename + '_' + f.ifuslot + amp + '_' + f.type + '.fits') for f in frames]
    
    command = 'extractfits -r %s' % (trimsec)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command ( command, 0 )
    
    [f.addbase ('e',amp=amp) for f in frames] 
    
    return command

def rmcosmicfits(frames, amp):

    filenames = [op.join ( f.origloc, f.actionbase[amp] + f.basename + '_' + f.ifuslot + amp + '_' + f.type + '.fits') for f in frames]
    
    for i in xrange(len(filenames)):

        print (filenames[i])

        array, header = cosmics.fromfits(filenames[i])
        # array is a 2D numpy array

        im_gain = header['GAIN']
        im_RN = header['RDNOISE']

        # Build the object :
        c = cosmics.cosmicsimage(array, gain=im_gain, readnoise=im_RN, satlevel = 65535.0, sigclip = config.sigclip, sigfrac = config.sigfrac, objlim = config.objlim)
        # There are other options, check the manual...

        # Run the full artillery :
        c.run(maxiter = 4)

        # Write the cleaned image into a new FITS file, conserving the original header :
        cosmics.tofits(filenames[i], c.cleanarray, header)

        # If you want the mask, here it is :
        #cosmics.tofits("s20160909T093737.7_066LL_sci_mask.fits", c.mask, header)
    
    
def ccdcombine(frames):
    
    for i in xrange(len(frames)):
        command = 'ccdcombine ' + op.join ( frames[i].origloc, frames[i].actionbase["LL"] + frames[i].basename + '_' + frames[i].ifuslot ) + ' ' + frames[i].type   
        run_cure_command( command, 0 )
    
    return command
    

def addphotonnoise(frames, side):

    command = 'addphotonnoise --gain_key GAIN'
    
    filenames = [op.join ( f.origloc, f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in frames]

    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    [f.addbase ('p', side=side) for f in frames] 
    
    return command
 
    
def dividepixelflat(frames, opts, side, uid):

    filenames = [op.join ( f.origloc, f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in frames if f.specid is uid]

    for i in filenames:
        command = 'dividefits %s %s' % (opts,filenames[i])
        run_cure_command( command, 0 )
    
    [f.addbase ('d', side=side) for f in frames if f.specid is uid] 
    
    return command 
    
    
def dividefits(filename, opts):

    command = 'dividefits %s %s' % (opts,filename)
    
    run_cure_command( command, 0 )
        
    return command
    
    
def addfits(filename, opts):

    command = 'addfits %s %s' % (opts,filename)
    
    run_cure_command( command, 0 )
        
    return command 
    

def addkeys(filename):

    command = 'addkeys %s' % (filename)
    
    command = command + ' << EOFK\n EXPTIME=1\n EOFK\n'
    
    run_cure_command( command, 0 )
    
    return command 
    
def headfits(filename,opts):

    command = 'headfits %s %s' % (opts,filename)
        
    run_cure_command( command, 0 )
    
    return command 
    
    
def meanlampfits(side, specid, lamp, dest_dir, basename , opts, frames):
    
    filenames = [op.join ( f.origloc, f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in frames]

    mastername = op.join ( dest_dir , basename + '_' + lamp + '_' + specid + '_' + side + '.fits')
    
    command = 'meanfits %s -o %s' % (opts,mastername)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    return command
    
    
def meantracefits(side, specid, dest_dir, basename , opts, frames):
    
    filenames = [op.join ( f.origloc, f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in frames]

    mastername = op.join ( dest_dir , basename + '_' + specid + '_' + side + '.fits')
    
    command = 'meanfits %s -o %s' % (opts,mastername)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    return command
    

def deformer(mastertrace,masterarc,linesfile,wave_range,ref_line,opts):
    
    command = 'deformer %s -s %s -W %s -o \"%s\" -l %s -a %s %s' % (opts,ref_line,wave_range,redux_dir,linesfile,masterarc,mastertrace)  

    run_cure_command( command, 0 )

    return command
    
    
def subtractsky(frames,side,distmodel,fibermodel,opts,skymaster=""):
    
    for f in frames:
        
        command = 'subtractsky %s %s -d %s -f %s %s' % (opts,skymaster,distmodel,fibermodel,f)  

        run_cure_command( command, 0 )
        
    #[f.addbase('S', side) for f in frames]

    return command


def subtractskyframe(sciframe,skyframe,side,skyscale,distmodel,fibermodel,opts):
    
    scifile = op.join ( redux_dir, sci_dir, sciframe.object, 'pses' + sciframe.basename + '_' + sciframe.ifuslot + '_' + sciframe.type + '_' + side + '.fits') 
    skyfile = op.join ( redux_dir, sci_dir, skyframe.object, 'pses' + skyframe.basename + '_' + skyframe.ifuslot + '_' + skyframe.type + '_' + side + '.fits') 

    skymaster = '-X ' + skyfile 

    command = 'subtractsky %s --x-sky-scaling %s -X %s -D %s -F %s -d %s -f %s %s' % (opts,skyscale,skyfile, distmodel, fibermodel, distmodel,fibermodel,scifile)  

    run_cure_command( command, 0 )
        
    #[f.addbase('S', side) for f in frames]

    return command
    

def fibextract_Resample(frames,distmodel,fibermodel,wave_range,dw,opts):

    for f in frames:
    
        command = 'fiberextract %s -p FeR -W %s -w %s -d %s -f %s %s' %(opts,wave_range,dw,distmodel,fibermodel,f)
        
        run_cure_command( command, 0 )

    return command

def fibextract(frames,distmodel,fibermodel,opts):

    for f in frames:
    
        command = 'fiberextract %s -x -d %s -f %s %s' %(opts,distmodel,fibermodel,f)
        
        run_cure_command( command, 0 )

    return command

def mkcube(ifufile,ditherfile,outname,diffatmref,opts):

    #the default for the DAR correction on CURE is True so adding -d turns off DAR correction
    if diffatmref:
        dar = ""
    else:
        dar = "-d"
    
    command = 'mkcube %s %s -o %s -i %s %s' %(opts,dar,outname,ifufile,ditherfile)
        
    run_cure_command( command, 0 )

    return command

def extend_trace_start(data,start_col,window):
    #defines a region where signal starts to show up - defined by start_col
    first_reg = data[:,start_col:start_col+window]
    #takes the average value of every row in that region
    first_vals = np.average(first_reg,axis=1)

    #extends that average value from the star_col to the x=0 edge of chip
    for i in range(start_col):
        data[:,i] = first_vals

    #returns the chip array with this correction 
    return data          

def initial_setup ( DIR_DICT = None, sci_objects = None, redux_dir = None):
    '''
    Running the initial setup which includes:
    1) Building a standard reduction folder structure
    2) Copying files from specified location into the structure
    3) Checking header to determine if old data procedures are needed (pre LRS2-B swap)
        a)If old LRS2-B sets ucam to old specid (501) else sets ucam to new specid (503)
        b)If reducing old red channel data: copies long exposure Qth frames from lrs2_config to flt for far-red reduction
    4) If reducing red channel data: copies long exposure FeAr frames to cmp for far-red reduction
    5) Creating class variable, VirusFrame, for each file which records
    keywords for other functions as well as time of the observation and basename
    '''

    #########################
    # Build redux directory #
    #########################
    if redux_dir is None:
        print ( "Please provide a reduction directory \"redux_dir\"." )
        return None
    else:
        if not op.exists ( redux_dir ):
            os.mkdir ( redux_dir )
        else:
            if RESTART_FROM_SCRATCH:
                shutil.rmtree ( redux_dir )
                os.mkdir ( redux_dir )

    print ('    ++++++++++++++++++++++++++++++++++++++++++++++++++')
    print ('    + Searching and Sorting All Files in Date Folder +')
    print ('    ++++++++++++++++++++++++++++++++++++++++++++++++++')

    ####################################################
    # Build VIRUS frames for all images in date folder #
    ####################################################
    aframes = [] # will fill this list with VirusFrame class objects for each image
    #make a list of all fits files in the redux directory 
    date_ims = glob.glob(op.join(config.date_folder,'lrs2/lrs2*/exp*/lrs2/*.fits'))
    for f in date_ims:            
        temp, temp1, temp2 = op.basename ( f ).split('_')
        amp                = temp1[3:5]
        if amp == "LL":
            a = VirusFrame( f ) 
            aframes.append(copy.deepcopy(a)) 

    #################
    # Define SPECID #
    #################
    #first frames to pull basic header information 
    a1 = aframes[0]

    #Finds the date of the data taken to know if it is the new or old LRS2-B
    #The swap time is also the time proper calibration data is taken
    #Before the swap time far-red Qth exposures will be taken from config file
    #if LRS2-B finds the date of the data taken to know if it is the new or old LRS2-B - stored as first_run
    data_time      = datetime(int(a1.year), int(a1.month), int(a1.day)) 
    lrs2B_swapTime = datetime(2016, 11, 15)
    if data_time > lrs2B_swapTime:
        first_run = False
    else:
        first_run = True
        print ("WARNING:You are running reduction on shared-risk LRS2_data - calibration data set may not be ideal")

    #if data was taken before 2016, 03, 15 then the only long exposure LDLS flats were taken and are saturated in the orange channel
    #if the data was taken before this data, flats from the orange channel are taken from the config folder
    #if the date is before the calibratio script change second_run = True
    lrs2B_calChange = datetime(2017, 03, 15)
    if data_time > lrs2B_calChange:
        second_run = False
    else:
        second_run = True

    #checks which unit to reduce and sets ucam (specid) - If B decides if old or new specid based of first_run found in intial_setup()
    if config.LRS2_spec == 'R':
        ucam = "502"
    elif config.LRS2_spec == 'B':
        #if LRS2-B finds the date of the data taken to know if it is the new or old LRS2-B
        if first_run:
            ucam = "501"
            print ("You are reducing data with the old LRS2-Blue - You must disregard data from the UV channel")
        else:
            ucam = "503"

    #define cal lamps used based on SPECID   
    if ucam == '501':
        LAMP_DICT = {0:'Hg',1:'Cd'}
        FLT_LAMP = 'ldls'
    if (ucam == '502'):
        FLT_LAMP = 'Qth'
        LAMP_DICT = {0:'Hg',1:'FeAr'}
    if (ucam == '503'):
        FLT_LAMP = 'ldls'
        LAMP_DICT = {0:'Hg',1:'Cd', 2:'FeAr'}

    #This function is for finding the closest folder date to the data taken
    #this is used for finding the best folder to pull special cal data from
    def close_cal_date(fold_list,data_time):
        #makes a list of all of the names of all the date folders
        long_folds = glob.glob(fold_list + '/*')
        folds = []
        #creates datetime object for each date folder 
        for l in long_folds:
            year  = int(l.split('/')[-1][0:4])
            month = int(l.split('/')[-1][4:6])
            day   = int(l.split('/')[-1][6:8])
            folds.append(datetime(year,month,day))
        #finds the closest date/folder to the date the data was taken 
        close_date      = min(folds, key=lambda x:abs(x-data_time))
        close_date_name = str(close_date.year).zfill(4)+str(close_date.month).zfill(2)+str(close_date.day).zfill(2)

        #returns a list of all of the files names of the files in that date folder 
        long_files = glob.glob(fold_list+'/'+close_date_name+'/exp*/lrs2/*.fits')
        return long_files


    ############################################
    # Find calibration data and science images #
    ############################################
    #from the aframes makes lists of files vframes according to type and specid (ucam)
    tframes  = [a for a in aframes if (a.specid == ucam) and (a.cal_side == None or a.cal_side == config.LRS2_spec)] # gives all frames 

    #------------#
    # zro frames #
    #------------#
    zframes_orig   = [t for t in tframes if t.type == "zro" ] # gives just "zro" frames

    if len(zframes_orig) == 0:
        sys.exit("No biases were found for this night")
    print ('Found '+str(len(zframes_orig))+' zro frames')

    #------------#
    # flt frames #
    #------------#

    #if old first run data and LRS2-R - need to use long Qth exposures in config
    if ucam == '502':
        if first_run:
            fframes_orig = []
            longQthR_folds  = config.configdir+'/longExpCals/long_Qth_R'
            longQthR_files = close_cal_date(longQthR_folds,data_time)

            num = 0
            for f in longQthR_files:            
                temp, temp1, temp2 = op.basename ( f ).split('_')
                amp                = temp1[3:5]
                if amp == "LL":
                    a = VirusFrame( f ) 
                    if a.specid == ucam:
                        fframes_orig.append(copy.deepcopy(a)) 
                        num = num + 1

            print ('Including '+str(num)+' long exposure Qth flats for far-red channel reduction')

        else:
            fframes_orig   = [t for t in tframes if t.type == "flt" and t.object == 'Qth'] # gives just "flt" frames

    #if old second run data and LRS2-B - need to use short LDLS exposures in config for orange channel
    elif (ucam == '501') or (ucam == '503' and second_run):
        fframes_orig = []
        shortLDLS_folds  = config.configdir+'/short_OrgFlts'
        shortLDLS_files = close_cal_date(shortLDLS_folds,data_time)

        num = 0
        for f in shortLDLS_files:            
            temp, temp1, temp2 = op.basename ( f ).split('_')
            amp                = temp1[3:5]
            if amp == "LL":
                a = VirusFrame( f ) 
                if a.specid == ('501' or '503'):
                    fframes_orig.append(copy.deepcopy(a)) 
                    num = num + 1

        print ('Including '+str(num)+' short exposure LDLS flats for orange channel reduction')

    else:
        fframes_orig   = [t for t in tframes if t.type == "flt" and t.object == 'ldls_short'] # gives just "flt" frames

    print ('Found '+str(len(fframes_orig))+' '+FLT_LAMP+' flt frames')

    if len(fframes_orig) == 0:
        sys.exit("No "+FLT_LAMP+" flat lamp exposures were found for this night")

    #-----------#
    # Hg frames #
    #-----------#
    hgframes  = [t for t in tframes if t.type == "cmp" and t.object == "Hg"] # gives just "Hg" frames

    if len(hgframes) == 0:
        sys.exit("No Hg lamp exposures were found for this night")
    print ('Found '+str(len(hgframes))+' Hg frames')

    #----------------#
    # 501 cmp frames #
    #----------------#
    if ucam == '501':
        cdframes = [t for t in tframes if t.type == "cmp" and t.object == "Cd"] # gives just "Cd" frames

        if len(cdframes) == 0:
            sys.exit("No Cd lamp exposures were found for this night")
        print ('Found '+str(len(cdframes))+' Cd frames')

        lframes_orig  = hgframes + cdframes # gives just "cmp" frames

    #----------------#
    # 502 cmp frames #
    #----------------#
    if ucam == '502':
        faframes = [t for t in tframes if t.type == "cmp" and t.object == "FeAr"] # gives just "FeAr" frames
         
        print ('Found '+str(len(faframes))+' FeAr frames')

        #if LRS2-R need to include long FeAr cmps for far-red channel
        longFeArR_folds  = config.configdir+'/longExpCals/long_FeAr_R'
        longFeArR_files = close_cal_date(longFeArR_folds,data_time)

        num = 0
        for f in longFeArR_files:            
            temp, temp1, temp2 = op.basename ( f ).split('_')
            amp                = temp1[3:5]
            if amp == "LL":
                a = VirusFrame( op.join( f )) 
                if a.specid == ucam:
                    faframes.append(copy.deepcopy(a)) 
                    num = num + 1

        print ('Including '+str(num)+' long exposure FeAr comps for far-red channel reduction')

        if len(faframes) == 0:
            sys.exit("No FeAr lamp exposures were found for this night")

        lframes_orig  = hgframes + faframes # gives just "cmp" frames

    #----------------#
    # 503 cmp frames #
    #----------------#
    if ucam == '503':
        cdframes = [t for t in tframes if t.type == "cmp" and t.object == "Cd"] # gives just "Cd" frames

        if len(cdframes) == 0:
            sys.exit("No Cd lamp exposures were found for this night")
        print ('Found '+str(len(cdframes))+' Cd frames')

        faframes = [t for t in tframes if t.type == "cmp" and t.object == "FeAr"] # gives just "FeAr" frames

        print ('Found '+str(len(faframes))+' FeAr frames')

        #if LRS2-B need to include long FeAr cmps for UV channel
        longFeArB_folds  = config.configdir+'/longExpCals/long_FeAr_B'
        longFeArB_files = close_cal_date(longFeArB_folds,data_time)

        num = 0
        for f in longFeArB_files:            
            temp, temp1, temp2 = op.basename ( f ).split('_')
            amp                = temp1[3:5]
            if amp == "LL":
                a = VirusFrame( op.join( f )) 
                if a.specid == ucam:
                    faframes.append(copy.deepcopy(a))
                    num = num + 1

        print ('Including '+str(num)+' long exposure FeAr comps for UV channel reduction')

        if len(faframes) == 0:
            sys.exit("No FeAr lamp exposures were found for this night")

        lframes_orig  = hgframes + cdframes + faframes # gives just "cmp" frames

    #------------#
    # drk frames #
    #------------#
    if config.subDarks:
        dframes_orig  = [t for t in tframes if t.type == "drk" ] # gives dark frames
        print ('Found '+str(len(dframes_orig))+' drk frames')
        drk_exptime = list(set([float(d.exptime) for d in dframes_orig]))
    else:
        dframes_orig  = []

    #------------#
    # sci frames #
    #------------#
    allsframes  = [a for a in aframes if a.type == "sci" and (a.specid == ucam)]
    sci_obj_names = [s.object for s in allsframes]
    sci_obj_names = list(set(sci_obj_names))

    if len(config.sci_objects) == 0:
        config.sci_objects = sci_obj_names
        sys.exit("No science objects choosen. You must choose science objects to be reduced")

    sframes_lis = []
    spframes_lis = []
    for s in config.sci_objects:
        spfr = [t for t in tframes if t.type == "sci" and t.object == s ]
        sfr  = [a for a in aframes if a.type == "sci" and (a.specid == ucam) and a.object == s]
        spframes_lis.append(spfr)
        sframes_lis.append(sfr)
        print ("There were "+str(len(sfr))+" science frames found for "+s)
        if len(sfr) == 0:
            print ("There were no science frames with object name: "+s)
            sys.exit("These are the object names found for "+op.basename(config.date_folder)+": "+str(sci_obj_names))

    spframes_orig = [j for i in spframes_lis for j in i] # gives just "sci" frames with correct config.LRS2_spec pointing
    sframes_orig  = [j for i in sframes_lis  for j in i] # gives just "sci" frames with any pointing
    sci_exptime = list(set([float(s.exptime) for s in sframes_orig])) #finds exposure time for all the sci frames

    #Finds sky frames if the use choose to use sky frames in the config file
    if config.use_sky_frames:  
        if config.LRS2_spec == 'R':
            sky_side = 'B'
        else:
            sky_side = 'R'

        #automatically searches for sky frames if the user did not give paths to sky frames
        if len(config.user_skyframes) == 0:
            print("Automatically searching for sky frames")
            if first_run == True:
                print ("Automatic search for sky frames cannot be done with shared-risk LRS2_data")
                sys.exit("You must give sky frame paths in user_skyframes in lrs2_config.py")

            #looks for frames with same ucam as specid but the objects names show pointing to other side. These will be sky frames
            allskyframes = [a for a in aframes if a.type == "sci" and (a.specid == ucam) and (a.cal_side == sky_side)]
            #Check that sky frames were found 
            if len(allskyframes) == 0:
                sys.exit("NO SKY FRAMES FOUND: enter a path to a sky frame in user_skyframes in lrs2_config.py")

            skyexptime = [a.exptime for a in allskyframes] #finds exposure time for all of the sky frames

            skyframe_index = []
            #loops through and finds the sky frame with the clostest exposure time to each science frame expsoure time
            for s in sci_exptime:
                closest_index = min(range(len(skyexptime)), key=lambda i: abs(skyexptime[i]-s))
                skyframe_index.append(closest_index)

            #remove duplicate incdecies to make sure there are not duplciate images in the skyframe list
            skyframe_index = list(set(skyframe_index))

            #appends the skyframe and the skyframe object name to lists (excluding duplicates)
            skyframes_orig = []
            for s in skyframe_index:
                skyframe = allskyframes[s]
                skyframes_orig.append(skyframe)

        #if they user provided sky frame paths use these sky frames instead of searching for frames
        elif len(config.user_skyframes) > 0:
            print("Finding user's sky frames")
            #need find the files corresponding to the observation folders the user provided 
            user_sky_list = []
            for s in config.user_skyframes:
                user_sky_list.extend(glob.glob(op.join(s,'lrs2/*.fits')))

            #Check that sky frames were found in the user provided path    
            if len(user_sky_list) == 0:
                sys.exit("NO SKY FRAMES FOUND: Check your the path you provided. It needs to point to an exposure folder.")

            #build virus frames for each sky frame found 
            skyframes_first = [] 
            for f in user_sky_list:    
                temp, temp1, temp2 = op.basename ( f ).split('_')
                amp                = temp1[3:5]
                if amp == "LL":
                    a = VirusFrame( f ) 
                    skyframes_first.append(copy.deepcopy(a)) 

            skyframes_orig = [a for a in skyframes_first if a.type == "sci" and (a.specid == ucam) and (a.cal_side == sky_side or a.cal_side == None)]
            
        skyframe_names  = [(a.basename + '_' + a.ifuslot + '_' + a.type) for a in skyframes_orig]
        skyframe_objs  = [a.object for a in skyframes_orig]
        skytimes       = [s.exptime for s in skyframes_orig]

        print ("There were "+str(len(skyframes_orig))+" sky frames found")
        print ("    Sky frames found: "+str(skyframe_names))

        #defines the frames that are just the science frames before the sky frames are added to the sci list
        only_sframes = sframes_orig

        #now the sky frames are added to the science frames for reduction
        sframes_orig = sframes_orig + skyframes_orig
        #The sky objects are added to the sci object list so the files are sorted properly 
        config.sci_objects = config.sci_objects + skyframe_objs
        #make sure there are not duplicate object names 
        config.sci_objects = list(set(config.sci_objects))

    else:
        #because sky frames are not added this is the same as just the sci frames list
        #variable made because it has to be returned in the function
        only_sframes = sframes_orig
        sky_side = None
        skytimes = None
        skyframes_orig = None

    #Check that data is correct
    if len(spframes_orig) == 0:
        print ("WARNING: Science frames were found for your science objects but not with LRS2-"+config.LRS2_spec+" pointings - these may just be sky frames")

    if all_copy:
        print ('    +++++++++++++++++++++++++++++++++++++++++++')
        print ('    + Copying Frames Into Directory Structure +')
        print ('    +++++++++++++++++++++++++++++++++++++++++++')
    #########################################
    # Copying frames to directory structure #
    #########################################
    vframes = []
    #dictionary of data type folders 
    file_loc_dir = [  zframes_orig,   dframes_orig,   lframes_orig,   fframes_orig,   sframes_orig ] # Should match order of dictionary, otherwise there will be mismatches

    # Loop through the file location directories     
    for i in xrange ( len ( file_loc_dir ) ):
        # If the file_loc_dir[i] is None, then nothing is done
        if file_loc_dir[i] is not None:
            # If the reduction location exists, don't re-make the directory (also, any files in that directory remain)
            if not op.exists ( op.join ( redux_dir, DIR_DICT[i] ) ):
                os.mkdir ( op.join ( redux_dir, DIR_DICT[i] ) )

            # Loop through the retrieved files names to copy to new structure
            # Create a VirusFrame class for each frame that can maintain info for each original frame
            # The VirusFrame is critical for the rest of the reduction pipeline
            # Only gather VirusFrame objects for LL frames as a single VirusFrame serves for all four amps
            if all_copy:
                #must copy all file to directories first
                for f in file_loc_dir[i]:  
                    for amp in SPEC: 
                        shutil.copy( op.join ( f.origloc, f.actionbase[amp] + f.basename + '_' + f.ifuslot + amp + '_' + f.type + '.fits'), op.join ( redux_dir, DIR_DICT[i] ) )
                for f in file_loc_dir[i]:        
                    a = VirusFrame( op.join( redux_dir, DIR_DICT[i], op.basename ( f.filename ) ) ) 
                    vframes.append(copy.deepcopy(a))
            else: 
                for f in file_loc_dir[i]:
                    a = VirusFrame( f.filename ) 
                    vframes.append(copy.deepcopy(a))
                        
    return vframes, first_run, second_run, ucam, LAMP_DICT, FLT_LAMP, sky_side, only_sframes, skyframes_orig, skytimes


def basicred(DIR_DICT, sci_objects, redux_dir, basic = False, dividepf = False,
              normalize = False, masterdark = False, masterarc = False, mastertrace = False):
    '''
    Running the basic reduction which includes:
    1) Overscan subtract and trim all frames
    2) Create blank error frame with readnoise in it
    3) Remove cosmic rays from sci frames (user option)
    4) Create master bias frame for each amp 
    5) Subtract master bias from cmp/flt/sci frames
    6) Create master dark for each amp make scaled frame for each sci exposure time (user option)
    7) Subtract master dark from each science frame (user option)
    8) Ccdcombine all frames which puts units in e-
    9) Add photon noise to combined frames
    10) Divide pixelflat from cmp/flt/sci frames (user option)
    11) Normalize cmps and flts 
    12) Combine cmps and flts into masterarc and mastertrace
    '''

    print ('*************************')
    print ('* BUILDING IMAGE FRAMES *')
    print ('*************************')

    #holds the VIRUS frames for all of the data 
    vframes, first_run, second_run, ucam, LAMP_DICT, FLT_LAMP, sky_side, only_sframes, skyframes_orig, skytimes = initial_setup ( DIR_DICT, config.sci_objects, redux_dir )

    #inital reference frame
    f1 = vframes[0]

    zframes  = [v for v in vframes if v.type == "zro" ] # gives just "zro" frames
    dframes  = [v for v in vframes if v.type == "drk" ] # gives just "drk" frames
    lframes  = [v for v in vframes if v.type == "cmp" ] # gives just "cmp" frames
    fframes  = [v for v in vframes if v.type == "flt" ] # gives just "flt" frames
    sframes  = [v for v in vframes if v.type == "sci" ] # gives just "sci" frames
    cframes  = fframes + lframes # gives "flt" and "cmp" frames
    oframes  = cframes + sframes + dframes # gives "flt", "drk", "cmp", and "sci" frames (basically just not "zro")

    #make a copy of the lsr2_config file to be added to your directory
    #if the file exists - remove the file and replace it.
    if os.path.isfile(redux_dir+'/lrs2_config_'+redux_dir.split('/')[-1]+'_copy.py') == True:
        os.remove(redux_dir+'/lrs2_config_'+redux_dir.split('/')[-1]+'_copy.py')
        shutil.copy ( os.path.dirname(os.path.realpath(__file__))+'/lrs2_config.py', redux_dir+'/lrs2_config_'+redux_dir.split('/')[-1]+'_copy.py' )
    else:
        shutil.copy ( os.path.dirname(os.path.realpath(__file__))+'/lrs2_config.py', redux_dir+'/lrs2_config_'+redux_dir.split('/')[-1]+'_copy.py' )

    # Run basic reduction
    if config.basic:
        for sp in SPEC:
            trimsec = f1.trimsec # Trimsec assumed to be the same for all frames of a given amp
            biassec = f1.biassec # Biassec assumed to be the same for all frames of a given amp
            print ('**************************')
            print ('* MAKE ERROR FRAME FOR '+sp+' *')
            print ('**************************')
            mkerrorframe ( vframes, sp )               # for all frames

            print ('***************************')
            print ('* SUBTRACT OVERSCAN FOR '+sp+' *')
            print ('***************************')
            subtractoverscan( biassec, vframes, sp )   # for all frames

            print ('*****************************')
            print ('* EXTRACT DATA REGION FOR '+sp+' *')
            print ('*****************************')
            extractfits ( trimsec, vframes, sp )       # for all frames

            # Remove cosmic rays using L.A.cosmic
            if rmcosmics:
                print ('******************************************')
                print ('* REMOVE COSMIC RAYS (SCI IMAGES) FOR '+sp+' *')
                print ('******************************************')
                #if reducing old LRS2-Blue data: do not run rmcosmics on UV data(LL+LU), only orange(RL+RU)
                if (config.LRS2_spec == 'B') and (first_run):
                    if (sp == 'RL') or (sp == 'RU'):
                        rmcosmicfits ( sframes, sp )       # for sci frames - because this is slow
                else:
                    rmcosmicfits ( sframes, sp )       # for sci frames - because this is slow                    
            
            print ('**************************')
            print ('* BUILD MASTERBIAS FOR '+sp+' *')
            print ('**************************')
            meanbiasfits   ( sp, ucam, redux_dir, 'masterbias', meanfitsopts, zframes ) # meanfits for masterbias for unique specid 

            print ('*****************************')
            print ('* SUBTRACT MASTERBIAS FOR '+sp+' *')
            print ('*****************************')
            masterbiasname = op.join ( redux_dir, 'masterbias' + '_' + ucam + '_' + sp + '.fits' ) 
            subtractbias   ( oframes, masterbiasname, sp) # for all frames

            # Create Master Dark Frames
            if masterdark:
                print ('********************************************')
                print ('* BUILDING MASTERDARK and SUBTRACTING FROM SCI FRAMES FOR '+sp+' *')
                print ('********************************************')
                for e in sci_exptime:
                    close_drk = min(drk_exptime, key=lambda x:abs(x-e))
                    dframesselect = [d for d in dframes if d.exptime == close_drk] 
                    meandarkfits ( sp, ucam, redux_dir, 'masterdark_'+close_drk+'sec', meanfitsopts, dframesselect ) # meanfits for masterdark for unique specid
               
                    #Subtracts Master Dark from Science Frames  
                    scale_drk = e/close_drk
                    sframesselect = [s for s in sframes if s.exptime == e] 
                    masterdarkname = op.join ( redux_dir, 'masterdark_'+close_drk+'sec' + '_' + ucam + '_' + sp + '.fits' )
                    subtractdark ( sframesselect, masterdarkname, scale_drk, sp) # for sci frames 

                
        print ('*******************************************')
        print ('* COMBINING CCD AMPS - MULTIPYING BY GAIN *')
        print ('*******************************************')

        ccdcombine ( oframes ) # Combine amps and multiply by gain for all frames other than zro's
        for o in oframes:
            o.actionbase["L"] = o.actionbase["LL"]
            o.actionbase["R"] = o.actionbase["RL"]
            if o.clean:
                for sp in SPEC:
                    filename   = op.join ( o.origloc,        o.actionbase[sp] + o.basename + '_' + o.ifuslot + sp + '_' + o.type + '.fits')
                    filename_e = op.join ( o.origloc, 'e.' + o.actionbase[sp] + o.basename + '_' + o.ifuslot + sp + '_' + o.type + '.fits')
                    os.remove ( filename )
                    os.remove ( filename_e )
        for side in SPECBIG:    
            addphotonnoise ( oframes , side )

    # Dividing by Pixel Flat
    if dividepf:
        print ('*************************')
        print ('* DIVDING BY PIXEL FLAT *')
        print ('*************************')
        for side in SPECBIG:
            pflat = op.join( pixflatdir, "pixelflat_cam{:03d}_{:s}.fits".format( ucam, side.upper() ) )
            opt = "--file {:s}".format(pflat)
            dividepixelflat(oframes, opt, side, ucam)
    
    # Normalizing Calibration Frames by Exposure Time
    if normalize:
        print ('***************************************************')
        print ('* NORMALIZING CALIBRATION IMAGES BY EXPOSURE TIME *')
        print ('***************************************************')
        for cal in cframes:
            for side in SPECBIG:
                opt     = '-c {:0.1f}'.format(cal.exptime)
                filename = op.join ( cal.origloc, cal.actionbase[side] + cal.basename + '_' + cal.ifuslot + '_' + cal.type + '_' + side + '.fits' )
                dividefits ( filename, opt )
                cal.addbase ('d', side = side)
                filename = op.join ( cal.origloc, cal.actionbase[side] + cal.basename + '_' + cal.ifuslot + '_' + cal.type + '_' + side + '.fits' )
                headfits ( filename , headfitsopts )
    
    # Make Master Arc Frames           
    if masterarc:        
        print ('**********************')
        print ('* BUILDING MASTERARC *')
        print ('**********************')
        for side in SPECBIG:
            for lamp in LAMP_DICT.values():
                #Creates a masterarc frame for each arc lamp in LAMP_DICT
                lframesselect = [l for l in lframes if l.object == lamp]
                #If more than one image for that lamp take median image 
                if len(lframesselect)>1:
                    meanlampfits(side, ucam, lamp, redux_dir, 'masterarc' , arcopts, lframesselect) 
                #If only one frame for the lamp make that frame the master frame for that lamp 
                elif len(lframesselect) == 1:
                    filename = [op.join ( f.origloc, f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in lframesselect]
                    mastername = op.join ( redux_dir , 'masterarc' + '_' + lamp + '_' + ucam + '_' + side + '.fits')
                    shutil.copy ( filename[0], mastername )
                    efilename = [op.join ( f.origloc, 'e.' + f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in lframesselect]
                    emastername = op.join ( redux_dir , 'e.masterarc' + '_' + lamp + '_' + ucam + '_' + side + '.fits')
                    shutil.copy ( efilename[0], emastername )
                #If there are zero frames found then the arc lamp information provided is wrong
                else:
                    print ("You did not provide the correct arc lamp data")
                    sys.exit( "For LRS2-"+config.LRS2_spec+" You must provide "+lamp+" data")

            #Combine each arc master frame for each lamp in LAMP_DICT into one masterarc
            #for the UV channel you need to add the Hg and FeAr (the second two) lamps together 
            if ucam == '503' and side == 'L':
                opt = "--file {:s}".format( op.join ( redux_dir, 'masterarc' + '_' + LAMP_DICT[0] + '_' + ucam + '_' + side + '.fits' ) )
                filename = op.join ( redux_dir, 'masterarc' + '_' + LAMP_DICT[2] + '_' + ucam + '_' + side + '.fits' )
                addfits ( filename, opt)
                shutil.copy( op.join ( redux_dir, 'smasterarc' + '_' + LAMP_DICT[2] + '_' + ucam + '_' + side + '.fits' ), 
                             op.join ( redux_dir, 'masterarc' + '_' + ucam + '_' + side + '.fits' ) )
                shutil.copy( op.join ( redux_dir, 'e.smasterarc' + '_' + LAMP_DICT[2] + '_' + ucam + '_' + side + '.fits' ), 
                             op.join ( redux_dir, 'e.masterarc' + '_' + ucam + '_' + side + '.fits' ) )
                os.remove  ( op.join ( redux_dir, 'smasterarc' + '_' + LAMP_DICT[2] + '_' + ucam + '_' + side + '.fits' ) )

            #for all other channels you just need to add two lamps together (for orange it add Cd and Hg (the first two lamps) together)
            else: 
                opt = "--file {:s}".format( op.join ( redux_dir, 'masterarc' + '_' + LAMP_DICT[0] + '_' + ucam + '_' + side + '.fits' ) )
                filename = op.join ( redux_dir, 'masterarc' + '_' + LAMP_DICT[1] + '_' + ucam + '_' + side + '.fits' )
                addfits ( filename, opt)
                shutil.copy( op.join ( redux_dir, 'smasterarc' + '_' + LAMP_DICT[1] + '_' + ucam + '_' + side + '.fits' ), 
                             op.join ( redux_dir, 'masterarc' + '_' + ucam + '_' + side + '.fits' ) )
                shutil.copy( op.join ( redux_dir, 'e.smasterarc' + '_' + LAMP_DICT[1] + '_' + ucam + '_' + side + '.fits' ), 
                             op.join ( redux_dir, 'e.masterarc' + '_' + ucam + '_' + side + '.fits' ) )
                os.remove  ( op.join ( redux_dir, 'smasterarc' + '_' + LAMP_DICT[1] + '_' + ucam + '_' + side + '.fits' ) )

        #clean intermediate frames in cmp if clean
        for l in lframes:
            if l.clean:
                for side in SPECBIG:
                    filename   = op.join ( l.origloc,        l.actionbase[side] + l.basename + '_' + l.ifuslot + '_' + l.type + '_' + side + '.fits')
                    filename_e = op.join ( l.origloc, 'e.' + l.actionbase[side] + l.basename + '_' + l.ifuslot + '_' + l.type + '_' + side + '.fits')
                    os.remove ( filename )
                    os.remove ( filename_e )
         
    # Make Master Trace Frames
    if mastertrace:
        print ('************************')
        print ('* BUILDING MASTERTRACE *')
        print ('************************')
        for side in SPECBIG:
            #If there is more than one frame take a median from of the ones provided for the mastertrace
            if len ( fframes ) > 1:
                meantracefits(side, ucam, redux_dir, 'mastertrace', traceopts, fframes)
            #If there was only one frame provided that frame becomes the mastertrace
            else:
                filename = [op.join ( f.origloc, f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in fframes]
                mastername = op.join ( redux_dir , 'mastertrace' + '_' + ucam + '_' + side + '.fits')
                shutil.copy ( filename[0], mastername )
                efilename = [op.join ( f.origloc, 'e.' + f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in fframes]
                emastername = op.join ( redux_dir , 'e.mastertrace' + '_' + ucam + '_' + side + '.fits')
                shutil.copy ( efilename[0], emastername )

        #Clean intermediate frames in flt if clean        
        for f in fframes:
            if f.clean:
                for side in SPECBIG:
                    filename   = op.join ( f.origloc,        f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits')
                    filename_e = op.join ( f.origloc, 'e.' + f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits')
                    os.remove ( filename )
                    os.remove ( filename_e )

    #Extend traces to detector edges for the far red channel. 
    #This will create a new mastertrace_502_R.fits and save the old one as mastertrace_502_R_orig.fits
    if fix_chan and (config.LRS2_spec == 'R'):
        print ('**************************************')
        print ('* FIXING FAR RED CHANNEL MASTERTRACE *')
        print ('**************************************')
        #X value to start with to extend trace to the edge of chip
        start_col = 50

        im = pyfits.open(redux_dir+'/mastertrace_'+str(ucam)+'_R.fits')
        hdr = pyfits.getheader(redux_dir+'/mastertrace_'+str(ucam)+'_R.fits')
        dat =  im[0].data

        dat = extend_trace_start(dat,start_col,10)

        #before writing the new file it moves the old mastertrace to a file with _orig appended to save original file
        shutil.move(redux_dir+'/mastertrace_'+str(ucam)+'_R.fits',redux_dir+'/mastertrace_'+str(ucam)+'_R_orig.fits')
        pyfits.writeto(redux_dir+'/mastertrace_'+str(ucam)+'_R.fits', data = dat, header = hdr, clobber=True)

        #fix far-red channel masterarc

    #Extend traces to detector edges for the orange channel. 
    #This will create a new mastertrace_501_R.fits and save the old one as mastertrace_501_R_orig.fits
    if fix_chan and (config.LRS2_spec == 'B'):
        print ('*************************************')
        print ('* FIXING ORANGE CHANNEL MASTERTRACE *')
        print ('*************************************')
        #X value to start with to extend trace to the edge of chip
        start_col = 100

        im = pyfits.open(redux_dir+'/mastertrace_'+str(ucam)+'_R.fits')
        hdr = pyfits.getheader(redux_dir+'/mastertrace_'+str(ucam)+'_R.fits')
        dat =  im[0].data

        dat = extend_trace_start(dat,start_col,10)

        #before writing the new file it moves the old mastertrace to a file with _orig appended to save original file
        shutil.move(redux_dir+'/mastertrace_'+str(ucam)+'_R.fits',redux_dir+'/mastertrace_'+str(ucam)+'_R_orig.fits')
        pyfits.writeto(redux_dir+'/mastertrace_'+str(ucam)+'_R.fits', data = dat, header = hdr, clobber=True)

    if sort_sci:
        print ('**************************************************')
        print ('* SORTING SCIENCE FRAMES INTO OBJECT DIRECTORIES *')
        print ('**************************************************') 
        #Make a directory in sci folder for each science object    
        for s in config.sci_objects:
            os.mkdir ( op.join( redux_dir, sci_dir, s ))
        for s in sframes:
            for side in SPECBIG:
                sname = 'pses'+s.basename+'_'+s.ifuslot+'_'+s.type+'_'+side+'.fits'
                shutil.move (op.join(s.origloc, sname ), op.join( s.origloc, s.object ) )
                shutil.move (op.join(s.origloc, 'e.'+sname ), op.join( s.origloc, s.object ) )            

    # Run Deformer
    if config.run_deformer:
        print ('*************************************************************************')
        print ('* RUNNING DEFORMER TO BUILD DISTORTION SOLUTION AND WAVELENGTH SOLUTION *')
        print ('*************************************************************************')
        #check that basic has been run 
        trace_files = glob.glob(op.join(redux_dir,'mastertrace*'))
        arc_files   = glob.glob(op.join(redux_dir,'masterarc*'))
        if len(trace_files) == 0 or len(arc_files) == 0:
            sys.exit("You must run basic reduction before you can run deformer")

        for side in SPECBIG:  
            #selects wavelength range and ref arc line for each channel
            if (config.LRS2_spec == 'B') and (side == 'L'):
                wave_range = '[3600,4700]'
                ref_line = 6
            if (config.LRS2_spec == 'B') and (side == 'R'):
                wave_range = '[4600,7000]'
                ref_line = 7
            if (config.LRS2_spec == 'R') and (side == 'L'):
                wave_range = '[6500,8500]'
                ref_line = 3
            if (config.LRS2_spec == 'R') and (side == 'R'):
                wave_range = '[8000,10500]'
                ref_line = 1
            #copy the lines file used to this directory 
            shutil.copy ( op.join(linesdir,'lines' + '_' + side + '_' + ucam +'.par'), op.join(redux_dir,'lines' + '_' + side + '_' + ucam +'.par' ) )
            #build the names of the files given to deformer 
            mastertrace = op.join ( redux_dir, 'mastertrace' + '_' + ucam + '_' + side + '.fits' )
            masterarc   = op.join ( redux_dir, 'masterarc' + '_' + ucam + '_' + side + '.fits' )
            linefile    = op.join ( redux_dir, 'lines' + '_' + side + '_' + ucam +'.par' )
            deformer ( mastertrace, masterarc, linefile, wave_range, ref_line, deformeropts)
    
    # Run sky subtraction            
    if config.subsky:  
        print ('************************************************')
        print ('* PERFORMING SKY SUBTRACTION ON SCIENCE FRAMES *')
        print ('************************************************')
        #check that deformer has been run 
        dist_files = glob.glob(op.join(redux_dir,'*.dist'))
        if len(dist_files) == 0:
            sys.exit("You must run deformer before you can run sky subtraction")

        if config.use_sky_frames:
            print ('    +++++++++++++++++++')
            print ('    + Using Sky Frame +')
            print ('    +++++++++++++++++++')
            for side in SPECBIG:            
                #Sfiles = glob.glob(redux_dir + "/" + sci_dir + "/*/" + "p*_sci_"+side+".fits")
                #skyframes = [s for s in sframes if s.cal_side == sky_side]
                distmodel = op.join ( redux_dir, 'mastertrace' + '_' + ucam + '_' + side + '.dist' )
                fibermodel = op.join ( redux_dir, 'mastertrace' + '_' + ucam + '_' + side + '.fmod' )
                for s in only_sframes:
                    #find the sky frame with the closest exposure time 
                    closest_index = min(range(len(skytimes)), key=lambda i: abs(float(skytimes[i])-float(s.exptime)))
                    skyframe = skyframes_orig[closest_index]
                    #define the scale factor to scale up sky exposure to to deal with different expsosure times between sky and sci frames
                    #scale factor for the sky frame is the sci exposure time divided by they sky frame exposure time
                    #this scales by exposure time and additionally by a factor the user chooses 
                    skyscale = float(s.exptime)/float(skytimes[closest_index]) * config.sky_scaling
                    subtractskyframe(s,skyframe,side,skyscale,distmodel,fibermodel,subskyopts)

        else:
            for side in SPECBIG:
                distmodel = op.join ( redux_dir, 'mastertrace' + '_' + ucam + '_' + side + '.dist' )
                fibermodel = op.join ( redux_dir, 'mastertrace' + '_' + ucam + '_' + side + '.fmod' )
                Sfiles = glob.glob(redux_dir + "/" + sci_dir + "/*/" + "p*_sci_"+side+".fits")
                subtractsky(Sfiles,side,distmodel,fibermodel,subskyopts)

    # Run fiberextract
    if config.fiberextract:  
        print ('****************************************')
        print ('* EXTRACTING SPECTRA IN SCIENCE FRAMES *')
        print ('****************************************')
        #check that deformer has been run 
        dist_files = glob.glob(op.join(redux_dir,'*.dist'))
        if len(dist_files) == 0:
            sys.exit("You must run deformer before you can run fiber extract")

        if config.wl_resample:
            for side in SPECBIG:
                #for each channel selects correct wavlength range and dw
                if (config.LRS2_spec == 'B') and (side == 'L'):
                    wave_range = '3643,4658'
                    dw = '0.49'
                if (config.LRS2_spec == 'B') and (side == 'R'):
                    wave_range = '4600,7000'
                    dw = '1.2'
                if (config.LRS2_spec == 'R') and (side == 'L'):
                    wave_range = '6432,8451'
                    dw = '0.978'
                if (config.LRS2_spec == 'R') and (side == 'R'):
                    wave_range = '8324,10565'
                    dw = '1.129'

                #finds if there are sky subtracted files. If so it uses those.
                Sfiles = glob.glob(redux_dir + "/" + sci_dir + "/*/" + "Sp*_sci_"+side+".fits")
                if len(Sfiles) == 0:
                    Sfiles = glob.glob(redux_dir + "/" + sci_dir + "/*/" + "p*_sci_"+side+".fits")

                distmodel = redux_dir + "/mastertrace_" + ucam + "_" + side + ".dist"
                fibermodel = redux_dir + "/mastertrace_" + ucam + "_" + side + ".fmod"
                fibextract_Resample(Sfiles,distmodel,fibermodel,wave_range,dw,fibextractopts) 
        else:
            print ('    +++++++++++++++++++++++++++++++++')
            print ('    + Extraction Without Resampling +')
            print ('    +++++++++++++++++++++++++++++++++')
            for side in SPECBIG:   
                #finds if there are sky subtracted files. If so it uses those.
                Sfiles = glob.glob(redux_dir + "/" + sci_dir + "/*/" + "Sp*_sci_"+side+".fits")
                if len(Sfiles) == 0:
                    Sfiles = glob.glob(redux_dir + "/" + sci_dir + "/*/" + "p*_sci_"+side+".fits")

                distmodel = redux_dir + "/mastertrace_" + ucam + "_" + side + ".dist"
                fibermodel = redux_dir + "/mastertrace_" + ucam + "_" + side + ".fmod"
                fibextract(Sfiles,distmodel,fibermodel,fibextractopts)

    #CURE saves these files from deformer outside of the redux directory for some reason.
    #This moves them inside of the redux directory.
    left_files = glob.glob('*.log') + glob.glob('*.residuals')
    if len(left_files) > 0:
        for l in left_files:
            shutil.move(l, op.join(redux_dir,l))

    # Run mkcube
    if config.makecube:
        print ('***********************')
        print ('* BUILDING DATA CUBES *')
        print ('***********************')
        for s in config.sci_objects:
        #cd inside of the science directory 
            location_prefix = redux_dir + "/" + sci_dir + "/" + s + "/"
            os.chdir(location_prefix)

            #builds a list of files to build a data cube 
            #checks for if there are wavelength resampled fiber extracted and sky subtracted files
            Fefiles = glob.glob( "FeRS*_sci_*.fits")
            #if not then checks for fiber extraced and sky subtracted files 
            if len(Fefiles) == 0:
                Fefiles = glob.glob("FeS*_sci_*.fits")
            #if not checks for any type of fiber extracted files (resampled or not)
            if len(Fefiles) == 0:
                Fefiles = glob.glob("Fe*_sci_*.fits")

            if len(Fefiles) == 0:
                sys.exit("No fiber extracted files found - You must run fiberextract before building data cube")

            for f in Fefiles:
                im  = pyfits.open(f)
                hdr = im[0].header
                #extracting header information for to know what channel file corresponds to and for dither file information 
                uca     = hdr['SPECID']
                side    = hdr['CCDPOS']
                airmass = hdr['AIRMASS']

                ditherfile = 'dither_LRS2_' + side + '.txt'

                #choose correct mapping file for the channel you are using
                if   (config.LRS2_spec == 'B') and (side == 'L'):
                    IFUfile = mapdir+'LRS2_B_UV_mapping.txt'
                elif (config.LRS2_spec == 'B') and (side == 'R'):
                    IFUfile = mapdir+'LRS2_B_OR_mapping.txt'
                elif (config.LRS2_spec == 'R') and (side == 'L'):
                    IFUfile = mapdir+'LRS2_R_NR_mapping.txt'
                elif (config.LRS2_spec == 'R') and (side == 'R'):
                    IFUfile = mapdir+'LRS2_R_FR_mapping.txt'

                psf      = 1.5
                basename = f[2:-7] + '_' + side 
                outname  = f[0:-5] 

                #call writeDither to build the dither file needed for mkcube
                ditherf = open(ditherfile, 'w')
                ditherinfo.writeHeader(ditherf)
                ditherinfo.writeDither(ditherf,basename,"../mastertrace_"+str(uca)+'_'+side, 0.0, 0.0, psf, 1.00, airmass)

                mkcube(IFUfile,ditherfile,outname,config.diffAtmRef,cubeopts) 

            #cd back into the reduction directory 
            os.chdir('../../../')

    # Run collapse cube
    if config.collapseCube:
        print ('***************************************')
        print ('* COLLAPSING DATA CUBE TO BUILD IMAGE *')
        print ('***************************************')

        numCufiles = glob.glob(redux_dir + "/" + sci_dir + "/*/" + "CuFeR*_sci_*.fits")
        #makes sure there are actually data cubes made from wavelength resampled, fiber extracted data in the sci directory 
        #If data cubes were made from fiber extracted fibers that do not wl resample they do not contain WCS info needed
        if len(numCufiles) == 0:
            sys.exit("You must build data cubes from wavelength resampled, fiber extracted data before running collapse cube")

        for s in config.sci_objects:
            #cd into the science directory 
            location_prefix = redux_dir + "/" + sci_dir + "/" + s + "/" 
            Cufiles = glob.glob(location_prefix + "CuFeR*_sci_*.fits")

            #user defined wavelength range to collapse cube 
            low_wave  = config.col_wave_range[0]
            high_wave = config.col_wave_range[1]

            #track number of cubes used in order to inform user if their values fall out of bounds and no images made
            num_cubes = 0

            min_wave_set = []
            max_wave_set = []
            #iterate through cube files 
            for c in Cufiles:
                im  = pyfits.open(c)
                hdr = im[0].header
                dat = im[0].data

                #read header for wavelength solution information
                lenZ  = np.shape(dat)[0]
                CRVAL = hdr['CRVAL3']                                        
                CDELT = hdr['CDELT3']
                Side  = hdr['CCDPOS']

                #Find the name of the channel for this data cube
                if   (config.LRS2_spec == 'B') and (Side == 'L'):
                    spec_chan = 'UV'
                elif (config.LRS2_spec == 'B') and (Side == 'R'):
                    spec_chan = 'orange'
                elif (config.LRS2_spec == 'R') and (Side == 'L'):
                    spec_chan = 'red'
                elif (config.LRS2_spec == 'R') and (Side == 'R'):
                    spec_chan = 'far-red'

                #build wavelength solution and find min and max wavelength of that solution
                wave_sol = np.add(np.arange(0,(lenZ*CDELT)+1,CDELT),CRVAL)
                max_wave = np.amax(wave_sol)
                min_wave = np.amin(wave_sol)

                #append the values for each frame to inform user of bounds of this data set if their values are out of bounds
                max_wave_set.append(max_wave)
                min_wave_set.append(min_wave)

                #If they choose to collapse entire cube ([0,0]) all data cubes are collapsed into images
                if low_wave == 0 and high_wave == 0:
                    num_cubes = num_cubes + 1 
                    print('Building image from '+spec_chan+' channel cube: '+op.basename(c))
                    print('Collapsing entire cube')

                    sum_image  = np.sum(dat, axis=0) #sums data cube in z direction
                    pyfits.writeto( location_prefix+'Col'+op.basename(c), sum_image, header=hdr, clobber=True)
                    print('\n')

                #If the user choose a wavelength range check if it in range of this cube 
                if (low_wave >= min_wave) and (high_wave <= max_wave):
                    num_cubes = num_cubes + 1 
                    print('Building image from '+spec_chan+' channel cube: '+op.basename(c))
                    print('Collapsing cube from '+str(low_wave)+' to '+str(high_wave))

                    #find what index these wavelengths most closely correspond to. 
                    low_ind  = (np.abs(wave_sol-low_wave)).argmin()
                    high_ind = (np.abs(wave_sol-high_wave)).argmin()

                    #find wavelength value at this index - not exactly users choice so want to print value
                    low_val  = str(wave_sol[low_ind]).split('.')[0]
                    high_val = str(wave_sol[high_ind]).split('.')[0]

                    #Build image from that data cube 
                    dat_region = dat[low_ind:high_ind,:,:] #build region from low to high z - include all x,y
                    sum_image  = np.sum(dat_region, axis=0) #sum the image in the z direction
                    pyfits.writeto( location_prefix+'Col_'+low_val+'_'+high_val+'_'+op.basename(c), sum_image, header=hdr, clobber=True)
                    print('\n')

            #if num_cubes is 0: all cubes out of wavelength range of users choice 
            if len(Cufiles) > 0 and num_cubes == 0:
                print ("Wavelength range you choose for collapse cube is out of range")
                sys.exit("This LRS2-"+config.LRS2_spec+" data set ranges between "+str(np.amin(min_wave_set))+" and "+str(np.amax(max_wave_set))+" Angstroms")

    return vframes
    
def main():
    frames = basicred( DIR_DICT, config.sci_objects, redux_dir, basic = config.basic, dividepf = dividepf,
                      normalize = normalize, masterdark = masterdark, masterarc = masterarc, mastertrace = mastertrace )                 
    
if __name__ == '__main__':
    main()  
