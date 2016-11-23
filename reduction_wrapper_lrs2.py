# -*- coding: utf-8 -*-
"""
Created on Friday April  8 02:10:02 2016

Reduces LRS2 data for either the blue or red channels
This script reads user parameters from lrs2_config.py 

@author: gregz, brianaindahl
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
from lrs2_config import * 

#blu 370-470nm
#org 460-700nm
#red 650-842nm
#frd 818-1050nm

#################################
# Defining which unit to reduce #
#################################

#specifying LRS2 unit to reduce 
if LRS2_spec == 'B':
    redux_dir   = redux_dir+'_B'
    print ('#########################################')
    print ('## RUNNING DATA REDUCTION ON LRS2-BLUE ##')
    print ('#########################################')
elif LRS2_spec == 'R':
    redux_dir   = redux_dir+'_R'
    print ('########################################')
    print ('## RUNNING DATA REDUCTION ON LRS2-RED ##')
    print ('########################################')
else:
    sys.exit('You need to choose either R or B for LRS2_spec')

##################################################
# Setting CUREBIN and check LRS2 defined in CURE #
##################################################

#setting CUREBIN
CUREBIN = None
if not CUREBIN:
    CUREBIN = environ.get('CUREBIN')
if not CUREBIN:
    sys.exit("Please set CUREBIN as environment variable")

#checking that LRS2 is defined in specconf.h 
cureversion = os.popen(op.join(CUREBIN, 'cureversion')).readlines()
spec_define = cureversion[4].split(' ')[1]
instrument = spec_define.rstrip('\n')

if instrument == 'LRS2':
    print ('CURE is set for LRS2 reduction')
else: 
    sys.exit('You need to update specconf.h in CURE to define LRS2')

###################################
# Defining which functions to run #
###################################

#dark folder does not get used if there are not darks
#if there are no darks it sets drk_folder to zro_folder for the script to run.
#nothing is done with these fake dark files
if len(drk_folder) == 0:
    ifdarks = False
    drk_folder = zro_folder
else: 
    ifdarks == True

#if basic reduction is run need to specify specific routines to run 
# divide pixel flat and masterdark are not being used now
if basic:
    rmcosmics       = rmCosmics 
    fix_chan        = True
    dividepf        = dividePixFlt
    normalize       = False
    masterdark      = ifdarks
    subtractdark    = ifdarks
    masterarc       = True  
    mastertrace     = True 
else:
    rmcosmics       = False
    fix_chan        = False
    dividepf        = False
    normalize       = False  
    masterdark      = False
    masterarc       = False  
    mastertrace     = False

# This makes sure that the redux folder is only overwritten if the user chooses to run basic reduction
# If you user only wants to run deformer, skysubtract, fiberextract, or mkcube it used the data in redux 
if basic:
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
subskyopts      = "-J --output-both -w "+str(window_size)+" -k "+str(sky_kappa[0])+','+str(sky_kappa[1])+" -m "+str(smoothing)+" -T "+str(sn_thresh)
fibextractopts  = "-P"
cubeopts        = "-a "+str(sky_sampling)+" -k "+str(max_distance)+" -s "+str(cube_sigma)

#########################
# Defining data folders #
#########################

#reads folder names from config file
zro_file_loc = []
for z in zro_folder:
    zro_file_loc.append(date_folder+"/"+z)

drk_file_loc = []
for d in drk_folder:
    drk_file_loc.append(date_folder+"/"+d)

flt_file_loc = []
for f in flt_folder:
    flt_file_loc.append(date_folder+"/"+f)

sci_file_loc = []
for s in sci_folder:
    sci_file_loc.append(date_folder+"/"+s)

Hg_file_loc = []
for h in Hg_folder:
    Hg_file_loc.append(date_folder+"/"+h)

Cd_file_loc = []
for c in Cd_folder:
    Cd_file_loc.append(date_folder+"/"+c)

FeAr_file_loc = []
for a in FeAr_folder:
    FeAr_file_loc.append(date_folder+"/"+a)

if   LRS2_spec == 'B':
    cmp_file_loc = Hg_file_loc + Cd_file_loc
elif LRS2_spec == 'R':
    cmp_file_loc = Hg_file_loc + FeAr_file_loc  

#specify folders where data types are stored
zro_dir  = "zro"
flt_dir  = "flt"
sci_dir  = "sci"
cmp_dir  = "cmp"
drk_dir  = "drk"

##########################
# Building Spec Libaries #
##########################

#dictionary of data type folders 
file_loc_dir = [ zro_file_loc, drk_file_loc, cmp_file_loc, flt_file_loc, sci_file_loc ] # Should match order of dictionary, otherwise there will be mismatches
DIR_DICT     = {    0:zro_dir,    1:drk_dir,    2:cmp_dir,    3:flt_dir,    4:sci_dir } # Dictionary for looping through directory names

#Detector Amps and spectrograph side lists 
SPEC = ["LL","LU","RL","RU"]
SPECBIG = ["L","R"]

#############################
# Define config directories #
#############################

#specifies directories for lines and mapping/cen files for LRS2 in the LRS2 config directory 
linesdir    = configdir + '/lines_files'
mapdir      = configdir + '/mapping_files/'
pixflatdir  = configdir + '/pixel_flats/'

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
            self.clean                  = CLEAN_AFTER_DONE

            ###### READ IN THE HEADER AND NECESSARY KEYWORDS ######
            self.trimsec    = {}
            self.biassec    = {}
            self.actionbase = {}
            for amp in SPEC:    
                rootname          = op.join (self.origloc, self.basename + '_' + self.ifuslot + amp + '_' + self.type + '.fits' )
                hdulist           = pyfits.open ( rootname )
                #self.trimsec[amp] = "\"" + re.split('[\[ \] ]',hdulist[0].header['TRIMSEC'])[1] + "\""
                #self.biassec[amp] = "\"" + re.split('[\[ \] ]',hdulist[0].header['BIASSEC'])[1] + "\""
                self.trimsec[amp] = "2:2065,1:1032" 
                self.biassec[amp] = "2066:2128,1:1032"      
                self.actionbase[amp] = ''            

            self.actionbase["L"] = initial_base  
            self.actionbase["R"] = initial_base  
            
            self.specid      = str(hdulist[0].header['SPECID']) 
            self.orggain     = hdulist[0].header['GAIN']
            self.orgrdnoise  = hdulist[0].header['RDNOISE']
            self.exptime     = hdulist[0].header['EXPTIME']

            #for cal data: object will have R or B appended for the liquid light guide used 
            #for sci data: object will have R or B appended for the unit that was pointed at the sci target
            #this is not true for older data and is sometimes inconsistent so if no tag cal_side = None
            if len(hdulist[0].header['OBJECT'].split('_')) == 1:
                self.object   = hdulist[0].header['OBJECT']
                self.cal_side = None
            elif self.type != 'sci':
                self.object   = hdulist[0].header['OBJECT'].split('_')[0]
                self.cal_side = hdulist[0].header['OBJECT'].split('_')[1]
            else:
                #sci objects can have multiple '_' so have to check if the last one is a tag R or B
                end_term = hdulist[0].header['OBJECT'].split('_')[-1]
                if (end_term == 'R') or (end_term == 'B'):
                    self.object   = hdulist[0].header['OBJECT'][0:-2]
                    self.cal_side = end_term
                else:
                    self.object   = hdulist[0].header['OBJECT']
                    self.cal_side = None

            #print (self.type, self.object, self.cal_side)

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
        return os.system(op.join(CUREBIN, command) +' 1>>output.log  2>> error.log')
    else:
        return os.system(op.join(CUREBIN, command))
        
        
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

def subtractdark(frames, masterdarkname, amp):
        
    filenames = [op.join ( f.origloc, f.actionbase[amp] + f.basename + '_' + f.ifuslot + amp + '_' + f.type + '.fits') for f in frames]
    exptimes  = [ f.exptime for f in frames ]

    command_2 = 'subtractfits -p '' -f %s' % ( mastername ) 
        
    for i in xrange(len(filenames)):
        command_1 = 'multiplyfits -p d%s_ -c %s' % ( exptimes[i],exptimes[i] ) 
        command_1 = command_1 + ' ' + masterdarkname
        
        run_cure_command( command_1, 0 )

        command_2 = command_2 + ' ' + filenames[i]
        run_cure_command( command_2, 0 )

    return command
    

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
        c = cosmics.cosmicsimage(array, gain=im_gain, readnoise=im_RN, sigclip = 7.0, sigfrac = 0.3, objlim = 7.0)
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

    command = 'addphotonnoise'
    
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
    
    filenames = [(redux_dir + '/' + sci_dir + '/pses' + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in frames]

    for f in filenames:
        
        command = 'subtractsky %s %s -d %s -f %s %s' % (opts,skymaster,distmodel,fibermodel,f)  

        run_cure_command( command, 0 )
        
    [f.addbase('S', side) for f in frames]

    return command
    
    
def fibextract_Resample(frames,base,side,distmodel,fibermodel,wave_range,dw,opts):

    filenames = [(redux_dir + '/' + sci_dir + '/' + base + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in frames]

    for f in filenames:
    
        command = 'fiberextract %s -p FeR -W %s -w %s -d %s -f %s %s' %(opts,wave_range,dw,distmodel,fibermodel,f)
        
        run_cure_command( command, 0 )

    return command

def fibextract(frames,base,side,distmodel,fibermodel,opts):

    filenames = [(redux_dir + '/' + sci_dir + '/' + base + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in frames]

    for f in filenames:
    
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

def initial_setup ( file_loc_dir = None, redux_dir = None, DIR_DICT = None):
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
    vframes = [] # will fill this list with VirusFrame class objects for each image
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

    if file_loc_dir is None:        
        print ( "Please provide a file location directory \"file_loc_dir\"." )
        print ( "This is a list of the location directories for each file type." )
        return None
        
    if DIR_DICT is None:        
        print ( "Please provide a directory dictionary \"DIR_DICT\"." )
        print ( "The directory dictionary order should match the file location directory." )
        return None
    
    # Loop through the file location directories     
    for i in xrange ( len ( file_loc_dir ) ):
        # If the file_loc_dir[i] is None, then nothing is done
        if file_loc_dir[i] is not None:
            # If the reduction location exists, don't re-make the directory (also, any files in that directory remain)
            if not op.exists ( op.join ( redux_dir, DIR_DICT[i] ) ):
                os.mkdir ( op.join ( redux_dir, DIR_DICT[i] ) )
            # Run through the list of ObsID's, exposures, or other structures to gather all desired images for a particular type (e.g., zro)
            for file_loc in file_loc_dir[i]:
                # Trying to figure out what part of the directory structure was given to properly copy it
                dir1 = op.basename ( file_loc )
                if len ( dir1 ) > 5:
                    fit_path  = "/exp*/lrs2/*.fits"
                    filenames = glob.glob ( file_loc + fit_path )
                elif len ( dir1 ) == 5:
                    if dir1[0:3] == "exp": 
                        fit_path  = "/lrs2/*.fits"                   
                        filenames = glob.glob ( file_loc + fit_path )                    
                    else:      
                        fit_path  = "/*.fits"             
                        filenames = glob.glob ( file_loc + fit_path )                   
                else:               
                    print ( "Did not recognize the " + DIR_DICT[i] + " basename as \"lrs2XXXXXXX\", \"expXX\", or \"lrs2\"." )
                    return None

                #this will perform a header check of the data when looping through the flats
                #checks if it is early data and special cals are needed to be copied from the config file
                if (file_loc_dir[i] == flt_file_loc) and (file_loc == file_loc_dir[i][0]):
                    print ('******************************')
                    print ('* Header Check for Data Date *')
                    print ('******************************')
                    #the header from the first flat frame is used 
                    file1 = filenames[0] 
                    a1 = VirusFrame( file1 )
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
                        if (LRS2_spec == "R") and (all_copy):
                            print ("Calibration data is not optimal - some reference flats and comps will be taken from lrs2_config")
                            #if first_run (old data) and reducing Red - copies long Qth exposures from config file to flt directory
                            longQth_files = glob.glob(configdir+'/FR_longCals/long_Qth'+fit_path)
                            print ('Copying long exposure Qth flats for far-red channel reduction')
                            for q in longQth_files:
                                shutil.copy ( q , op.join ( redux_dir, DIR_DICT[i] ) )

                #If reducing Red data: copy long FeAr exposures to cmp directory for far-red reduction                
                if (LRS2_spec == "R") and (all_copy):
                    if (file_loc_dir[i] == cmp_file_loc) and (file_loc == file_loc_dir[i][0]):
                        longFeAr_files = glob.glob(configdir+'/FR_longCals/long_FeAr'+fit_path)
                        print ('Copying long exposure FeAr comps for far-red channel reduction')
                        for fa in longFeAr_files:
                            shutil.copy ( fa , op.join ( redux_dir, DIR_DICT[i] ) )

                # Loop through the retrieved files names to copy to new structure
                # Create a VirusFrame class for each frame that can maintain info for each original frame
                # The VirusFrame is critical for the rest of the reduction pipeline
                # Only gather VirusFrame objects for LL frames as a single VirusFrame serves for all four amps
                if all_copy:
                    #must copy all file to directories first
                    for f in filenames:  
                        shutil.copy ( f, op.join ( redux_dir, DIR_DICT[i] ) )
                    #once files copied it builds the virus frames for each
                    for f in filenames:            
                        temp, temp1, temp2 = op.basename ( f ).split('_')
                        amp                = temp1[3:5]
                        if amp == "LL":
                            a = VirusFrame( op.join( redux_dir, DIR_DICT[i], op.basename ( f ) ) ) 
                            vframes.append(copy.deepcopy(a))
                else:
                    #if not all_copy does not copy files to directories but builds virus frames for each
                    for f in filenames:            
                        temp, temp1, temp2 = op.basename ( f ).split('_')
                        amp                = temp1[3:5]
                        if amp == "LL":
                            a = VirusFrame(  f  ) 
                            vframes.append(copy.deepcopy(a))
                        
    return vframes, first_run

def basicred( file_loc_dir, redux_dir, DIR_DICT, basic = False, dividepf = False,
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
    vframes, first_run = initial_setup ( file_loc_dir, redux_dir, DIR_DICT )

    #first frames to pull basic header information 
    f1 = vframes[0]

    #checks which unit to reduce and sets ucam (specid) - If B decides if old or new specid based of first_run found in intial_setup()
    if LRS2_spec == 'R':
        ucam = "502"
    elif LRS2_spec == 'B':
        #if LRS2-B finds the date of the data taken to know if it is the new or old LRS2-B
        if first_run:
            ucam = "501"
            print ("You are reducing data with the old LRS2-Blue - You must disregard data from the UV channel")
        else:
            ucam = "503"

    #from the vframes makes lists of files vframes according to type and specid (ucam)
    tframes  = [v for v in vframes if (v.specid == ucam) and (v.cal_side == None or v.cal_side == LRS2_spec)] # gives all frames 

    oframes  = [t for t in tframes if t.type != "zro" ] # gives "flt", "drk", "cmp", and "sci" frames (basically just not "zro")
    dframes  = [t for t in tframes if t.type == "drk" ] # gives dark frames
    cframes  = [t for t in tframes if (t.type == "flt" or t.type == "cmp") ] # gives "flt" and "cmp" frames
    lframes  = [t for t in tframes if t.type == "cmp" ] # gives just "cmp" frames
    fframes  = [t for t in tframes if t.type == "flt" ] # gives just "flt" frames
    zframes  = [t for t in tframes if t.type == "zro" ] # gives just "zro" frames
    spframes = [t for t in tframes if t.type == "sci" ] # gives just "sci" frames with correct LRS2_spec pointing

    sframes  = [v for v in vframes if v.type == "sci" and (v.specid == ucam)] # gives just "sci" frames with any pointing

    #Check that data is correct
    if len(lframes) == 0:
        print ("You did not provide the correct comp frames for this data set")
        sys.exit("Either the flt_folder you provided are not flats or these are not for LRS2-"+LRS2_spec)
    if len(fframes) == 0:
        print ("You did not provide the correct flat frames for this data set")
        sys.exit("Either the Hg/Cd/FeAr_folder you provided are not comps or these are not for LRS2-"+LRS2_spec)
    if len(zframes) == 0:
        print ("You did not provide zero frames for this data set")
        sys.exit("Check the data type of the zro_folder you provided to make sure they are zro images")
    if len(sframes) == 0:
        print ("You did not provide science frames for this data set")
        sys.exit("Check the data type of the sci_folder you provided to make sure they are sci images")
    if len(spframes) == 0:
        print ("WARNING: You did not provide sci frames with LRS2-"+LRS2_spec+" pointings - these may just be sky frames")

    #make a copy of the lsr2_config file to be added to your directory
    #if the file exists - remove the file and replace it.
    if os.path.isfile(redux_dir+'/lrs2_config_'+redux_dir+'_copy.py') == True:
        os.remove(redux_dir+'/lrs2_config_'+redux_dir+'_copy.py')
        shutil.copy ( os.path.dirname(os.path.realpath(__file__))+'/lrs2_config.py', redux_dir+'/lrs2_config_'+redux_dir+'_copy.py' )
    else:
        shutil.copy ( os.path.dirname(os.path.realpath(__file__))+'/lrs2_config.py', redux_dir+'/lrs2_config_'+redux_dir+'_copy.py' )

    if basic:
        for sp in SPEC:
            trimsec = f1.trimsec["LL"] # Trimsec assumed to be the same for all frames of a given amp
            biassec = f1.biassec["LL"] # Biassec assumed to be the same for all frames of a given amp
            print ('**************************')
            print ('* MAKE ERROR FRAME FOR '+sp+' *')
            print ('**************************')
            mkerrorframe ( tframes, sp )               # for all frames

            print ('***************************')
            print ('* SUBTRACT OVERSCAN FOR '+sp+' *')
            print ('***************************')
            subtractoverscan( biassec, tframes, sp )   # for all frames

            print ('*****************************')
            print ('* EXTRACT DATA REGION FOR '+sp+' *')
            print ('*****************************')
            extractfits ( trimsec, tframes, sp )       # for all frames

            if rmcosmics:
                print ('******************************************')
                print ('* REMOVE COSMIC RAYS (SCI IMAGES) FOR '+sp+' *')
                print ('******************************************')
                #if reducing old LRS2-Blue data: do not run rmcosmics on UV data(LL+LU), only orange(RL+RU)
                if (LRS2_spec == 'B') and (first_run):
                    if (sp == 'RL') or (sp == 'RU'):
                        rmcosmicfits ( sframes, sp )       # for sci frames - because this is slow
                else:
                    rmcosmicfits ( sframes, sp )       # for sci frames - because this is slow                    
            
            print ('**************************')
            print ('* BUILD MASTERBIAS FOR '+sp+' *')
            print ('**************************')
            vframesselect  = [v for v in tframes if v.type == "zro" and v.specid == ucam] 
            meanbiasfits   ( sp, ucam, redux_dir, 'masterbias', meanfitsopts, vframesselect ) # meanfits for masterbias for unique specid 

            print ('*****************************')
            print ('* SUBTRACT MASTERBIAS FOR '+sp+' *')
            print ('*****************************')
            masterbiasname = op.join ( redux_dir, 'masterbias' + '_' + ucam + '_' + sp + '.fits' ) 
            oframesselect  = [o for o in oframes if o.specid == ucam]
            subtractbias   ( oframesselect, masterbiasname, sp) # for all frames

            # Create Master Dark Frames
            if masterdark:
                print ('********************************************')
                print ('* BUILDING MASTERDARK FOR SCI FRAMES FOR '+sp+' *')
                print ('********************************************')
                dframesselect = [d for d in dframes if d.specid == ucam] 
                meandarkfits ( sp, ucam, redux_dir, 'masterdark', meanfitsopts, dframesselect ) # meanfits for masterdark for unique specid
               
                #Subtracts Master Dark from Science Frames  
                print ('***************************************')
                print ('* SUBTRACTING DARKS FROM SCI FRAMES '+sp+' *')
                print ('***************************************')
                sframesselect = [s for s in sframes if s.specid == ucam] 
                masterdarkname = op.join ( redux_dir, 'masterdark' + '_' + ucam + '_' + sp + '.fits' )
                subtractdark ( sframesselect, masterdarkname, sp) # for sci frames 

                
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
    
    if normalize:
        print ('***************************************************')
        print ('* NORMALIZING CALIBRATION IMAGES BY EXPOSURE TIME *')
        print ('***************************************************')
        # Normalizing Calibration Frames by Exposure Time
        for cal in cframes:
            for side in SPECBIG:
                opt     = '-c {:0.1f}'.format(cal.exptime)
                filename = op.join ( cal.origloc, cal.actionbase[side] + cal.basename + '_' + cal.ifuslot + '_' + cal.type + '_' + side + '.fits' )
                dividefits ( filename, opt )
                cal.addbase ('d', side = side)
                filename = op.join ( cal.origloc, cal.actionbase[side] + cal.basename + '_' + cal.ifuslot + '_' + cal.type + '_' + side + '.fits' )
                headfits ( filename , headfitsopts )
                
    if masterarc:        
        # Making Master Arc Frames
        print ('**********************')
        print ('* BUILDING MASTERARC *')
        print ('**********************')
        if ucam == '501':
            LAMP_DICT = {0:'Hg',1:'Cd'}
        if (ucam == '502') or (ucam == '503'):
            LAMP_DICT = {0:'Hg',1:'FeAr'}
        for side in SPECBIG:
            for lamp in LAMP_DICT.values():
                #Creates a masterarc frame for each arc lamp in LAMP_DICT
                lframesselect = [l for l in lframes if l.specid == ucam and l.object == lamp]
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
                    sys.exit( "For LRS2-"+LRS2_spec+" You must provide "+lamp+" data")

            #Combine each arc master frame for each lamp in LAMP_DICT into one masterarc 
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
         
    # Making Master Trace Frames
    if mastertrace:
        print ('************************')
        print ('* BUILDING MASTERTRACE *')
        print ('************************')

        if LRS2_spec == 'B':
            FLT_LAMP = 'ldls'
        if LRS2_spec == 'R':
            FLT_LAMP = 'Qth'
        for side in SPECBIG:
            #Check that the proper flat lamp images were choosen for the spec being reduced 
            if fframes[0].object.split('_')[0] != FLT_LAMP:
                print ('LRS2-'+LRS2_spec+' requires '+FLT_LAMP+' flats: You choose '+fframes[0].object.split('_')[0]+' flats')
                print ('Please update config file: flt_folder should be folder containing '+FLT_LAMP+' flats')
                return None 
            #Build a list of flat frames matching the correct lamp
            fframesselect = [f for f in fframes if f.specid == ucam and f.object.split('_')[0] == FLT_LAMP] 
            #If there is more than one frame take a median from of the ones provided for the mastertrace
            if len ( fframesselect ) > 1:
                meantracefits(side, ucam, redux_dir, 'mastertrace', traceopts, fframesselect)
            #If there was only one frame provided that frame becomes the mastertrace
            else:
                filename = [op.join ( f.origloc, f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in fframesselect]
                mastername = op.join ( redux_dir , 'mastertrace' + '_' + ucam + '_' + side + '.fits')
                shutil.copy ( filename[0], mastername )
                efilename = [op.join ( f.origloc, 'e.' + f.actionbase[side] + f.basename + '_' + f.ifuslot + '_' + f.type + '_' + side + '.fits') for f in fframesselect]
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
    if fix_chan and (LRS2_spec == 'R'):
        print ('**************************************')
        print ('* FIXING FAR RED CHANNEL MASTERTRACE *')
        print ('**************************************')
        #X value to start with to extend trace to the edge of chip
        start_col = 50

        im = pyfits.open(redux_dir+'/mastertrace_502_R.fits')
        hdr = pyfits.getheader(redux_dir+'/mastertrace_502_R.fits')
        dat =  im[0].data

        dat = extend_trace_start(dat,start_col,10)

        #before writing the new file it moves the old mastertrace to a file with _orig appended to save original file
        shutil.move(redux_dir+'/mastertrace_502_R.fits',redux_dir+'/mastertrace_502_R_orig.fits')
        pyfits.writeto(redux_dir+'/mastertrace_502_R.fits', data = dat, header = hdr, clobber=True)

        #fix far-red channel masterarc

    #Extend traces to detector edges for the orange channel. 
    #This will create a new mastertrace_501_R.fits and save the old one as mastertrace_501_R_orig.fits
    if fix_chan and (LRS2_spec == 'B'):
        print ('*************************************')
        print ('* FIXING ORANGE CHANNEL MASTERTRACE *')
        print ('*************************************')
        #X value to start with to extend trace to the edge of chip
        start_col = 100

        im = pyfits.open(redux_dir+'/mastertrace_501_R.fits')
        hdr = pyfits.getheader(redux_dir+'/mastertrace_501_R.fits')
        dat =  im[0].data

        dat = extend_trace_start(dat,start_col,10)

        #before writing the new file it moves the old mastertrace to a file with _orig appended to save original file
        shutil.move(redux_dir+'/mastertrace_501_R.fits',redux_dir+'/mastertrace_501_R_orig.fits')
        pyfits.writeto(redux_dir+'/mastertrace_501_R.fits', data = dat, header = hdr, clobber=True)
            
    # Run Deformer
    if run_deformer:
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
            if (LRS2_spec == 'B') and (side == 'L'):
                wave_range = '[3700,4700]'
                ref_line = 6
            if (LRS2_spec == 'B') and (side == 'R'):
                wave_range = '[4600,7000]'
                ref_line = 7
            if (LRS2_spec == 'R') and (side == 'L'):
                wave_range = '[6500,8500]'
                ref_line = 3
            if (LRS2_spec == 'R') and (side == 'R'):
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
    if subsky:  
        print ('************************************************')
        print ('* PERFORMING SKY SUBTRACTION ON SCIENCE FRAMES *')
        print ('************************************************')
        #check that deformer has been run 
        dist_files = glob.glob(op.join(redux_dir,'*.dist'))
        if len(dist_files) == 0:
            sys.exit("You must run deformer before you can run sky subtraction")

        for side in SPECBIG:
            distmodel = op.join ( redux_dir, 'mastertrace' + '_' + ucam + '_' + side + '.dist' )
            fibermodel = op.join ( redux_dir, 'mastertrace' + '_' + ucam + '_' + side + '.fmod' )
            sframesselect = [s for s in sframes if s.specid == ucam]
            subtractsky(sframesselect,side,distmodel,fibermodel,subskyopts)

    # Run fiberextract
    if fiberextract:  
        print ('****************************************')
        print ('* EXTRACTING SPECTRA IN SCIENCE FRAMES *')
        print ('****************************************')
        #check that deformer has been run 
        dist_files = glob.glob(op.join(redux_dir,'*.dist'))
        if len(dist_files) == 0:
            sys.exit("You must run deformer before you can run fiber extract")

        #finds if there are sky subtracted files. If so it uses those.
        Sfiles = glob.glob(redux_dir + "/" + sci_dir + "/Sp*_sci_*.fits")
        if len(Sfiles) == 0:
            base = 'pses'
        else:
            base = 'Spses'

        if wl_resample:
            print ('    ++++++++++++++++++++++++++')
            print ('     Resampling in Wavelength ')
            print ('    ++++++++++++++++++++++++++')
            for side in SPECBIG:
                #for each channel selects correct wavlength range and dw
                if (LRS2_spec == 'B') and (side == 'L'):
                    wave_range = '3643,4668'
                    dw = '0.496'
                if (LRS2_spec == 'B') and (side == 'R'):
                    wave_range = '4600,7000'
                    dw = '1.2'
                if (LRS2_spec == 'R') and (side == 'L'):
                    wave_range = '6432,8451'
                    dw = '0.978'
                if (LRS2_spec == 'R') and (side == 'R'):
                    wave_range = '8324,10565'
                    dw = '1.129'

                sframesselect = [s for s in sframes if s.specid == ucam]

                distmodel = redux_dir + "/mastertrace_" + ucam + "_" + side + ".dist"
                fibermodel = redux_dir + "/mastertrace_" + ucam + "_" + side + ".fmod"
                fibextract_Resample(sframesselect,base,side,distmodel,fibermodel,wave_range,dw,fibextractopts) 
        else:
            print ('    +++++++++++++++++++++++++++++++')
            print ('     Extraction Without Resampling ')
            print ('    +++++++++++++++++++++++++++++++')
            for side in SPECBIG:   
                sframesselect = [s for s in sframes if s.specid == ucam]

                distmodel = redux_dir + "/mastertrace_" + ucam + "_" + side + ".dist"
                fibermodel = redux_dir + "/mastertrace_" + ucam + "_" + side + ".fmod"
                fibextract(sframesselect,base,side,distmodel,fibermodel,fibextractopts)

    #CURE saves these files from deformer outside of the redux directory for some reason.
    #This moves them inside of the redux directory.
    left_files = glob.glob('*.log') + glob.glob('*.residuals')
    if len(left_files) > 0:
        for l in left_files:
            os.rename(l, op.join(redux_dir,l))

    #Run mkcube
    if makecube:
        print ('***********************')
        print ('* BUILDING DATA CUBES *')
        print ('***********************')
        #cd inside of the science directory 
        location_prefix = redux_dir + "/" + sci_dir + "/" 
        os.chdir(location_prefix)

        #builds a list of files to build a data cube 
        #checks for if there are wavelength resampled fiber extracted and sky subtracted files
        Fefiles = glob.glob("FeRS*_sci_*.fits")
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
            if   (LRS2_spec == 'B') and (side == 'L'):
                IFUfile = mapdir+'LRS2_B_UV_mapping.txt'
            elif (LRS2_spec == 'B') and (side == 'R'):
                IFUfile = mapdir+'LRS2_B_OR_mapping.txt'
            elif (LRS2_spec == 'R') and (side == 'L'):
                IFUfile = mapdir+'LRS2_R_NR_mapping.txt'
            elif (LRS2_spec == 'R') and (side == 'R'):
                IFUfile = mapdir+'LRS2_R_FR_mapping.txt'

            psf      = 1.5
            basename = f[2:-7] + '_' + side 
            outname  = f[0:-5] 

            #call writeDither to build the dither file needed for mkcube
            ditherf = open(ditherfile, 'w')
            ditherinfo.writeHeader(ditherf)
            ditherinfo.writeDither(ditherf,basename,"../mastertrace_"+str(uca)+'_'+side, 0.0, 0.0, psf, 1.00, airmass)

            mkcube(IFUfile,ditherfile,outname,diffAtmRef,cubeopts) 

        #cd back into the reduction directory 
        os.chdir('../../')

    #Run mkcube
    if collapseCube:
        print ('***************************************')
        print ('* COLLAPSING DATA CUBE TO BUILD IMAGE *')
        print ('***************************************')
        #cd into the science directory 
        location_prefix = redux_dir + "/" + sci_dir + "/" 
        Cufiles = glob.glob(location_prefix + "CuFeR*_sci_*.fits")

        #makes sure there are actually data cubes made from wavelength resampled, fiber extracted data in the sci directory 
        #If data cubes were made from fiber extracted fibers that do not wl resample they do not contain WCS info needed
        if len(Cufiles) == 0:
            sys.exit("You must build data cubes from wavelength resampled, fiber extracted data before running collapse cube")

        #user defined wavelength range to collapse cube 
        low_wave  = col_wave_range[0]
        high_wave = col_wave_range[1]

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
            if   (LRS2_spec == 'B') and (Side == 'L'):
                spec_chan = 'UV'
            elif (LRS2_spec == 'B') and (Side == 'R'):
                spec_chan = 'orange'
            elif (LRS2_spec == 'R') and (Side == 'L'):
                spec_chan = 'red'
            elif (LRS2_spec == 'R') and (Side == 'R'):
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
        if num_cubes == 0:
            print ("Wavelength range you choose for collapse cube is out of range")
            sys.exit("This LRS2-"+LRS2_spec+" data set ranges between "+str(np.amin(min_wave_set))+" and "+str(np.amax(max_wave_set))+" Angstroms")

    return vframes
    
def main():
    frames = basicred( file_loc_dir, redux_dir, DIR_DICT, basic = basic, dividepf = dividepf,
                      normalize = normalize, masterdark = masterdark, masterarc = masterarc, mastertrace = mastertrace )                 
    
if __name__ == '__main__':
    main()  
