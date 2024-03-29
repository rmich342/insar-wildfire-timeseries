#!/usr/bin/env python
"""Module for utility functions

A. Chen, 2011
"""

import os
import sys
import string
import stat
import time

import numpy as np

def run_command(com,verbose=True):
    if(verbose):
       print("\n ***** "+com)
    os.system(com)
    if(verbose):
       print("\n")
    sys.stdout.flush()

def run_matlab(fin):
    """fin is name of Matlab script to run"""

    matlab_path = "/home/grpsw/bin/matlab -nosplash -nodesktop 2> /dev/null < "
    print("\n *** MATLAB "+fin)
    os.system(matlab_path+fin)
    print("\n")
    
def numlines(filename,linelength,type):
    """numlines(filename,linelength,type)
      type can be:
        'ri' for complex pixel-interleaved
        'byte' for 1-byte unsigned int
        'float' for 4-byte IEEE float
        'int' for 4-byte IEEE int
        'int*2' for 2-byte IEEE int
    """
    # get command line arguments
    filename = string.strip(filename)
    linelength = int(linelength)
    type = string.strip(type)

    # get file size in bytes
    filesize = os.stat(filename)[stat.ST_SIZE]

    pixel_bytes = 0
    if type == 'ri':
        pixel_bytes = 8
    if type == 'byte':
        pixel_bytes = 1
    if type == 'float':
        pixel_bytes = 4
    if type == 'int':
        pixel_bytes = 4
    if type == 'int*2':
        pixel_bytes = 2
        
    # check if there are a round number of lines
    if filesize%(linelength*pixel_bytes) != 0:
        print('Error: fractional number of lines!\n')
        sys.exit(1)
        
    # calculate number of lines
    numlines = filesize/(linelength*pixel_bytes)
    return numlines
                                                                            
def image_corners(rgbins): 
   lat = np.fromfile('lat',dtype='float32',count=-1,sep="")\
         .reshape((-1,rgbins),order='C')
   lon = np.fromfile('lon',dtype='float32',count=-1,sep="")\
         .reshape((-1,rgbins),order='C')
   corners = np.array([ [lon[0,0],lat[0,0]],
                        [lon[0,-1],lat[0,-1]],
                        [lon[-1,-1],lat[-1,-1]],
                        [lon[-1,0],lat[-1,0]] ])
   return corners

########### UTILITIES FOR READING FILES ##############

def read_cpx_file(name,rgbins):
   """Reads a pixel-interleaved (PIL) format complex file.
      name is file name 
      rgbins is number of range bins"""
   s = np.fromfile(name,dtype='complex64',count=-1,sep="").\
      reshape((-1,rgbins),order='C')
   return s

def read_cor_file(name,rgbins):
    """Writes the float arrays amp and hgt to file.
       name is name of file.
    """
    s = np.fromfile(name,dtype='float32',count=-1,sep="").\
      reshape((2*int(rgbins)),-1,order='C')
    amp=s[0:(int(rgbins)),:]
    cor=s[int(rgbins):,:]
    return (amp,cor)

def read_hgt_file(name,rgbins):
   """Reads a line-interleaved format float file.
      name is file name 
      rgbins is number of range bins"""
   s = np.fromfile(name,dtype='float32',count=-1,sep="").\
      reshape((-1,2*int(rgbins)),order='C')
   amp = s[:,:int(rgbins)]
   hgt = s[:,int(rgbins):]
   return (amp,hgt)

def read_int_file(name,rgbins):
   """Reads a int16 file.
      name is file name 
      rgbins is number of range bins"""
   s = np.fromfile(name,dtype='int16',count=-1,sep="").\
      reshape((-1,rgbins),order='C')
   return s

def read_byte_file(name,rgbins):
   """Reads a uint8 file.
      name is file name 
      rgbins is number of range bins"""
   s = np.fromfile(name,dtype='uint8',count=-1,sep="").\
      reshape((-1,rgbins),order='C')
   return s

def read_float_file(name,rgbins):
   """Reads a float32 file.
      name is file name 
      rgbins is number of range bins"""
   s = np.fromfile(name,dtype='float32',count=-1,sep="").\
      reshape((-1,rgbins),order='C')
   return s

def read_params2dict(filename):
    """For a file called "filename", this function reads its lines
    one at a time. For any line that is not blank, it reads the first
    word as a parameter name, and the second word as a parameter value.
    This function returns a dictionary of name:value pairs.
    """
    proj_params = {}

    fid = open(filename, 'r')
    for line in fid:
       L = line.split()
       if( len(L) > 0 ):
          proj_params[ L[0] ] = L[1]
    fid.close()

    return proj_params

########### UTILITIES FOR WRITING FILES ##############

def write_hgt_file(name,amp,hgt):
    """Writes the float arrays amp and hgt to file.
       name is name of file.
    """
    s = np.zeros((2*amp.shape[0],amp.shape[1]),dtype='float32',order='C')
    s[0::2,:] = amp
    s[1::2,:] = hgt
    s.tofile(name)


def write_float_file(name,amp):
    """Writes the float arrays amp and hgt to file.
       name is name of file.
    """
    s = np.zeros((amp.shape[0],amp.shape[1]),dtype='float32',order='C')
    s = amp
    s.tofile(name)

def write_int_file(name,amp):
    """Writes the float arrays amp and hgt to file.
       name is name of file.
    """
    s = np.zeros((amp.shape[0],amp.shape[1]),dtype='int16',order='C')
    s = amp
    s.tofile(name)

def write_cor_file(name,amp,length):
    """Writes the float arrays amp and hgt to file.
       name is name of file.
    """
    s = np.zeros((2*amp.shape[0],amp.shape[1]),dtype='float32',order='C')
    s[0:(int(length)),:] = amp
    s[int(length):,:] = amp
    s.tofile(name)

def write_uavsar_cor_file(name,amp,length):
    """Writes the float arrays amp and hgt to file.
       name is name of file.
    """
    s = np.zeros((amp.shape[0],2*amp.shape[1]),dtype='float32',order='C')
    s[:,0:(length)] = amp
    s[:,length:] = amp
    s.tofile(name)

def write_cpx_file(name,amp,length):
    """Writes the cpx file out in PIL format.
       name is name of file.
    """
    s = np.zeros((amp.shape[0],amp.shape[1]),dtype='complex64',order='C')
    s=amp
    s.tofile(name)


########## GMT UTILITIES ##########

gmt_path = '/home/acchen/stdsw/GMT/bin/'
Rgreenland = ' -R-55.798/59.199/7.548/80.072r '

def GMT_settings(width=7,length=10):
    """Initializes GMT using the GMT command gmtset
    usage: GMT_settings(width=7,length=10)
       width is width of paper in inches
       length is length of paper in inches
    """
    run_command( \
       gmt_path+ \
       'gmtset PAPER_MEDIA Custom_'+str(width)+'ix'+str(length)+'i '+\
       'VERBOSE true ' + \
       'OBLIQUE_ANNOTATION 6 ' + \
       'BASEMAP_TYPE plain ' + \
       'HEADER_FONT_SIZE 18 ' + \
       'LABEL_FONT_SIZE 12 ' + \
       'BASEMAP_AXES SWne ' + \
       'PLOT_DEGREE_FORMAT ddd:mm:ss ')
    print(gmt_path)

# ############ ProgressBar class ###########
# class ProgressBar:
#     """ Creates a text-based progress bar. Call the object with 
#     the simple `print'command to see the progress bar, which looks 
#     something like this:

#     [=======> 22% ]

#     You may specify the progress bar's width, min and max values on init.
    
#     .. note::
#         Taken from GIAnT tsio.py by P.Agram, with
#         code originally from http://code.activestate.com/recipes/168639/"""

#     def __init__(self, minValue = 0, maxValue = 100, totalWidth=None):
#         self.progBar = "[]" # This holds the progress bar string
#         self.min = minValue
#         self.max = maxValue
# 	self.span = maxValue - minValue

# 	if totalWidth is None:
# 	    self.width = 50
# 	else:
# 	    self.width = totalWidth

#         self.reset()

#     def reset(self):
#         '''Reset the counter to zero.'''
#         self.start_time = time.time()
#         self.amount = 0 # When amount == max, we are 100% done
#         self.updateAmount(0) # Build progress bar string

#     def updateAmount(self, newAmount = 0):
#         """ Update the progress bar with the new amount (with min and max
#         values set at initialization; if it is over or under, it takes the
#         min or max value as a default. """
#         if newAmount < self.min:
#             newAmount = self.min
#         if newAmount > self.max:
#             newAmount = self.max
#         self.amount = newAmount

#         # Figure out the new percent done, round to an integer
#         diffFromMin = np.float(self.amount - self.min)
#         percentDone = (diffFromMin / np.float(self.span)) * 100.0
#         percentDone = np.int(np.round(percentDone))

#         # Figure out how many hash bars the percentage should be
#         allFull = self.width - 2 - 8
#         numHashes = (percentDone / 100.0) * allFull
#         numHashes = np.int(np.round(numHashes))

#         # Build a progress bar with an arrow of equal signs; special cases for
#         # empty and full
#         if numHashes == 0:
#             self.progBar = '[>%s]' % (' '*(allFull-1))
#         elif numHashes == allFull:
#             self.progBar = '[%s]' % ('='*allFull)
#         else:
#             self.progBar = '[%s>%s]' % ('='*(numHashes-1),
#                             ' '*(allFull-numHashes))
#             # figure out where to put the percentage, roughly centered
#             percentPlace = (len(self.progBar) / 2) - len(str(percentDone))
#             percentString = ' ' + str(percentDone) + '% '
#             elapsed_time = time.time() - self.start_time
#             # slice the percentage into the bar
#             self.progBar = ''.join([self.progBar[0:percentPlace], percentString,
#                     self.progBar[percentPlace+len(percentString):], ])
#             if percentDone > 0:
#                 self.progBar += ' %6ds' % (int(elapsed_time))

#     def update(self, value, every=1):
#         """ Updates the amount, and writes to stdout. Prints a
#          carriage return first, so it will overwrite the current
#           line in stdout."""
#         if value % every == 0 or value >= self.max:
#             self.updateAmount(newAmount=value)
#             sys.stdout.write('\r' + self.progBar)
#             sys.stdout.flush()

#     def close(self):
#         """Prints a blank space at the end to ensure proper printing
#         of future statements."""
#         print ' '

# ##############################  End of progress bar class  ##################################

