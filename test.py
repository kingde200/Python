#!/usr/bin Python

# import all necessary modules

import struct
import math
from   array import *
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches 
from   matplotlib.widgets import Button, Slider
import ConfigParser
import string
import time
import calendar
#import scipy.signal
import numpy as np
#import scipy.io.matlab.byteordercodes
import pylab

import wave
import glob
import os
import subprocess
import gydatos as gy

epsilon = 0.0001



# Read configuration parameters

config = ConfigParser.RawConfigParser()
config.read('test.cfg')

verbose = config.getint('control_parameters', 'verbose')
time_format = config.get('time_interval','format') 
start_time = config.get('time_interval','start_time')
end_time = config.get('time_interval','end_time')
wfdisc = config.get('waveforms', 'wfdisc')
triad = config.get('waveforms','Triad')

station1 = triad+'1'

# Parse time variables and get start and end epoch time

stime = time.strptime(start_time, time_format)
etime = time.strptime(end_time, time_format)
start_e = calendar.timegm(stime)
end_e   = calendar.timegm(etime)

if verbose > 1:
	print "Start time:", stime,'\n'
	print "Epoch start:", start_e
	print "End time:", etime,'\n'
	print "Epoch end:", end_e

wf1=array('f',[])

samprate, null, nsamp = gy.ReadWfm(wf1,wfdisc, start_e, end_e, station1, verbose)

t=array('f',[])
for j in range (len(wf1)):
	t.append(float(j)/samprate)
ax1=plt.subplot(111)
ax1.set_xlim([0.,0.5])
plt.xlim(70.,130.)
w1, = plt.plot(t, wf1)
fm = string.Formatter()
label = fm.format('Waveform {0}',station1)
plt.title(label)
plt.show()
