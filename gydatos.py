#!/usr/bin/env Python

# This module contains a number of general-use utilities 
# to read CSS data and do limited, and extendable, signal 
# processing.
#    
# 
# Developped February 2013 - Ronan Le Bras - Gydatos LLC
#

# import all necessary modules

import struct
import math
from   array import *
import matplotlib.pyplot as plt
from   matplotlib.widgets import Button, Slider
import string
import time
import calendar
import scipy as sc
import scipy.sparse
import numpy as np
import scipy.io.matlab.byteordercodes
import pylab
#import pyaudio
import wave
import glob
import os
import collections
import subprocess
#import matplotlib.mpl as mp
import matplotlib.cm as cm
import cx_Oracle

fm = string.Formatter()

Site = collections.namedtuple('Site', ['sta','lat','lon'], verbose=True)
Detection = collections.namedtuple('Detection',
['sta','Det_num','time','dt_cc','duration','freq','whale_disc','pow_15_25','pow_30_50','pow_50_60','spectim'], 
verbose=True)
wfdisc_tuple = collections.namedtuple('wfdisc',
['sta','chan', 'time', 'wfid', 'chanid', 'jdate', 'endtime', 'nsamp', 'samprate', 'calib', 'calper', 'dir',  'dfile', 'foff'],verbose=True)
KMperDEG = 111.32
LARGE = 9999999.
epsilon = 0.001


def CheckFreqRatio(signal,samprate,frange,ratio):
	fft1 = abs(scipy.fft(signal))
	freqs = scipy.fftpack.fftfreq(len(signal),1./samprate)
	print frange
	band1 = 0.
	band2 = 0.
	band3 = 0.
	nf1 = 0
	nf2 = 0
	nf3 = 0
	for f in range(len(freqs)):
				
		if (freqs[f] > frange[0] and freqs[f] < frange[1]) :
			nf1 += 1
			band1 += fft1[f]
		if (freqs[f] > frange[2] and freqs[f] < frange[3]) :
			nf2 += 1
			band2 += fft1[f]
		if (freqs[f] > frange[4] and freqs[f] < frange[5]) :
			nf3 += 1
			band3 += fft1[f]
	band1 /= float(nf1)
	band2 /= float(nf2)
	band3 /= float(nf3)
	whale = False
	
	if (band2/band1 > ratio):
		 whale=True
	return whale, band1, band2, band3


def HydroGroupSTRefine(hag,w1,w2,w3,sfile,triad,win_front,win_back,samprate,low, high, plot_stack, plot_spectim, percent, test_case):
	'''
	Refine the delta-t times obtained from LTA/STA using crosscorrelation 
        and get an updated direction from these.
	'''
	
	lhag = len(hag)
	win = win_back+win_front
	if(plot_spectim):

		print("Plot Spectim If", int(low), int(high))	
		width = int(high-low)		
		dim = (3*width, lhag*win)
		spim1 = np.zeros(dim)			
				
	
	for i in range(lhag):		
#
# For each group, get the time window, envelope, and cross-correlation
#
		t_start = hag[i].Det1.time-win_front
		t_end = hag[i].Det1.time+win_back
		wf1 = w1[int(t_start*samprate):int(t_end*samprate)]
		wf2 = w2[int(t_start*samprate):int(t_end*samprate)]		
		wf3 = w3[int(t_start*samprate):int(t_end*samprate)]

	
		cmap=mp.colors.ListedColormap([[0.,0.,0.],[0.,0.,0.05],[0.,0.05,0.05],[0.05,0.05,0.05],[0.05,0.05,.1],[.05,0.1,.1],[.1,0.1,.1],[.1,0.1,.2],[.1,0.2,0.2],[0.2,0.2,0.2],[0.2,0.2,0.3],[.2,0.3,.3],[.3,0.3,.3],[.3,0.3,.4],[.3,0.4,.4],[0.4,0.4,0.4],[1.,0.,0.]])
		Pxx1,freqs,bins,im1 = pylab.specgram(
		w1[int((hag[i].Det1.time-win_front)*samprate):int((hag[i].Det1.time+win_back)*samprate)],NFFT=512,Fs=samprate,noverlap=511,cmap=cmap)

		norm=mp.colors.Normalize(vmin=0.,vmax=25.)

		plt.show()
		Pxx2,freqs,bins,im2 = pylab.specgram(
		w2[int((hag[i].Det1.time-win_front)*samprate):int((hag[i].Det1.time+win_back)*samprate)],NFFT=512,Fs=samprate,noverlap=511,cmap=cmap)
		plt.show()
		Pxx3,freqs,bins,im3 = pylab.specgram(
		w3[int((hag[i].Det1.time-win_front)*samprate):int((hag[i].Det1.time+win_back)*samprate)],NFFT=512,Fs=samprate,noverlap=511,cmap=cmap)
		plt.show()
#
# Compute optimum direction based on maximum of stacked power
#		
		
		OptDirection, OptVelocity, theta, radii, pp1, pp2, pp3 = StackOptDirection(Pxx1,Pxx2,Pxx3,bins,freqs,sfile,triad,percent)
		if(plot_spectim):
												
			for m in range((len(pp1)/3)):
				spim1[int(pp1[3*m])-int(low),math.floor(pp1[3*m+1])+i*win] = 1.

			for m in range((len(pp2)/3)):
				spim1[width+int(pp2[3*m])-int(low),math.floor(pp2[3*m+1])+i*win] = 1. 

			for m in range((len(pp3)/3)):
				spim1[2*width+int(pp3[3*m])-int(low),math.floor(pp3[3*m+1])+i*win] = 1.
									
		
		hag[i]=hag[i]._replace(direction_stk = OptDirection)
		optfit = max(radii)
		hag[i]=hag[i]._replace(stkfit = optfit)
		print("Direction and fit from spectrogram stack:",OptDirection,optfit)
		if(plot_stack):
			fig = plt.figure(figsize=(6,6))
			ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
			N=len(radii)
			bwidth = 2.*np.pi/N
			bars = ax.bar(theta, radii, width=bwidth, bottom=0.)
			for r,bar in zip(radii, bars):
				bar.set_facecolor( cm.jet(r/10.))
				bar.set_alpha(0.6)
		
			label = fm.format('HAG nb {0} - Beam Stacking', i)
			label = triad+' '+label

			plt.title(label)
			plt.show()

	if(plot_spectim):
		cmap=mp.colors.ListedColormap([[0.,0.,0.],[1.,0.,0.],[1.,0.,0.],[1.,0.,0.],[1.,0.,0.],[1.,0.,0.],[1.,0.,0.],[1.,0.,0.]])
		plt.imshow(spim1,cmap=cmap)
		plt.show()
		
def WhaleDiscriminant(hag,w1,w2,w3,sfile,win_front,win_back,samprate):
	'''
	Compute the discriminant and update the hag structures
	'''	
	
	

	for i in range(len(hag)):		
#
# For each group, get the time window, envelope, and cross-correlation
#
		centime = hag[i].Det1.time
		pretime = centime-win_front
		postime = centime+win_back
		print('Times of window:', pretime, centime, postime, samprate)
		if(pretime < 0.):
			pretime=0.
		if(int(postime*samprate)  > len(w1)):
			postime = float(len(w1)*samprate)  
		wf1 = w1[int(pretime*samprate):int(postime*samprate)]
		wf2 = w2[int(pretime*samprate):int(postime*samprate)]		
		wf3 = w3[int(pretime*samprate):int(postime*samprate)]

		frange = 15,25,30,50,50,60
		ratio = 1.05
		whale = False
		whale1, pow11, pow12, pow13 = CheckFreqRatio(wf1, samprate, frange, ratio)
		if (whale1):
			whale = True
		whale2, pow21, pow22, pow23 = CheckFreqRatio(wf2, samprate, frange, ratio)
		if (whale2):
			whale = True
		whale3, pow31, pow32, pow33 = CheckFreqRatio(wf3, samprate, frange, ratio)
		if (whale3):
			whale = True
		print ('CCRefine - Whale discriminant',whale1, whale2, whale3, whale)

		Detnew = hag[i].Det1._replace(whale_disc=whale1, pow_15_25=pow11, pow_30_50=pow12, pow_50_60=pow13)
		hag[i] = hag[i]._replace(Det1=Detnew)

		Detnew = hag[i].Det2._replace(whale_disc=whale2, pow_15_25=pow21, pow_30_50=pow22, pow_50_60=pow23)
		hag[i] = hag[i]._replace(Det2=Detnew)		

		Detnew = hag[i].Det3._replace(whale_disc=whale3, pow_15_25=pow31, pow_30_50=pow32, pow_50_60=pow33)
		hag[i] = hag[i]._replace(Det3=Detnew)
		whale = False
		if((math.log(pow12)+math.log(pow22)+math.log(pow32))/(math.log(pow11)+math.log(pow21)+math.log(pow31)) > ratio) :
			whale = True
		hag[i] = hag[i]._replace(whale_disc=whale)		


def HydroGroupCCRefine(hag,w1,w2,w3,sfile,triad,win_front,win_back,samprate,tcor,envelope,refine_crossco,halfwidth,plot_waveform,plot_crossco,plot_detect,plot_spect):
	'''
	Refine the delta-t times obtained from LTA/STA using crosscorrelation 
        and get an updated direction from these.
	'''	
#
# Compute either the crosscorrelation of the signal or the crosscorrelation of the envelopes of the signal in a window set to be win_front seconds before the STA-LTA pick and win_back seconds after.
#
	
	
	
	
	label = fm.format('Weighted crossco. envelope - {0} detections', len(hag))
	label = triad+' '+label
	
	

	for i in range(len(hag)):		
#
# For each group, get the time window, envelope, and cross-correlation
#
		centime = hag[i].Det1.time
		pretime = centime-win_front
		postime = centime+win_back
		print('Times of window:', pretime, centime, postime, samprate)
		wf1 = w1[int(pretime*samprate):int(postime*samprate)]
		wf2 = w2[int(pretime*samprate):int(postime*samprate)]		
		wf3 = w3[int(pretime*samprate):int(postime*samprate)]

		if(plot_waveform):
			time=[]
			for j in range(len(wf1)):
				time.append(float(j-len(wf1)/2)/samprate)
			figc=plt.figure(figsize=(12,6))
			plt.title(label)
			subc1=figc.add_subplot(3,1,1)
			subc2=figc.add_subplot(3,1,2, sharex=subc1)
			subc3=figc.add_subplot(3,1,3, sharex=subc1)
			subc1.plot(time,wf1,'red')
			if (i == 1) :
				label=fm.format('{0}',sta1)
				subc1.text(10.,0.7*ref1,label)
				label=fm.format('{0}',sta2)
				subc2.text(10.,0.7*ref2,label)
				label=fm.format('{0}',sta3)
				subc3.text(10.,0.7*ref3,label)
			subc2.plot(time,wf2,'green')
			subc3.plot(time,wf3,'blue')			
			plt.show()

		

		if(plot_spect):
			figfft = plt.figure(figsize=(10,6))
#			ax = figfft.add_axes([0.1, 0.1, 0.8, 0.8])
			subf1=figfft.add_subplot(1,1,1)						
	                fft1 = abs(scipy.fft(wf1))
			freqs = scipy.fftpack.fftfreq(len(wf1),1./samprate)
			band2050 = 0.
			band5080 = 0.
			nf2050 = 0
			nf5080 = 0
			for f in range(len(freqs)):
				
				if (freqs[f] > 20. and freqs[f] < 50.) :
					nf2050 += 1
					band2050 += fft1[f]
				if (freqs[f] > 50. and freqs[f] < 80.) :
					nf5080 += 1
					band5080 += fft1[f]
			band2050 /= float(nf2050)
			band5080 /= float(nf5080)
			label = fm.format('ratio 5080/2050 is {0}',band5080/band2050)
			plt.title(label)
			plt.xlim(18.,22.)
			subf1.plot(freqs,20.*scipy.log10(fft1),'red')
			


		
#
# Compute Cross-Correlations
#
		lcor = int(tcor*samprate)
		
		crossco0 = array('f',[])
		crossco1 = array('f',[])
		crossco2 = array('f',[])
		crossco3 = array('f',[])
		crossco0 = np.correlate(wf1,wf1,'full')
		crossco1 = np.correlate(wf2,wf1,'full')
#
#
#
		print ("Crosscorrelation 1")
		crossco2 = np.correlate(wf3,wf1,'full')
		crossco3 = np.correlate(wf3,wf2,'full')
		t0 = TimeOfMax(crossco0)
		t1 = TimeOfMax(crossco1)
		t2 = TimeOfMax(crossco2)
		t3 = TimeOfMax(crossco3)
		dt1 = (t1-t0)/samprate
		dt2 = (t2-t0)/samprate
		dt3 = (t3-t0)/samprate
#
# Hilbert transform
#
		if(envelope):
			hwf1 = scipy.signal.hilbert(crossco1)
			print ("Hilbert 1")
			hwf2 = scipy.signal.hilbert(crossco2)
			print ("Hilbert 2")
			hwf3 = scipy.signal.hilbert(crossco3)
			print ("Hilbert 3")
#
# Envelope
#
			wf1 = np.absolute(hwf1)
			wf2 = np.absolute(hwf2)
			wf3 = np.absolute(hwf3)
#
# Compute spectrogram
#
		t1 = TimeOfMax(wf1)
		t2 = TimeOfMax(wf2)
		t3 = TimeOfMax(wf3)
		dt1 = (t1-t0)/samprate
		dt2 = (t2-t0)/samprate
		dt3 = (t3-t0)/samprate
		maxvalue = wf1.max()
		sumvalue = sum(wf1)/samprate
		ref1=maxvalue/sumvalue
		print ("dts",dt1,dt2,dt3)
		for j in range (len(wf1)):
			wf1[j] = wf1[j]/sumvalue
		n_it=0
		dtdiff = len(wf1)/samprate
		while (abs(dtdiff) > .5/samprate and n_it < 100):
			t1,std1 = FitLaplace(wf1, samprate, halfwidth, dt1)
			dt1prev=dt1
			dt1 = t1 -dt1 
			dtdiff = dt1-dt1prev
			n_it += 1
		maxvalue = wf2.max()
		sumvalue = sum(wf2)/samprate
		ref2=maxvalue/sumvalue
		for j in range (len(wf2)):
			wf2[j] = wf2[j]/sumvalue	
		n_it=0
		dtdiff = len(wf1)/samprate
		while (abs(dtdiff) > .5/samprate and n_it < 100):
			t2,std2 = FitLaplace(wf2, samprate, halfwidth, dt2)
			dt2prev=dt2
			dt2 = t2 -dt2 
			dtdiff = dt2-dt2prev
			n_it += 1
		maxvalue = wf3.max()
		sumvalue = sum(wf3)/samprate
		ref3=maxvalue/sumvalue
		for j in range (len(wf3)):
			wf3[j] = wf3[j]/sumvalue
		n_it=0
		dtdiff = len(wf1)/samprate
		while (abs(dtdiff) > .5/samprate and n_it < 100):
			t3,std3 = FitLaplace(wf3, samprate, halfwidth, dt3)
			dt3prev=dt3
			dt3 = t3 -dt3 
			dtdiff = dt3-dt3prev
			n_it += 1
		print ("HAG number, HAG:", i, hag[i])
		
		Detnew = hag[i].Det1._replace(dt_cc=hag[i].Det1.time)
		hag[i] = hag[i]._replace(Det1=Detnew)
		Detnew = hag[i].Det2._replace(dt_cc=hag[i].Det2.time+dt1)
		hag[i] = hag[i]._replace(Det2=Detnew)		
		Detnew = hag[i].Det3._replace(dt_cc=hag[i].Det3.time+dt2)
		hag[i] = hag[i]._replace(Det3=Detnew)
		sta1 = hag[i].Det1.sta
		sta2 = hag[i].Det2.sta
		sta3 = hag[i].Det3.sta
		dt = dt1, dt2, dt3, std1, std2, std3
						
		OptDirection, OptVelocity, theta, radii, maxdiff = CompOptDirection(dt,sfile,triad)
		print ("New direction and fit quality:", OptDirection, maxdiff)
	#		print ("Values:", radii)
		hag[i]=hag[i]._replace(direction_cc=OptDirection)
		hag[i]=hag[i]._replace(relfit=maxdiff)
		print ("HAG number, New HAG:", i, hag[i])
		if(plot_crossco):
				time=[]
				for j in range(len(wf1)):
					time.append(-t0/samprate+float(j)/samprate)
				figc=plt.figure(figsize=(12,6))
				plt.title(label)
				subc1=figc.add_subplot(3,1,1)
				subc2=figc.add_subplot(3,1,2, sharex=subc1)
				subc3=figc.add_subplot(3,1,3, sharex=subc1)
				subc1.plot(time,wf1,'red',dt1,ref1,'bo')
				if (i == 1) :
					label=fm.format('{0}-{1}',sta1,sta2)
					subc1.text(10.,0.7*ref1,label)
					label=fm.format('{0}-{1}',sta1,sta3)
					subc2.text(10.,0.7*ref2,label)
					label=fm.format('{0}-{1}',sta2,sta3)
					subc3.text(10.,0.7*ref3,label)
				subc2.plot(time,wf2,'green',dt2,ref2,'ro')
				subc3.plot(time,wf3,'blue' ,dt3,ref3,'ro')	
				plt.show()		

		if(plot_detect):
			fig = plt.figure(figsize=(6,6))
			ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
			N=len(radii)
			width = 2.*np.pi/N
			bars = ax.bar(theta, radii, width=width, bottom=0.)
			for r,bar in zip(radii, bars):
				bar.set_facecolor( cm.jet(r/10.))
				bar.set_alpha(0.6)
		
			label = fm.format('HAG nb {0} - Refined', i)
			label = triad+' '+label

			plt.title(label)


def FitLaplace(dist, samprate, halfwidth, dt0):
	'''
	Get the best two parameters (average and standard deviation that fits the curve
	First search every coarse time step (10 times step), then refine and interpolate to sub-sampling
	Search for standard deviation steps around the computed standard deviation for each time step
	'''
	
	length = float((len(dist)-1))/samprate
#
# Compute average [sum(x*p(x)*dx)] and stdv [sum(x^2*p(x)*dx] where p(x) is dist (use it as a distribution)
#	
	start = int((-1.*halfwidth+dt0+length/2.)*samprate)
	end = int((halfwidth+dt0+length/2.)*samprate)
	sumd = sum(dist[start:end])/samprate
	avg = 0.
	stdv = halfwidth
	if(sumd != 0.):
	
		for i in range(int(2.*halfwidth*samprate)+1):
			x = -1.*halfwidth+dt0+float(i)/samprate
			j = int((x+length/2.)*samprate)
			avg += x*dist[j]/sumd/samprate
		for i in range(int(2.*halfwidth*samprate)+1):
			x = -1.*halfwidth+dt0+float(i)/samprate
			j=int((x+length/2.)*samprate)
			stdv += (x-avg)*(x-avg)*dist[j]/samprate/sumd

	stdv = math.sqrt(stdv)

#	print("Sum, Average and standard deviation:",sumd,avg+dt0,stdv)
        
	moment = 0.
			
	return avg+dt0, stdv


def ReadSite(sfile, triad, verbose):
	'''
	ReadSite reads a file containing CSS3.0 format site tuples for a hydroacoustic triad.
	 
	Input:
		sfile: name of the ASCII file containing the tuples
		triad: name of the triad being processed
		verbose: verbosity level, the higher the more verbose (max=4) 
	Output:
		Returns two directional vectors in the cordinates system East-North
	'''

	s = open(sfile,"r")
	while True:
		line = s.readline()
		if not line:
			break
		pass
#
# Get rid of all double white space in the string
#
		while '  ' in line:
			line=line.replace("  "," ")
	
		if verbose > 1:
			print (line)
		values = line.split(" ")
#
# Print the whole tuple if the verbosity is higher than 2
#
		if verbose > 2:
	    		print ("Sta=",values[0],'\n', "jdate_on=",values[1], "jdate=",values[2],"\n","Latitude=",values[3],"\n","Longitude=",values[4],"\n","elevation=",values[5],"\n")
	
#
# For now build the name of stations from the name of triad. 
# This has some built-in assumptions.
# Should be smarter in the future.
#
		if (    values[0] != triad+'1' and 
			values[0] != triad+'2' and 
			values[0] != triad+'3'):
			continue

		coord = [values[0],values[3],values[4]]
		if (values[0] == triad+'1'): 
			s1 = Site._make(coord)	 			
		if (values[0] == triad+'2'): 
			s2 = Site._make(coord)	
		if (values[0] == triad+'3'): 
			s3 = Site._make(coord)	

	print (Site)
	y1=KMperDEG*(float(s2.lat)-float(s1.lat))
	x1=KMperDEG*(float(s2.lon)-float(s1.lon))
	y2=KMperDEG*(float(s3.lat)-float(s1.lat))
	x2=KMperDEG*(float(s3.lon)-float(s1.lon))
	
	direction1 = [x1,y1]
	direction2 = [x2,y2]
	s.close()
	return direction1, direction2

def GetTopAmp(p, bins, freqs, percent):
	pp = array('f',[])
	n = 0
#	maxarr=math.log(p.max())
	maxarr=p.max()
	frac = maxarr*percent/100.
	print('Maximum of spectrum:',maxarr)
	for i in range(len(freqs)):
#		print("i=",i)
		for j in range(len(bins)):	
#			print('Value:',math.log(p[i,j]))		
			if(p[i,j] > frac):
				n=n+1
				pp.append(freqs[i])
				pp.append(bins[j])
				pp.append((p[i,j]-frac)/(maxarr-frac))
	print('Number of elements above',percent,'% :',n)
	return pp


def StackOptDirection(p1,p2,p3,bins,freqs,sfile,triad,percent):
	'''
	Computes the best fit direction given three power spectra in time and frequency
	'''
	d1,d2 = ReadSite(sfile, triad, 0)
	'''
	Get the maximum possible travel time for this triad (normalization purpose)
	
	'''	
	radii = array('f',[])
	ang   = array('f',[])
	maxpow = 0.
	minpow = 999999999999999999.
	OptDirection = 0.
	OptVelocity  = 0.
	pp1  = GetTopAmp(p1, bins, freqs, percent)
	pp2  = GetTopAmp(p2, bins, freqs, percent)
	pp3  = GetTopAmp(p3, bins, freqs, percent)
	max1 = p1.max()*len(pp1)/3.
	max2 = p2.max()*len(pp2)/3.
	max3 = p3.max()*len(pp3)/3.
	max0 = math.sqrt(max1*max1+max2*max2+max3*max3)
	samprate = 1./(bins[1]-bins[0])
	for WaterVelocity in ([1.48]):
		power=0.
		for i in range (360):
			theta = 2.*np.pi*float(i)/360.		
			dt1,dt2 = CompDt(d1,d2,theta,WaterVelocity)
			
			
#			print("theta, dt1, dt2, samprate, power=", theta, dt1, dt2, samprate, power)
			power = 0.
			power_max= 0.
			ang.append(theta)
			for pt in range(len(pp1)/3):
				k1=int(pp1[3*pt+1]*samprate)
				k2=k1+int(dt1*samprate)
				k3=k1+int(dt2*samprate)

				power+=pp1[3*pt+2]
				
				for pt2 in range(len(pp2)/3):
					if(int(pp2[3*pt2+1]*samprate) == k2 and pp2[3*pt2] == pp1[3*pt]):
						power+=pp2[3*pt2+2]
#						print("k1, k2, k3=", k1,k2,k3, pp2[3*pt2])
				for pt3 in range(len(pp3)/3):
					if(int(pp3[3*pt3+1]*samprate) == k3 and  pp3[3*pt3] == pp1[3*pt]):
						power+=pp3[3*pt3+2]
#						print("k1, k2, k3=", k1,k2,k3, pp3[3*pt3])
				if(power > power_max):
					power_max = power
			if (power <= minpow):
				minpow = power	
#				print "Direction min, angle (deg), angle (rad), velocity, max, min ", i, theta, WaterVelocity, maxpow, minpow
			if (power > maxpow):
				maxpow = power
				OptDirection = theta
				OptVelocity = WaterVelocity
#				print "Direction max, angle (deg), angle (rad), velocity, max, min ", i, theta, WaterVelocity, maxpow, minpow
			radii.append(power)
			if(maxpow == minpow):
				minpow -= 1.
#		print(radii)
		for r in range(len(radii)):				
			radii[r]=1.-math.erf(1.-(radii[r]-minpow)/((maxpow-minpow)+1.))
	return OptDirection, OptVelocity, ang, radii, pp1,pp2,pp3	

def CompOptDirection(dt,sfile,triad):
	'''
	Computes the best fit direction given three time differences
	'''
	d1,d2 = ReadSite(sfile, triad, 0)
	'''
	Get the maximum possible travel time for this triad (normalization purpose)
	
	'''	
	max0 = (dt[3]+dt[4]+dt[5])/3.
	if(max0 == 0.):
		max0 = 0.0000000000001
	radii = array('f',[])
	ang = array('f',[])
	maxdiff = 0.
	OptDirection=0.
	OptVelocity= 99999.
	print ("time differences:",dt)
	for WaterVelocity in (1.42,1.44,1.46,1.48,1.50,1.52,1.54,1.56,1.58,1.60,1.62,1.64,1.66,1.68,1.70,1.72,1.74,1.76,1.78):
#		print('Velocity refine:',WaterVelocity)
		for i in range (720):
			theta = 2.*np.pi*float(i)/720.		
			dt1,dt2 = CompDt(d1,d2,theta,WaterVelocity)
			dt3 = dt2-dt1
			diff = math.sqrt((dt1-dt[0])*(dt1-dt[0])+(dt2-dt[1])*(dt2-dt[1])+(dt3-dt[2])*(dt3-dt[2]))

			diff = 1-math.erf(diff/max0)
#			print(diff)
			if (diff > maxdiff):
				maxdiff = diff
				OptDirection = theta
				OptVelocity = WaterVelocity
#				print "Direction", theta, WaterVelocity, maxdiff
	for i in range (720):
		theta = 2.*np.pi*i/720.		
		dt1,dt2 = CompDt(d1,d2,theta,OptVelocity)
		dt3 = dt2-dt1
		diff = math.sqrt((dt1-dt[0])*(dt1-dt[0])+(dt2-dt[1])*(dt2-dt[1])+(dt3-dt[2])*(dt3-dt[2]))
		diff = 1-math.erf(diff/max0)
		radii.append(diff)
		ang.append(theta)
	return OptDirection, OptVelocity, ang, radii, maxdiff

def CompmaxDt(sfile, triad):	
	'''
	Computes the maximum expected absolute propagation time difference
	between sta1 and sta2 (dt1), and sta1 and st3 (dt2) in the triad 
	formed by sta1, sta2, and sta3
	''' 

	d1,d2 = ReadSite(sfile, triad, 0)
	dt1max=0.
	dt2max=0.
	WaterVelocity=1.42
	for i in range (720):
#		print("i=",i)
		theta=2.*np.pi*i/720.
		dt1,dt2=CompDt(d1,d2,theta,WaterVelocity)
		if abs(dt1) > dt1max:
			dt1max=abs(dt1)
		if abs(dt2) > dt2max:
			dt2max=abs(dt2) 
	return dt1max, dt2max 

def CompDt(d1, d2, d, WaterVelocity):
	'''
	Computes two time differences given two input vectors giving the triad geometry in East-North coordinate system.
	The vectors are (sta2-sta1) and (st3-sta1), where sta1, sta2, and sta3 are the hydrophone elements
	of the triad. 
	The time differences are the difference in propagation time of an acoustic wave in the water given a velocity of 1.5km/s.
                       
                        d1										
	sta1 	O----------------->O sta2	
	        |
	     d2 |	        
	        |
		v
	sta3    O

	For instance if sta2 is one kilometer to the east of sta1, and sta3 is 0.5 kilometers to the south of sta1, 
	the first vector is (1.,0), and the second vector is (0.,-0.5)

	Input:
		d1: first vector (tuple with two floats)
		d2: second vector (tuple with two floats)
		d:  direction in radian from East. Counterclockwise is positive.
	Output:
		dt1, dt2: tuple containing the two time differences
	'''
#	print ("Watervelocity=", WaterVelocity, d1, d2, d)
	dot = d1[0]*math.cos(d)+d1[1]*math.sin(d)
	dt1 = dot/WaterVelocity
	dot = d2[0]*math.cos(d)+d2[1]*math.sin(d)
	dt2 = dot/WaterVelocity
	return dt1, dt2


def Semblance(theta, w1, w2, w3, samprate, d1, d2, verbose):
	semb = array('f',[])	
	for WaterVelocity in (1.46,1.48,1.50,1.52,1.54):
		for i in range(len(theta)):
			dot = d1[0]*math.cos(theta[i])+d1[1]*math.sin(theta[i])
			dt1 = dot*KMperDEG/WaterVelocity
			dot = d2[0]*math.cos(theta[i])+d2[1]*math.sin(theta[i])
			dt2 = dot*KMperDEG/WaterVelocity
			s = 0.
			span = len(w1) - int(samprate*max(abs(dt1),abs(dt2)))	
			if verbose > 2:
				print ("span=", span, "Length(w1) =", len(w1))
			for j in range(span):
				j1 = j + int(dt1*samprate)
				j2 = j + int(dt2*samprate)
				s = s + w1[j]+w2[j1]+w3[j2]
				if verbose > 2:
					print ("D1=", d1, "D2=",d2, "theta=", theta[i], "Dt=", dt1, dt2, s)
			semb.append(s/float(span))	
	#
# Normalize semblance
#
	avg = 0.
	for i in range(len(theta)):
		avg = avg + semb[i]
	avg = avg / float (len(theta))
	for i in range(len(theta)):
		semb[i] = semb[i] - avg
	m = 0.
	for i in range(len(theta)):
		if abs(semb[i]) > m : 
			m = abs(semb[i])
	for i in range(len(theta)):
		if m == 0. : 
			semb[i] = 0.
		elif m > 0. :
			semb[i] = 0.5*(semb[i]+m)/m

#
# Get bearing of line between first and second station from North
#
	if(d1[0] != 0.):
		theta0 = math.atan(d1[0]/d1[1])
		print ("Angle of first station, from North:", theta0)
	elif (d1[0] == 0.):
		theta0 = np.pi/2.
	if theta0 >= 0. :
		shft = int(theta0*len(theta)/2*np.pi)
	elif theta0 < 0.:
		shft = len(theta) - int(-1*theta0*len(theta)/2*np.pi)
	
	sem = shiftInPlace(semb, shft)
	return sem
		

def shiftInPlace(l, n):
	n = n % len(l)
	head = l[:n]
#    l[:n] = []
	l.extend(head)
	shl = l[n:]
	return shl


def CorrMaster(w, mast, length):
	''' 
	CorrMaster correlates a waveform w against a master waveform
	Input:
		w: input waveform to be correlated with the master waveform mast
		mast: master waveform
		length: length (in number of samples) of the correlation function to be computed.
			The total length will be twice this.
	Output:
		Correlation function of length 2*length seconds 	
	''' 
	lmast = len(mast)
	lw = len(w)
	crossco = array('f',[])
	
# Compute the energy in master waveform

	en_master = 0.
	for i in range (lmast):
		en_master += mast[i]*mast[i]
	en_master = math.sqrt(en_master)	
#		
# Compute the energy in waveform 
#
	en_w = 0.
	for j in range (len(w)):
		en_w += w[j]*w[j]
	en_w = math.sqrt(en_w)	
#		
# Compute the energy in waveform 
#	
	for k in range (length):
		cross = 0.
		if (k%100 == False): 
			print (k)
		for l in range((len(w)-k)):
			cross += w[k+l]*mast[l]
		cross = cross/en_w/en_master
		crossco.append(cross)
	array.reverse(crossco)
	for k in range (length):
		cross = 0.
		if (k%100 == False):
			print (k)
		for l in range(1,(len(w)-k)):
			cross += mast[k+l]*w[l]
		cross = cross/en_w/en_master
		crossco.append(cross)
	array.reverse(crossco)
	return crossco

def Across(w, length):
	''' 
	Across computes a normalized autocorrelation. The normalization allows to take into account the smaller correlation 
	length at larger lag times
	Input:
		w: input waveform to be correlated with the master waveform mast
		length: length (in number of samples) of the correlation function to be computed.
			The total length will be twice this.
	Output:
		AutoCorrelation function of length 2*length seconds 	
	''' 

	lw = len(w)
	across = array('f',[])
#		
# Compute the energy in waveform 
#	
	for k in range (length):
		cross = 0.
		if (k%100 == False): 
			print (k)
		for l in range(lw-k):
			cross += w[k+l]*w[l]
		across.append(cross)
	lag0 = across[0]
	for k in range (length):
		if (k%100 == False):
			print (k)
		across[k] = across[k]/lag0
	return across

def TimeOfMax(w):
	'''
	Returns the index of the absolute maximum of the input function
		Input:
			w: input waveform
		Output:
			maxindex: Index of the maximum absolute value
	'''
	maximum=0.
	maxindex=0
	for i in range(len(w)):
		if(abs(w[i]) > maximum):
			maximum=abs(w[i])
			maxindex=i
	return maxindex

def MakeDetect(Det, w, snr, samprate, t_group, sta, shift):
	det = []
	spectim=[]
	init=0
	det_num=0
	for i in range(1,len(w)):
		t0=float(i)/samprate
		if(w[i] >= snr and w[i] > w[i-1] and init==0):
			init=1
			dw = (snr-w[i-1])/(w[i]-w[i-1])
			if dw != 0. :
				t = t0+shift+dw/samprate
			ndet = len(det)
#			print "ndet", ndet, det
			if((ndet==0) or (ndet > 0 and (t-det[ndet-1].time) > t_group)):
				newdet=[sta,det_num,t,0.,0.,0.,False,0.,0.,0.,spectim]
				det_num+=1
				print (t0, w[i-1], w[i])
				print (newdet,'\n')
				det.append(Detection._make(newdet))	
		if(w[i] < snr ):
			init=0	
	
	return det	
	
def CompStaLta (L_ta, S_ta, samprate, w):
	'''
	Computes a short-term-average (STA) over long-term-average (LTA) function
		Input:
			L_sta: 		length in seconds of the long-term-average window
			S_sta:	 	length in seconds of the short-term-average
			samprate: 	sample rate in number of samples per second
			w:     		input waveform
		Output:
			returns the STA/LTA function of length len(w)-L_ta.
			The first sample is computed at time L_ta from the beginning of the w function
	'''

	StaLta = array('f',[]) 
	W = array('f',[])
	for j in range(len(w)):
		W.append(w[j]*w[j])
	Lta = 0.
	Sta = 0.
	Lshift = int(L_ta*samprate)
	Sshift = int(S_ta*samprate)
	L_length = float(L_ta*samprate)
	S_length = float(S_ta*samprate)
	for j in range(Lshift):
		Lta += W[j]
		if j > (Lshift-Sshift):
			Sta += W[j]

	Lta = Lta / L_length
	Sta = Sta / S_length
	
	StaLta.append(Sta/Lta)
	j=0
	while j <  (len(w)-Lshift):	 
		Lta = Lta+(W[Lshift+j] - W[j])/L_length
		Sta = Sta+(W[Lshift+j] - W[j+Lshift-Sshift+1])/S_length 
		StaLta.append(Sta/Lta)
		if(Sta/Lta < 0.) :
			print("Warning !! negative STA - LTA ratio",j, Sta,Lta)
		j+=1
	del W
	return StaLta


def ReadSrate(wfdisc, station, verbose):
	'''
	Reads sample rate from the tuples 
        
	'''

	fwd = open(wfdisc,"r")

#
# Read wfdisc tuples. Loop around the lines
#
	while True:
	
		line = fwd.readline()
		if not line:
			break
		pass #Parse Wfdisc tuple    
		if verbose > 3: 
			print (line)
#
# Get rid of all double white space in the string
#
		while '  ' in line:
			line=line.replace("  "," ")
		if verbose > 3 :
		   	print (line)
		values=line.split(" ")
#
# Print the whole tuple if the verbosity is higher than 2
#
		if station != values[0] :
				continue
		got_station = True
		if verbose > 2:
			print ("Sta=",values[0],'\n', "Chan=",values[1],'\n',"Time=",values[2],'\n',"Wfid=",values[3],"\n","Chanid=",values[4],"\n","jdate=",values[5],"\n","endtime=",values[6],"\n","nsamp=",values[7],"\n","samprate=",values[8],"\n","calib=",values[9],"\n","calper=",values[10],"\n","instype=",values[11],"\n","segtype=",values[12],"\n","datatype=",values[13],"\n","clip=",values[14],"\n","dir=",values[15],"\n","dfile=",values[16],"\n","Byte offet foff=",values[17],"\n","commid=",values[18],"\n","lddate=",values[19],'\n',"HH:MM:SS=",values[20],"\n")


		srate      = float(values[8])

	if(got_station):
		return srate
	else:
		return -999.
def Readflat(filename, skip, verbose):
	'''
	Reads waveforms in a flat ASCII file
	'''

	print(filename)
	fwd = open(filename, "r")

	wfm=array('f',[])
	time=array('f',[])
	samprate = 0.01
	for i in range(skip):
		line = fwd.readline()
		if verbose > 2 :
			print(line)

	print(line)
	tim = 0.
	while True:
		line = fwd.readline()
		if not line:
			break
#
# Get rid of all double white space in the string
#
		while '  ' in line:
			line=line.replace("  "," ")
		if(verbose > 2) :
			print(line)

		values = line.split(" ")

		if(len(values) > 1):
			if verbose > 2:
				print ("Time",values[1],'\n', "Count=",values[2],"\n")
			tim = float(values[1])
			amp = float(values[2])
		else:
			amp = float(values[0])
			tim += samprate
		if verbose > 2:
			print("Time=",tim, "Amplitude=", amp)	
		wfm.append(amp)
		time.append(tim)
#
# Remove average from the waveform. Compute average ignoring zero values
#
	Sum = 0.
	for i in range (len(wfm)):
			Sum = Sum + wfm[i]
	avg = 0.
	lzero = 0
	for i in range (len(wfm)):
			if( wfm[i] == 0.):
				lzero += 1
	if wfm:
		avg = Sum / (float(len(wfm)-lzero))
#
# Set zero values to average
#
	for i in range (len(wfm)):
		if(wfm[i] == 0.):
			wfm[i] = avg

	for i in range (len(wfm)):
		wfm[i] = wfm[i]-avg
	
	return time, wfm

def ReadWfm(wfm, wfdisc, start_epoch, end_epoch, station, verbose):
	'''
	Reads waveform in CSS o4 format between start_epoch and end_epoch times in epoch seconds
	
		Input:
			wfdisc: 		WFDISC tuples ASCII file 
			start_epoch:	 	start epoch time in seconds
			end_epoch:     		end epoch time in seconds
			station:		Station name
			verbose:		verbosity level. The higher the more verbose (max=4)
		Output:
			returns the sample rate in a single float, and the waveform in an array of floats.
	'''


	fwd = open(wfdisc,"r")

#
# Read wfdisc tuples. Loop around the lines
#

#	wfm=array('f',[])
	start_time = float(start_epoch)

	first_time = True
	while True:
		line = fwd.readline()
		if not line:
        		break
		pass #Parse Wfdisc tuple    
		if verbose > 3: 
			print (line)
#
# Get rid of all double white space in the string
#
		while '  ' in line:
			line=line.replace("  "," ")
		if verbose > 3 :
	    		print (line)
		values=line.split(" ")
	#
	# Print the whole tuple if the verbosity is higher than 2
	#
		if verbose > 2:
			print ("Sta=",values[0],'\n', "Chan=",values[1],'\n',"Time=",values[2],'\n',"Wfid=",values[3],"\n","Chanid=",values[4],"\n","jdate=",values[5],"\n","endtime=",values[6],"\n","nsamp=",values[7],"\n","samprate=",values[8],"\n","calib=",values[9],"\n","calper=",values[10],"\n","instype=",values[11],"\n","segtype=",values[12],"\n","datatype=",values[13],"\n","clip=",values[14],"\n","dir=",values[15],"\n","dfile=",values[16],"\n","Byte offet foff=",values[17],"\n","commid=",values[18],"\n","lddate=",values[19],'\n',"HH:MM:SS=",values[20],"\n")
	    
		if station != values[0] :
			continue

		tuptime    = float(values[2])
		tupendtime = float(values[6])
		foff       = int(values[17])
		nsamp      = int(values[7])
		wfid       = int(values[3])
		srate      = float(values[8])
		nint       = 1+int((end_epoch - start_epoch)*srate)
		N          = nsamp
		directory  = values[15]
		dfile      = values[16]
		datatype   = values[13]
		num_samp=0

	#   	
	# Set the array to zero if first time in
	#
		
			

		if nint < nsamp and nint > 0 :
			N = nint

		if verbose > 2:
			print ("Tuple Start Time=", tuptime,'\n', "Tuple End Time=",tupendtime,'\n')
			print ("nint=", nint, "nsamp=", nsamp)
		if  start_epoch <  tuptime and end_epoch <= tuptime:
			continue
		if  start_epoch >=  tupendtime and start_epoch > tupendtime and end_epoch > tupendtime:
			continue

		wfile = directory+"/"+dfile 
		f   = open(wfile,"rb")
		if  ((abs(start_epoch - tuptime) < 1./srate) or (start_epoch > tuptime)) and ((end_epoch < tupendtime) or (abs(end_epoch-tupendtime) < 1./srate)):

# All the requested waveform is included in this tuple.

			num_samp = int((start_epoch-tuptime)* srate)
			if(datatype == 's4'):
				foff = foff + 4*num_samp
			if(datatype == 's3'):
				foff = foff + 3*num_samp
			if (verbose > 1):
				print ("All requested Waveform in tuple:", wfid,"station:",station,'\n')
				print ("Reading", N, "samples",(N-1)/srate," seconds",'\n') 
				print ("Tuple Start Time = ", tuptime,'\n', "Tuple End Time = ",tupendtime,'\n')
				print ("start_epoch = ", start_epoch, "end_epoch = ", end_epoch)
				print ("Offset = ", foff)		
			
			readwfm = readNpoints(f,foff,N,verbose,datatype)
			for i in range(N):
				wfm[i]=readwfm[i]
			start_time = start_epoch
			break	

		if  ((abs(start_epoch - tuptime) < 1./srate) or (start_epoch > tuptime)) and start_epoch < tupendtime :

		    # The beginning of the waveform is in this tuple.
			num_samp = int((start_epoch-tuptime)* srate)
			if(datatype == 's4'):
				foff = foff + 4*num_samp
			if(datatype == 's3'):
				foff = foff + 3*num_samp
			N = int((tupendtime-start_epoch)*srate) 
			if (verbose > 1):
				print ("Beginning of Waveform in tuple:", wfid,"station:",station,'\n')
				print ("Reading", N, "samples",(N-1)/srate," seconds",'\n')
				print ("Tuple Start Time=", tuptime,'\n', "Tuple End Time=",tupendtime,'\n')
				print ("start_epoch=", start_epoch, "end_epoch=", end_epoch)
			readwfm = readNpoints(f,foff,N,verbose,datatype)
			for i in range(N):
				wfm[i] = readwfm[i]
			continue
	
		if  start_epoch < tuptime and end_epoch >= tuptime:

			# The end or continuation of the waveform is in this tuple.
			if(verbose > 1):
				print ("Station:", values[0], "Continuation of waveform in tuple:", wfid,"station:",station,'\n')
				print ("Start, Start_tuple, End",start_epoch, tuptime, end_epoch)
				print ("Tuple Start Time=", tuptime,'\n', "Tuple End Time=",tupendtime,'\n')
				print ("start_epoch=", start_epoch, "end_epoch=", end_epoch)
			N = nsamp
			num_samp += N
			if nint < nsamp and nint > 0 :
				N = nint 
			if end_epoch < tupendtime:
				N = int((end_epoch-tuptime)*srate)
		
			if (verbose > 1):
				print ("Reading", N, "samples",(N-1)/srate," seconds",'\n')
			start_time = tuptime
			readwfm = readNpoints(f,foff,N,verbose,datatype)
			start_samp = int(srate*(tuptime-start_epoch))
			for i in range (N):
				wfm[start_samp+i] = readwfm[i]

		f.close()
	fwd.close()
#
# Remove average from the waveform. Compute average ignoring zero values
#
	Sum = 0.
	for i in range (len(wfm)):
			Sum = Sum + wfm[i]
	avg = 0.
	lzero = 0
	for i in range (len(wfm)):
			if( wfm[i] == 0.):
				lzero += 1
	if wfm:
		avg = Sum / (float(len(wfm)-lzero))
#
# Set zero values to average
#
	for i in range (len(wfm)):
		if(wfm[i] == 0.):
			wfm[i] = avg

	for i in range (len(wfm)):
		wfm[i] = wfm[i]-avg

	return srate, start_time, num_samp 


def sext24(d):
	if ord(d[2]) & 0x80:
		return d+'\xff'
	else:
		return d+'\x00' 


def readNpoints(f, foff, N, verbose, datatype):
	'''
	Read N points from file at offset foff 
	Does the unpacking of the 4 byte string into two byte integers
	
	Input: 
		f:	Name of the file containing the four byte integer to be read 
		foff:	Offset in number of bytes from the beginning of the file where to start reading it 
	 	N:	Number of points to read 
	
	 Output: 
	 	wfm:	Array containing the N points to be read
	
	'''


#   	
# Set the array to zero if first time in
#
	first_time = True
	wfm = array('f',[])
	if(first_time):
		first_time = False
		for i in range(N):
			wfm.append(0.)	
	f.seek(foff,0)	
   
	
	for i in range(N) :
		if(datatype == 's4'):
			ndstring = f.read(4)
			ndbyte = np.ndarray(shape = (2,), dtype='>i2', buffer=ndstring)
			if verbose > 3:
				print (ndbyte[0], ndbyte[1])
			wfm[i] = float(ndbyte[1])		
			if verbose > 3:
				print ("Sample", i, wfm[i-1], "\n")
		if(datatype == 's3'):
			ndstring = f.read(3)
			ndbyte = sext24(ndstring)
			upack = struct.unpack(">%sl" % 1 , ndbyte) 
			if verbose > 3:
				print (upack[0])
			wfm[i]=float(upack[0])		
			if verbose > 3:
				print ("Sample", i, wfm[i-1], "\n")
	return wfm

def ReadWfm_DB(wfm, account, wfdisc, start_epoch, end_epoch, station, verbose):
	'''
	Reads waveform in CSS o4 or o3 format between start_epoch and end_epoch times in epoch seconds
	
		Input:  wfm:                    buffer containing waveform in output
                        account:                database account
		        wfdisc:                 name of WFDISC table
			start_epoch:	 	start epoch time in seconds
			end_epoch:     		end epoch time in seconds
			station:		Station name
			verbose:		verbosity level. The higher the more verbose (max=4)
		Output:
			
			returns the sample rate in a single float, and the waveform in an array of floats.
	'''

	con = cx_Oracle.connect(account)
	if verbose > 3: 
		print ('Version:',con.version)
	

#
# Read wfdisc tuples. Loop around the lines
#

	
	query = fm.format(' Select sta, chan, time, wfid, chanid, jdate, endtime, nsamp, samprate, calib, calper, dir,  dfile, foff from {0} where time < {1} and endtime > {2} and sta=\'{3}\' ', wfdisc, start_time, start_time, sta)
	cur.execute(query)

	while True:
		
		for result in cur:
			Ntup = Ntup+1
			if verbose > 3: 
				print (Ntup, result)
			sta, chan, tuptime, wfid, chanid, jdate, tupendtime, nsamp, srate, calib, calper, directory,  dfile, foff  = result
			nint       = 1+int((end_epoch - start_epoch)*srate)
			N          = nsamp
			
			num_samp=0

	#   	
	# Set the array to zero if first time in
	#
		
			

			if nint < nsamp and nint > 0 :
				N = nint

			if verbose > 2:
				print ("Tuple Start Time=", tuptime,'\n', "Tuple End Time=",tupendtime,'\n')
				print ("nint=", nint, "nsamp=", nsamp)
			if  start_epoch <  tuptime and end_epoch <= tuptime:
				continue
			if  start_epoch >=  tupendtime and start_epoch > tupendtime and end_epoch > tupendtime:
				continue

			wfile = directory+"/"+dfile 
			f   = open(wfile,"rb")
			if  ((abs(start_epoch - tuptime) < 1./srate) or (start_epoch > tuptime)) and ((end_epoch < tupendtime) or (abs(end_epoch-tupendtime) < 1./srate)):

# All the requested waveform is included in this tuple.

				num_samp = int((start_epoch-tuptime)* srate)
				if(datatype == 's4'):
					foff = foff + 4*num_samp
				if(datatype == 's3'):
					foff = foff + 3*num_samp
				if (verbose > 1):
					print ("All requested Waveform in tuple:", wfid,"station:",station,'\n')
					print ("Reading", N, "samples",(N-1)/srate," seconds",'\n') 
					print ("Tuple Start Time = ", tuptime,'\n', "Tuple End Time = ",tupendtime,'\n')
					print ("start_epoch = ", start_epoch, "end_epoch = ", end_epoch)
					print ("Offset = ", foff)		
			
				readwfm = readNpoints(f,foff,N,verbose,datatype)
				for i in range(N):
					wfm[i]=readwfm[i]
				start_time = start_epoch
				break	

			if  ((abs(start_epoch - tuptime) < 1./srate) or (start_epoch > tuptime)) and start_epoch < tupendtime :

		    # The beginning of the waveform is in this tuple.
				num_samp = int((start_epoch-tuptime)* srate)
				if(datatype == 's4'):
					foff = foff + 4*num_samp
				if(datatype == 's3'):
					foff = foff + 3*num_samp
				N = int((tupendtime-start_epoch)*srate) 
				if (verbose > 1):
					print ("Beginning of Waveform in tuple:", wfid,"station:",station,'\n')
					print ("Reading", N, "samples",(N-1)/srate," seconds",'\n')
					print ("Tuple Start Time=", tuptime,'\n', "Tuple End Time=",tupendtime,'\n')
					print ("start_epoch=", start_epoch, "end_epoch=", end_epoch)
				readwfm = readNpoints(f,foff,N,verbose,datatype)
				for i in range(N):
					wfm[i] = readwfm[i]
				continue
	
			if  start_epoch < tuptime and end_epoch >= tuptime:

			# The end or continuation of the waveform is in this tuple.
				if(verbose > 1):
					print ("Station:", values[0], "Continuation of waveform in tuple:", wfid,"station:",station,'\n')
					print ("Start, Start_tuple, End",start_epoch, tuptime, end_epoch)
					print ("Tuple Start Time=", tuptime,'\n', "Tuple End Time=",tupendtime,'\n')
					print ("start_epoch=", start_epoch, "end_epoch=", end_epoch)
				N = nsamp
				num_samp += N
				if nint < nsamp and nint > 0 :
					N = nint 
				if end_epoch < tupendtime:
					N = int((end_epoch-tuptime)*srate)
		
				if (verbose > 1):
					print ("Reading", N, "samples",(N-1)/srate," seconds",'\n')
				start_time = tuptime
				readwfm = readNpoints(f,foff,N,verbose,datatype)
				start_samp = int(srate*(tuptime-start_epoch))
				for i in range (N):
					wfm[start_samp+i] = readwfm[i]

	
#
# Remove average from the waveform. Compute average ignoring zero values
#
	Sum = 0.
	for i in range (len(wfm)):
			Sum = Sum + wfm[i]
	avg = 0.
	lzero = 0
	for i in range (len(wfm)):
			if( wfm[i] == 0.):
				lzero += 1
	if wfm:
		avg = Sum / (float(len(wfm)-lzero))
#
# Set zero values to average
#
	for i in range (len(wfm)):
		if(wfm[i] == 0.):
			wfm[i] = avg

	for i in range (len(wfm)):
		wfm[i] = wfm[i]-avg

	return srate, start_time, num_samp 


def sext24(d):
	if ord(d[2]) & 0x80:
		return d+'\xff'
	else:
		return d+'\x00' 


def readNpoints(f, foff, N, verbose, datatype):
	'''
	Read N points from file at offset foff 
	Does the unpacking of the 4 byte string into two byte integers
	
	Input: 
		f:	Name of the file containing the four byte integer to be read 
		foff:	Offset in number of bytes from the beginning of the file where to start reading it 
	 	N:	Number of points to read 
	
	 Output: 
	 	wfm:	Array containing the N points to be read
	
	'''


#   	
# Set the array to zero if first time in
#
	first_time = True
	wfm = array('f',[])
	if(first_time):
		first_time = False
		for i in range(N):
			wfm.append(0.)	
	f.seek(foff,0)	
   
	
	for i in range(N) :
		if(datatype == 's4'):
			ndstring = f.read(4)
			ndbyte = np.ndarray(shape = (2,), dtype='>i2', buffer=ndstring)
			if verbose > 3:
				print (ndbyte[0], ndbyte[1])
			wfm[i] = float(ndbyte[1])		
			if verbose > 3:
				print ("Sample", i, wfm[i-1], "\n")
		if(datatype == 's3'):
			ndstring = f.read(3)
			ndbyte = sext24(ndstring)
			upack = struct.unpack(">%sl" % 1 , ndbyte) 
			if verbose > 3:
				print (upack[0])
			wfm[i]=float(upack[0])		
			if verbose > 3:
				print ("Sample", i, wfm[i-1], "\n")
	return wfm

def WriteSoundFile(sound, prefix, samprate, N):
	'''
	Writes a .wav file containing audible compressed at N times the original sample rate
	'''
	soundfilename = prefix+'.wav'
		
	file = glob.glob(soundfilename)
	for filename in file:	
		os.remove(filename)
	soundfile = open(soundfilename,'w')
	soundfile.write('')
	soundfile.close()
	soundfile = open(soundfilename,'w')
	sf = wave.open(soundfile,'w')
	nchannels = 1
	sampwidth = 2
	framerate = int(samprate*N)
	nframes = len(sound)
	comptype = "NONE"
	compname = "not compressed"
	maxwf=0.
	for i in range(len(sound)):
		if abs(sound[i]) > maxwf:
			maxwf = abs(sound[i])
		continue	

	optimax = 30000./maxwf

	for i in range(len(sound)):
		sound[i] = int(float(sound[i]*optimax))

	sf.setparams((nchannels, sampwidth, framerate, nframes, comptype, compname))

	for i in range(len(sound)):
		sf.writeframes(struct.pack('h',sound[i]))
	sf.close()
	soundfile.close()	



