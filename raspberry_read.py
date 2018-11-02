# -*- coding: utf-8 -*-
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import matplotlib.pyplot as plt

start_1="2018-10-25 22:56:00"
starttime=UTCDateTime(start_1)
endtime=starttime+900
client = Client(base_url='https://fdsnws.raspberryshakedata.com/')
waveform0 = client.get_waveforms('AM','R8A6C','00','EHZ',starttime,endtime)
waveform1 = client.get_waveforms('AM','R898B','00','SHZ',starttime,endtime)
waveform2 = client.get_waveforms('AM','R655B','00','SHZ',starttime,endtime)

print(waveform2)
waveform0.filter("bandpass",freqmin=0.05,freqmax=2.)
waveform1.filter("bandpass",freqmin=0.05,freqmax=2.)
waveform2.filter("bandpass",freqmin=0.05,freqmax=2.)
start_p="2018-10-25 22:57:00"
startplot=UTCDateTime(start_p)
waveform0.plot(starttime=startplot,outfile='Zakinthos0.png')
waveform1.plot(starttime=startplot,outfile='Zakinthos1.png')
waveform2.plot(starttime=startplot,outfile='Zakinthos2.png')
