import math
import cmath
import numpy 
import matplotlib.pyplot as plt
import poles_zeros as pz

P = [complex(-0.037,-0.037),complex(-0.037,0.037),complex(-188.8279,195.5396),complex(-188.8279, -201.8228),complex(-259.2216,719.6446),complex(-259.2216,-719.6446)]
Z = [complex(0,0),complex(0,0)]
FF = [0.0001,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1, 0.5, 0.8,1.,2.,10.,20.,40.,50.,60.,70.,80.,90.,100.] 
A = pz.poles_zeros(P,Z,FF)
A_amp = []
ph = []
for i in range(len(A)):
   A_amp.append(abs(A[i])*4.395259e+10)
   ph.append(cmath.phase(A[i]))
plt.subplot(211)
plt.title('Frequency response of Guralp CMG3TB seismometer')
plt.plot(FF,A_amp)
plt.yscale('log')
plt.ylabel('Amplitude(V/m/s)')
plt.xscale('log')
plt.xlabel('Frequency(Hz)')
plt.subplot(212)
plt.plot(FF,ph)
plt.ylabel('Phase(radians)')
plt.xscale('log')
plt.xlabel('Frequency(Hz)')
plt.show()

