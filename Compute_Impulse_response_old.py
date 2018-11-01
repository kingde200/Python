import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
import poles_zeros as pz
from scipy.fftpack import fft,ifft

P = [complex(-0.037,-0.037),complex(-0.037,0.037),complex(-188.8279,195.5396),complex(-188.8279, -201.8228),complex(-259.2216,719.6446),complex(-259.2216,-719.6446)]
Z = [complex(0,0),complex(0,0)]


FF = []
TT = []
N = 1025
S = 20./N
T = 1./40.
T0 = 0.

for i in range(N):
   FF.append(i*S)
for i in range(N-1):
   TT.append(-25.6+i*T)
for i in range(N-1):
   TT.append(i*T)


#
#Compute Frequency domain response
#
A = pz.poles_zeros(P,Z,FF)
B = []
for i in range(N):
   arg = -2.*math.pi*complex(0.,1.)*T0*i*S
   shift = np.exp(arg)
   B.append(A[i]*shift)
for i in range(N-1):
   B.append(np.conj(B[N-i-1]))
Imp = ifft(B)
K = len(Imp);
Im = []
for i in range(K/2):
   Im.append(Imp[K/2+i]) 
for i in range(K/2):
   Im.append(Imp[i])
A_amp = []
ph = []
for i in range(len(A)):
   A_amp.append(abs(A[i])*4.395259e+10)
   ph.append(cmath.phase(A[i]))
plt.subplot(311)
plt.plot(TT,Im)
plt.title('Impulse response')
plt.subplot(312)
plt.plot(FF,A_amp)
plt.yscale('log')
plt.ylabel('Amplitude(V/m/s)')
plt.xscale('log')
plt.title('Amplitude Spectrum')
plt.subplot(313)
plt.plot(FF,ph)
plt.ylabel('Phase(radians)')
plt.xscale('log')
plt.xlabel('Frequency(Hz)')
plt.title('Phase Spectrum')

plt.show()
