import math
import cmath
import numpy 
import matplotlib.pyplot as plt
import poles_zeros as pz
import impulse_response as ir

P = [complex(1.,0)]
Z = [complex(1.,0)]
N = 1025
S = 20./N
T = 1./40
T0 = 0.
A = imp.impulse_response(P,Z,N)
A_amp = []
ph = []
for i in range(len(A)):
   A_amp.append(abs(A[i]))
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
