import math
import cmath
import numpy 
import matplotlib.pyplot as plt

def poles_zeros(P,Z):
    A = []
    FF = [0.0001,0.001, 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1, 0.5, 0.8,1.,2.,10.,20.,40.,50.,60.,70.,80.,90.,100.] 
    for k in range(len(FF)):
      F = FF[k]
      arg = complex(0.,1.)*math.pi*F
      print('ff')
      Num = 1.
      for i in range(len(Z)):
          zz = Z[i]
          print('Zero', i, zz.real)
          Num = Num*(arg-zz)
      Den = 1.
      for j in range(len(P)):
          print('Pole',j, P[j].real)
          Den = Den*(arg-P[j])
  #Compute H(F) 
      A.append(Num/Den)
  #Compute amplitude and phase
     

    return A,FF


P = [complex(-48.8255,306.345),complex(-48.8255,-306.345),complex(-222.351,118.067),complex(-222.351,-118.067),complex(-3.32882,1.98402),complex(-3.32882,-1.98402)]
Z = [complex(-4.51187,308.022),complex(-4.51187,-308.022),complex(0,0),complex(0,0)]

A, Freq = poles_zeros(P,Z)
A_amp = []
ph = []
for i in range(len(A)):
   A_amp.append(abs(A[i])*65354.1)
   ph.append(cmath.phase(A[i]))
plt.subplot(211)
plt.title('Frequency response of an instrument')
plt.plot(Freq,A_amp,)
plt.yscale('log')
plt.ylabel('Amplitude(V/m/s)')
plt.xscale('log')
plt.xlabel('Frequency(Hz)')
plt.subplot(212)
plt.plot(Freq,ph)
plt.ylabel('Phase(radians)')
plt.xscale('log')
plt.xlabel('Frequency(Hz)')
plt.show()
