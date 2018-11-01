import math
import numpy as np
from scipy.fftpack import fft,ifft

def poles_zeros(P,Z,FF):
    A = []
   
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
     

    return A


def impulse_response(A,ny,t0):

	B = []
        N = len(A)
        s = ny/(N-1) 
	for i in range(N):
	   arg = -2.*math.pi*complex(0.,1.)*t0*i*s
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

        return Im



