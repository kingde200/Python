import math
import numpy 


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
     

    return A,FF



