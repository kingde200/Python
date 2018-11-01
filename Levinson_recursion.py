import math
def levinson_recursion(R,N):
    A = []
    if (N < 3) : 
       print("Use a value larger than 2 for number of iterations\n")
       return A
# Compute A[0]
    A.append(1.)
    nr = len(R)

# Compute A[1] and v 
    A.append(-1.*R[1]/R[0])
    v = R[0]-R[1]*R[1]/R[0]
    nit = N-2
    for k in range(nit):
         e = 0.;
         for i in range(k+1):
             e = e + A[i]*R[2+k-i] 

         for i in range(k):
             A[i+1] = A[i+1]*(1.-e/v)
         A.append(-1.*e/v)      
         vp = v*(1.-e*e/v/v)    
         v = vp

    return A


R = [1.,0.98,0.95,0.90,0.8,0.7]
A = levinson_recursion(R,3)
for x in A:
  print(x)
