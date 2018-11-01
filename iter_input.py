def my_function(A,n):
    B = []
    for i in range(len(A)):
       B.append(0.)
    k = len(B)
    
    for i in range(k): 
         for it in range(n):
             print(it)
             if (it == 0):
               B[i] = 1./A[i]
             else:
               B[i] = B[i]/A[i]
    return B

A = [1.,2.,3.,4.]
for x in A:
    print(x)
C = my_function(A,4)
for x in C:
     print(x)
   
