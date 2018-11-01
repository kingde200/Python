def my_function(A):
    n = len(A)
    B = A
    for i in range(n):
      B[i]= A[i]/2
    for x in B:
      print(x)
    return B
A = [1.,2.,3.]
C = my_function(A)
for x in C:
   print(x)
