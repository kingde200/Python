import math
def my_function(A):
    n = len(A)
    B = []
    for i in range(n):
      B.append(math.cos(A[i]))
    for x in B:
      print(x)

    return B


A = [1.,2.,3.]
C = my_function(A)
for x in C:
   print(x)
