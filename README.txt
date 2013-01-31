Tom Lintern
11/1/2013

Have added in MAtt deMers plugin to create FunctionBasedBushing forces.

I have created function ArbNonLinear.h which is a function that returns y = Ax^n + B where A, n, B are free parameters.

I have set up an optimisation to run where

T[0] = 0 -> functions are zero by default.
T[1] = 0
T[2] = f(A,n,B,theta) + f(A,n,B,x) + f(A,n,B,y)
f[1] = f(A,n,B,theta) + f(A,n,B,x) + f(A,n,B,y)
f[2] = f(A,n,B,theta) + f(A,n,B,x) + f(A,n,B,y)    nb all A B and n's are different.
f[3] = 0

Code has also been deployed to hpc.

Plan to try optimisation without cross terms before including them to see if result improves.



