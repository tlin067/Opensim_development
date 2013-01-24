Tom Lintern
11/1/2013

v1b Have included coupledbushingforce.cpp. This allows user to specify full 6x6 stiffness and damping matrix. There are no coupling terms however. 

Need to add function based descriptions of these, as currently only linear coef's are allowed. 

v1c Have extended the stiffness matrix and q vector in "coupledBushing: class. Am working on creating properties which I can use to define different values within the stiffness matrix. Trying to use Matt deMers funciton bushing code to do so.

v1d Have created a functionBasedCoupled bushing force function.
i.e. ... 
T1 = 0
T2 = 0
T3 = f(theta) + f(x) + f(y)
Fx = f(theta) + f(x) + f(y)
Fy = f(theta) + f(x) + f(y)
Fz = 0

