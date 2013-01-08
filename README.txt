Tom Lintern
8/1/2013

Memory leak fixed when chanigng force properties in model.

Done by accessing existing model forces and casting them to be bushing or limit forces respectively.

Code is set up to reproduce results from CMBBE conference abstract.