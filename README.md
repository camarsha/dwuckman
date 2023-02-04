## Optical Model Calculations

This is a side project that hopefully will be fleshed out slowly over the coming months. The idea is to 
have a lightweight optical model code that can be called easily by other programs. Right now it is in the testing 
stage.

## To Do

1) Confirmation of cross sections over a broader range of energies/masses/etc.

2) Test for convergence to stop partial wave loop (might be difficult with parallel implementaton).

3) Refactor the code a bit so that phase shifts are the output of the main l-loop. Implement total cross sections.

4) Test for neutrons

5) Probably need to implement COULFG to replace the gsl library coulomb functions, which have known issues.

Test
