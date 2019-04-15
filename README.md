# MERK: Multirate Exponential Runge--Kutta methods # 

#### Rujeko Chinomona, Daniel R. Reynolds, and Vu Thai Luan ####
Department of Mathematics, Southern Methodist University

This repository contains MATLAB files for test problems highlighted in the paper, 
   "A new class of high-order methods for multirate differential equations" Vu Thai Luan, Rujeko Chinomona & Daniel R. Reynolds. 
   
Multirate Expontential Runge-Kutta (MERK) methods of orders three (`MERK3`), four (`MERK4`), and five (`MERK5`) are implemented on the additively split multirate problem: 

  ```y' = Ay + g(t,y)``` 

where ```Ay``` denotes the fast linear part (cheap) and ```g(t,y)``` denotes the slow and expensive part. In addition, we also run tests with the Multirate Infinitesimal Step method (`MIS-KW3`, J. Wensch, O. Knoth and A. Galant, BIT Numer. Math., 2009) for comparison.

The test problems are divided into two categories as discussed in the paper. The first category has problems that are concerned about stiffness in the fast component and an expensive slow component. For this set of problems, the micro step denoted *hfast* is held constant while the macro time step *h* and the time scale separation factor *m* vary. The brusselator test problem (`driver_brusselator.m`) together with a reaction diffusion (`driver_randd.m`) problem are considered. 

The second category is concerned with problems for which the fast component contributes significantly to the slow dynamics. We consider two test problems, one with bidirectional coupling (`driver_bidirectioncouple.m`) between fast and slow and another with only one directional coupling (`driver_onedirectioncouple.m`). For this set of problems, the "optimal" time scale separation factor *m* is held constant while the macro and micro time steps vary.

For all the tests, three figures are generated, one convergence plot and two efficiency plots. For the brusselator problem, the reaction diffusion problem and the bidirectional coupling problem reference solutions are provided in `.mat` files. This implementation focuses on the case where the inner ODE solvers used are explicit Runge-Kutta methods of the same order of convergence as the MERK method. Other explicit Runge-Kutta solvers from butcher.m can be used and with modifications to the MERK solvers implicit methods can be used as well. 
