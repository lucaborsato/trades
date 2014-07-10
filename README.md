trades
======

Fortran90 code for the orbital simulation of exoplanetary system with simultaneous fit of Transit Times and Radial Velocities.    
It searches for the best parameter configuration by:    

* Levenberg-Marquardt (local minimun) algorithm by MINPACK    
* Grid search (one body, parameters: mass, period/semi-major axis, eccentricity, and arument of pericenter)    
* Genetic Algorithm (PIKAIA)    
* Particle Swarm Optimization (PSO)    

Still in development :-)
