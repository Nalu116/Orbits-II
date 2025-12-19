# Orbits-II
This is not at all polished lol I just grabbed files from back in the day and threw them here for reference (I was a terrible student) 

Classical2Cart.m is a function where user may input classical orbital elements to convert to cartesian positions and first rates of change about all three axes




ephem.m and lambert_solver.m are given functions where:

ephem.m pulls and propagates spice data from a .mat file
lambert_solver solves Lambert's Problem using Battin's method where user inputs two position vectors, the estimated time of flight, grafitational parameters, and some boundry values for numerical solver.


