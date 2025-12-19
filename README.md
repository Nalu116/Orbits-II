# Orbits-II
This is not at all polished lol I just grabbed files from back in the day and threw them here for reference (I was a terrible student) 

Classical_to_Cartesian_Orbital_Element_Converter.m is a function where user may input classical orbital elements to convert to cartesian positions and first rates of change about all three axes.

Cartesian_to_Classical_Orbital_Element_Convertor.m is a script that converts user inputed cartesian values to classical orbital elements

Transfers.m is a script that parses through planetary ephemeris data and plots first order approximations for orbital transfers using lambers TOF equations. Script intakes Julian Start and Julian arrival dates for whatever timeframe mission is desired. This specific script only applies to an Earth/Mars Interplanetary Direct Transfer, using patched conic method about the three spheres of influence (Earth-Centric, Helio-Centric, Mars-Centric) and plots specific energies to accomidate planetary transfers. 


ephem.m and lambert_solver.m are given functions where:

ephem.m pulls and propagates spice data from a .mat file
lambert_solver solves Lambert's Problem using Battin's method where user inputs two position vectors, the estimated time of flight, grafitational parameters, and some boundry values for numerical solver.


