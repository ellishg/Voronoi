# Voronoi

Creates a voronoi diagram on the surface of a sphere. It uses Fortune's algorithm for speed.

Known bugs:

Bad seeds as shown in the main.cpp.

For 8000 sites one might find a bad seed roughly 1% of the time. Bad seeds are more common for 16000 sites. There is almost never an error for 4000 sites.

This is most definitely due to the fact that the calculations are very sensitive to rounding errors when the sites are so close together.


This algorithm can be run on one, two, or four threads.


Here are a few optimizations that I could possibly do:

There are tons of calculations in the make_circle, parabolic_intersection, and phi_to_point functions. Maybe I could find some clever way to speed these up.