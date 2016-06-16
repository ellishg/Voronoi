# Voronoi

Creates a voronoi diagram on the surface of a sphere. It uses Fortune's algorithm for speed.

Known bugs:

Bad seeds as shown in the main.cpp.



Here are a few optimizations that I could possibly do:

There are tons of calculations in the make_circle, parabolic_intersection, and phi_to_point functions. Maybe I could find some clever way to speed these up.