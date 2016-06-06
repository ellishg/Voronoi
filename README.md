# Voronoi

Creates a voronoi diagram on the surface of a sphere. It uses Fortune's algorithm for speed.

Known bugs:

Usually towards the bottom of the sphere there is a missing edge. I suspect that when the beach has only two arcs, they do not get connected. 
Everything usually looks fine for 100 or fewer sites, but with more sites there seems to be many extra edges. This may happen for a few special cases. I do not handle the possibility that two sites could occupy the same location. I also do not handle the possibility that three arcs intersect in the same location. 

Here are a few optimizations that I could possibly do:
The beachline is implimented using a circularly doubly linked skip list, but I do not take full advantage of this when scanning the beachline. 
There are tons of calculations in the make_circle, parabolic_intersection, and phi_to_point functions. Maybe I could find some clever way to speed these up.
