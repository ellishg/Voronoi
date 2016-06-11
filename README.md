# Voronoi

Creates a voronoi diagram on the surface of a sphere. It uses Fortune's algorithm for speed.

Known bugs:

Usually towards the bottom of the sphere there is a missing edge. I suspect that when the beach has only two arcs, they do not get connected. 


Everything usually looks fine for 100 or fewer sites, but with more sites there seems to be many extra edges. This may happen for a few special cases. I do not handle the possibility that two sites could occupy the same location. I also do not handle the possibility that three arcs intersect in the same location. But neither of these appear to be the problem. Somehow there are circle events that get created that are completely above the sweepline. Pushing these circle events onto the priority queue anyway usually solves the problem but not all the time. I suspect that some rounding errors come into play. When I switched from floats to doubles some of the problems went away. I think I also have memory leaks. After I have handled a circle event, I should delete it. But for some reason I get errors when I try to delete a invalid circle event. 



Here are a few optimizations that I could possibly do:

The beachline is implimented using a circularly doubly linked skip list, but I do not take full advantage of this when scanning the beachline. There are tons of calculations in the make_circle, parabolic_intersection, and phi_to_point functions. Maybe I could find some clever way to speed these up.
