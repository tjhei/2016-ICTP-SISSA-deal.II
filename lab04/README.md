#  Lab 04 - C++ Templates and Dimension Independence (step-4)
## Deal.II Users and Developers Training 
### SMR2909 - MHPC P2.5

**Timo Heister** <heister@clemson.edu> 
and
**Luca Heltai** <luca.heltai@sissa.it>

* * * * *

1.  See documentation of step-4 at
    <https://www.dealii.org/8.4.0/doxygen/deal.II/step_4.html>

2.  Write a member function `void mesh_info()` that prints the following
    information about the triangulation to the screen: 1) number of
    active cells, 2) number of active/used vertices, lines, quads, hexs
    (only if appropriate for the dimension!).

3.  Use `VectorTools::compute_mean_value` (see step-3) and verify the
    convergence order of the mean in 2d and 3d.

4.  Go back to step-1 and visualize the surface of the Torus by creating
    a `Triangulation<2,3>` (2=dimension of the cells, 3=dimension of the
    space) using `GridGenerator::torus` if you haven't done so.

5.  Try to use the function GridTools::rotate in make_grid() to rotate the
    mesh by 45 degrees (only in 2d!). Note that the function doesn't exist in
    3d (hint: function specialization).

6.  Change the mesh to an L-shape, only apply boundary values to the
    faces adjacent to the center (see `set_boundary_indicator()` in the
    step-3 description), change the boundary values to be $1+\|x\|_2$
    and the right-hand side to be 1. Finally, visualize your solutions
    in ParaView in 2d and 3d.

