#  Lab 05 - Error Computation
## Deal.II Users and Developers Training 
### SMR2909 - MHPC P2.5

**Timo Heister** <heister@clemson.edu> 
and
**Luca Heltai** <luca.heltai@sissa.it>

* * * * *

1.  The topic of this lab session is a modified version of step-4 made
    available for you
    <https://www.dealii.org/8.4.0/doxygen/deal.II/step_4.html>

2.  For more information about computing errors see step-7 (it is a bit
    more complicated though)
    <https://www.dealii.org/8.4.0/doxygen/deal.II/step_7.html>

3.  Run the program and check the graphical and text output.

4.  Where is the right-hand side defined and where do the boundary
    conditions come from?

5.  Fix the right-hand side and boundary conditions to get the
    manufactured solution $$u(x) = \sin(\pi x )\cdot\cos(\pi y)$$ and
    make sure the L2 errors are converging.

6.  Increase the polynomial degree of the finite element space and check
    the convergence of the L2 error.

7.  Implement the computation of the H1 error. For this you need to
    compute the gradient of the manufactured solution and implement it
    (see commented out code for a start).

8.  Implement a suitable 3d manufactured solution and test the
    convergence.


