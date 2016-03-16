#  Lab 06 - Higher Order Mappings
## Deal.II Users and Developers Training 
### SMR2909 - MHPC P2.5

**Timo Heister** <heister@clemson.edu> 
and
**Luca Heltai** <luca.heltai@sissa.it>

* * * * *


1.  The topic of this lab session is a modified version of step-4 made
    available for you
    <https://www.dealii.org/8.4.0/doxygen/deal.II/step_4.html>

2.  For more information in higher order mappings see step-10\
    <https://www.dealii.org/8.4.0/doxygen/deal.II/step_10.html>

3.  Run the program and check the graphical and text output.

4.  Adjust the right-hand side and solution to get the manufactured
    solution $$u(x,y) = r^2 \sin(3\theta)\cos(\frac{1}{2} \pi r)$$ and
    apply zero boundary conditions. You can use wolframalpha.com to
    compute $- \triangle u$. Make sure the L2 errors are converging.

5.  What mapping degree is required to get optimal convergence of the L2
    error based on the polynomial degree of the finite element space?

6.  Try getting H1 convergence to work correctly too.


