#  Lab 02 - Sparsity Patterns (step-2)
## Deal.II Users and Developers Training 
### SMR2909 - MHPC P2.5

**Timo Heister** <heister@clemson.edu> 
and
**Luca Heltai** <luca.heltai@sissa.it>

* * * * *

1.  See documentation at
    <https://www.dealii.org/8.4.0/doxygen/deal.II/step_2.html>

2.  Copy and run step-2. Look at the sparsity patterns in firefox.

3.  How does the pattern change if you increase the polynomial degree
    from 1 to 2 or to 3?

4.  How does the pattern change if you use a globally refined (say 3
    times) unit square?

5.  Are these patterns symmetric? Why/why not?

6.  How many entries per row in the sparsity pattern do you expect for a
    Q1 element (assuming four cells are around each vertex)? Check that
    this is true for the mesh in b) (look for `row_length(i)` and output
    them for each row). Can you construct a 2d mesh (without hanging
    nodes) that has a row with more entries?

7.  Print all entries for row 42 for the original renumbered sparsity
    pattern.

8.  Bonus: Compute and output statistics like the the number of
    unknowns, bandwidth of the sparsity pattern, average number of
    entries per row, and fill ratio.


