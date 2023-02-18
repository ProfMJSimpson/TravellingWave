# TravellingWave
These algorithms provide a range of numerical, exact and approximate solutions of the nonlinear reaction-diffusion model described in Hajek, McCue and Simpson (2023).  

Hajek_PDE_Julia.jl computes a numerical solution of the reaction-diffusion PDE using the method of lines, and compares the solution with the exact and approximate solutions.
This code is fully commented.

Hajek_PDE_Julia_LimitingCasea0.jl computes a numerical solution of the reaction-diffusion PDE using the method of lines in the limit $a=0$.  This numerical solution is
compared with the exact solution where $a=0.01$.  This second code is not thoroughly commented, but the structure of the code is identical to the structure in the previous 
code, so comments in that code apply here.
