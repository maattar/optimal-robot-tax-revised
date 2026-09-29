============================================================================

The Calibration Suite

============================================================================

Optimal Taxation of Robotic Capital with a Binding Minimum Wage
                      
M. Aykut Attar
Bilin Neyapti

============================================================================

Instructions:
1) Open << calibrate_parallel.m >> in the script editor.
2) Set SPEC \in {'baseline','equal_shares','auto33','acemoglu2020'}.
3) Run << calibrate_parallel.m >> from the command window. 

Notes: 
1) Install the Global Optimization and Parallel Computing toolboxes.
2) The calibration results used by the "Policy Suite" are available in the 
.../calibration-suite/ folder. These are the four *.mat files.
3) For any SPEC, a new run may converge to a different solution because 
the refinement stage (simulated annealing) is stochastic. Edit the code 
accordingly to save your optimum with a different name for each SPEC. 
4) Please see Online Appendix D for details.
5) Please send your comments or corrections to << maattar@gmail.com >>. 
