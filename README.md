# tower_gpmsa

### Main Operation
`runmcmc.m` is main code to run analysis, calls `sc.m` to read input data and then other code in "Code" folder to run MCMC routine.

`design` - LHS design table for simulation input (X and then Theta(s))
  
 `sim_outputs` - Output table corresponding to LHS table inputs to model
 
 `obs_outputs` - Inputs (Xs) and outputs of experimental testing
 
 Theta values (normalized 0-1) are exported in `theta` after calculation
