TODO List:

General Issues:

1. Compute SVD faster for the special structure of F without even forming F
2. Are the weights necessary for generalization? 
3. Test filtering functions in SL.
4. How to make ILC work with noisy errors?
5. How does aILC work perform on the robot classes?
6. Investigate LQR differences for different trajectories. 
7. Adapt ILC fully to Yanlong's table tennis task. Check the different functions in SL for table tennis.
8. Correct paper with new results, change methodology
9. Implement EM and then REPS, PI2 algorithms on MATLAB
10. Code DMP in SL, check ddmp code of Jens.
11. Try simple ILC on the robot with a whole trajectory
12. Fast ways to construct, parameterize LQR matrix K or ILC matrix F online?
13. Test end-point learning ILC
14. Make MATLAB experiments for generalization: DDP, regression, convex hull learning. 

Technical Issues:

- Symplectic Euler causes problems for convergence
- Why does the width (h) of the basis functions matter?
- Reference DMP (e.g. loaded dmp.txt) should start with zero velocity
- How can we learn with a smaller (e.g. 20) number of weights?
- Why are two tracking LQRs not exactly the same (at the end)? 
  [maybe R dependence is not correct, index difference?]
- Effects of error coupling on DMPs?
- Avoid large accelerations/large control inputs at start/end
  of DMP?
- How to apply downsampling with ILC to ensure convergence?
- Why does ILC learning with feedback not improve?