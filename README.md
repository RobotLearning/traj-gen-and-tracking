TODO List:

Major TODO:

- Try simple ILC on the robot with a whole trajectory
- Are the weights necessary for generalization? 
- Check ddmp code and other functions in SL
- How to make ILC work with noisy errors? 
- Adapt ILC fully to Yanlong's table tennis task. Check the different functions in SL for table tennis.
- Correct paper with new results, change methodology
- Implement REPS, PI2 algorithms on MATLAB
- Test end-point learning ILC
- Make MATLAB experiments for generalization: DDP, regression, convex hull learning. 

Minor TODO:

- Symplectic Euler causes problems for convergence in RR
- Why does the width (h) of the basis functions matter?
- Reference DMP (e.g. loaded dmp.txt) should start with zero velocity
- How can we learn with a smaller (e.g. 20) number of weights?
- Why are two tracking LQRs not exactly the same (at the end)? 
  [maybe R dependence is not correct, index difference?]
- Effects of error coupling on DMPs?
- Avoid large accelerations/large control inputs at start/end
  of DMP?
- Why does ILC learning with feedback not improve?
- How to compute recursive pseudoinverse of F from A,B,Ainv,Binv fast? Can we use inverse dynamics lin.?
- Does aILC work on robot classes?
- Fast ways to construct, parameterize LQR matrix K or ILC matrix F online?
- Investigate LQR differences for different trajectories.
- Why does only dmpILC show improvement on generalization?
- Test inverse dynamics in feedback mode
- Fix the indexing problem when up/downsampling with unom,K,F