TODO List:

- Why are two tracking LQRs not exactly the same (at the end)? 
  [maybe R dependence is not correct, index difference?]
- Effects of error coupling on DMPs?
- How to avoid large accelerations/large control inputs at start/end
  of trajectory?

- Test filtering functions in SL. How to make ILC work with filtering?
- Make MATLAB experiments for generalizing ILC with RR
- Correct DMP wILC with Barrett arm
- Is LQR different on SL for different trajectories? Can we construct LQR matrix K online?
- Adapt to Yanlong's table tennis scenario in SL
- Test 200 Hz learning on MATLAB and SL
- Correct paper, put simulation results + comparison with REPS,PI2
- Load weights of DMP in SL instead of traj
- Code DMP and wILC in SL
- Make MATLAB experiments for generalizing ILC with Barrett arm

- Symplectic Euler causes problems for convergence
- Why does the height of the basis functions matter?
- Reference DMP should start with zero velocity
- Why do we have to run kinematics for generalizing ILC, i.e. [~,refNew] = rr.kinematics(jNew);
- How can we learn with a small number of weights?