TODO List:

General Issues:
- Correct DMP wILC with Barrett arm
- Test filtering functions in SL. How to make ILC work with noisy errors? Check aILC?
- Is LQR different on SL for different trajectories? Can we construct LQR matrix K online?
- Adapt to Yanlong's table tennis scenario in SL
- Possible to reduce computation time by scaling down frequency of ILC/DMP? [e.g. to 200 Hz]
- Correct paper, put simulation results + comparison with REPS,PI2
- Load weights of DMP in SL instead of traj
- Code DMP and wILC in SL
- Make MATLAB experiments for generalizing ILC with Barrett arm

Technical Issues:
- Symplectic Euler causes problems for convergence
- Why does the width (h) of the basis functions matter?
- Reference DMP (e.g. loaded dmp.txt) should start with zero velocity
- Why do we have to run kinematics for generalizing ILC, i.e. [~,refNew] = rr.kinematics(jNew);
- How can we learn with a smaller (e.g. 20) number of weights?
- Why are two tracking LQRs not exactly the same (at the end)? 
  [maybe R dependence is not correct, index difference?]
- Effects of error coupling on DMPs?
- Avoid large accelerations/large control inputs at start/end
  of DMP?

reduce frequency of trajectories:

-load dmp
-downsample to 50 hz
-construct F in 50 hz
-traj.t should be in 50 hz
-update u's in 50 hz
-upsample to 500 hz

-dont change feedback freq