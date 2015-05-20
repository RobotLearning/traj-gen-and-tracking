TODO List:

- Test inverse dynamics and fb on dmp_task after ensuring dmps work.
- Merge Guilherme's PMP idea with this folder
- Copy kinematics (forward) from SL
- Make MATLAB experiments for generalization: regression, convex hull learning. 
  Are the weights necessary for generalization? 
- Test LQG on MATLAB 
- Correct paper with SL results, change methodology
- Implement REPS, PI2 algorithms on MATLAB
- Check out Yanlong's and Katharina's code
- Derive Kalman filter from minimum variance
- Correct initial conditions of zero-phase filtering.
- Load constraints umin and umax, joint limits from file and include in lift_constraints
- Can we make aILC work on robot classes?

Extensions:

- Optimize DMPs with minimum jerk criterion
- Evolving DMPs based on current position, as an oracle 
- Minimizing expectation with variance added for ILC to come up with a new update rule
- Total Least Squares implementation for ILC?
- Check DDP and explore fully the connection with ILC (and regression?)
- Recursive pseudoinverse feasible? Connection to IDM and to adjoint system?
- Prove ILC convergence : keep Fk bounded and show that the cost fnc is convex
- Fast ways to construct, parameterize LQR matrix K or ILC matrix F online?
- How to take inverse in end-goal learning in Mayer form
- Investigate LQR differences for different trajectories.
- Variational Bayes for estimating noise on positions and velocities?

Technical Issues:
- Try IDM as fb in MATLAB
- Articulate inverse dynamics in MATLAB does not match to SL! [in SL it is as good as NE]
implement the test function in MATLAB that checks for differences
- Symplectic Euler causes problems for convergence in RR
- Why does the width (h) of the basis functions matter?
- Why are two tracking LQRs not exactly the same (at the end)? 
  [maybe R dependence is not correct, index difference?]
- Effects of error coupling on DMPs?
- Why does ILC learning with feedback not improve?
- Check wILC and RDMP

Reading:
- Read some more ILC papers
- Read robotics book up to control chapters
- Read policy search review
- Read the Barrett WAM inertial specifications
- Read maximum principle chapter
