TODO List:

- Implement DMP in SL and create dmp_task. Load using weights only
- Go to initial posture slowly with a dmp (zero forcing)
- Try zero-phase filtering in SL as a noncausal offline filter. Test in MATLAB first.
- Review Kalman filter
- Does aILC work on robot classes? EM algorithm as an extension ?
- Read some more ILC papers

REAL ROBOT TODO:

- Test LQR and different feedback (LQG?) w/o learning
- Test different DMPs, with different taus and velocities w/o learning
- Incorporate extending horizon as worst case for experiments

General TODO (i.e. MATLAB):

- Make MATLAB experiments for generalization: regression, convex hull learning. 
  Are the weights necessary for generalization? 
- Optimize DMPs with minimum jerk criterion
- Test LQG on MATLAB 
- Correct paper with SL results, change methodology
- Implement REPS, PI2 algorithms on MATLAB

READING TODO:

- Read robotics book up to control chapters
- Read policy search review
- Read the Barrett WAM inertial specifications
- Read maximum principle chapter

Extensions:

- Evolving DMPs based on current position, as an oracle 
- Minimizing expectation with variance added for ILC to come up with a new update rule
- Total Least Squares implementation for ILC?
- Check DDP and explore fully the connection with ILC (and regression?)
- Recursive pseudoinverse feasible? Connection to IDM?
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