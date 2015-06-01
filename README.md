Paper List:

- Include new update laws: Mayer form and total least squares.
- Compare mILC with different parameters, aILC, bILC and only feedback-based ILC (current iteration)
- Regression on different trajectories?
- Include wILC ?

Real Robot Notes:

- We can test the initialization problem with polynomials.
- Test with f and goto_posture.
- Test bILC with p = 0.1. Test mILC with line search with beta < 1?

SL Notes:

- Can we run traj again after f and st without problems?
- Why does task servo feedback blow up in Robot computer?
- Check dmp_task. Why are there oscillations? 
- Clean-up: Check out Yanlong's and Katharina's code.

Matlab Notes:

- How does the mayerFormILC perform with/out Kalman filter with/out Convex Opt. with/out TLS?
- Does wILC work better with the new dmp?
- Make MATLAB experiments for generalization: regression, convex hull learning. 
  Are the weights necessary for generalization? 
- Merge Guilherme's PMP idea with this folder.

- Correct initial conditions of zero-phase filtering.
- Load constraints umin and umax, joint limits from sensor offsets file and include in lift_constraints
- Symplectic Euler causes problems for convergence in RR
- Try IDM as fb in MATLAB
- Test LQG on MATLAB
- Implement REPS, PI2 algorithms on MATLAB
- Copy kinematics (forward) from SL
- Articulate inverse dynamics in MATLAB does not match to SL! [in SL it is as good as NE]
implement the test function in MATLAB that checks for differences
- Why are two tracking LQRs not exactly the same (at the end)? 
  [maybe R dependence is not correct, index difference?]
- Why does ILC learning with feedback not improve?
- Check wILC and RDMP

Theoretical TODO:

- Derive Kalman filter from minimum variance.
- Read some more ILC papers
- Read robotics book up to control chapters
- Read policy search review
- Read the Barrett WAM inertial specifications
- Read maximum principle chapter

Extensions:
- Can we compute one K for all smooth robot trajectories?
- Optimize DMPs with minimum jerk criterion
- Evolving DMPs based on current position, as an oracle. Effects of error coupling on DMPs?
- Minimizing expectation with variance added for ILC to come up with a new update rule
- Check DDP and explore fully the connection with ILC (and regression?)
- Recursive pseudoinverse feasible? Connection to IDM and to adjoint system?
- Prove ILC convergence : keep Fk bounded and show that the cost fnc is convex
- Fast ways to parameterize LQR matrix K or ILC matrix F for different trajectories?
- Investigate LQR differences for different trajectories.
- Variational Bayes for estimating noise on positions and velocities?