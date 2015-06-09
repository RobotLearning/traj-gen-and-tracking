Theoretical TODO:

- Derive Kalman filter from minimum variance and derive seperation principle.
- Read some more ILC papers
- Read robotics book up to control chapters
- Read policy search review
- Read the Barrett WAM inertial specifications
- Read maximum principle chapter

Matlab List:

- check power implementation 
- check total least squares on the ball prediction dataset
- Why does ILC learning with feedback generate zero Finv?
- Merge Guilherme's PMP idea with this folder.

Paper List:

- Compare different mILC with aILC, bILC and only feedback-based ILC (current iteration)

SL Notes:

- Can we run traj again after f and st without problems?
- Why does task servo feedback blow up in Robot computer?
- Check dmp_task. Why are there oscillations? 
- Clean-up: Check out Yanlong's and Katharina's code.

Real Robot Notes:

- Check the initialization pos/vel/acc errors with PID/PD and different filter values.
- Test with f and goto_posture.
- Test bILC with p = 0.1. Test mILC with line search with beta < 1?

Minor:
- Correct initial conditions of zero-phase filtering.
- Load constraints umin and umax, joint limits from sensor offsets file and include in lift_constraints
- Try IDM as fb in MATLAB
- Test LQG on MATLAB
- Implement REPS, PI2 algorithms on MATLAB
- Copy kinematics (forward) from SL
- Articulate inverse dynamics in MATLAB does not match to SL! [in SL it is as good as NE]
implement the test function in MATLAB that checks for differences
- Why are two tracking LQRs not exactly the same (at the end)? 
  [maybe R dependence is not correct, index difference?]
- Check RDMP

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