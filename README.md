Humanoids TODO list:

- Plot some dmps in Cart space and include in paper.
- Implement RDMP
- Correct initial conditions of zero-phase filtering.
- Analyze the index problem for dmp trajectories. Why +1?
- Why does task servo feedback blow up in Robot computer?

Notes:
- why doesn't ball and robot trajectory coincide in simulations? check kinematics in MATLAB
- test performance of Fact and tls version on 0.5 sec trajectory
- check total least squares on the ball prediction dataset
- Why does ILC learning with feedback generate zero Finv?
- Merge Guilherme's PMP idea with this folder.
- Clean-up: Check out Yanlong's and Katharina's code.
- Load constraints umin and umax, joint limits from sensor offsets file and include in lift_constraints
- Try IDM as fb in MATLAB
- Test LQG on MATLAB
- Implement REPS, PI2 algorithms on MATLAB
- Articulate inverse dynamics in MATLAB does not match to SL! [in SL it is as good as NE]
implement the test function in MATLAB that checks for differences
- Why are two tracking LQRs not exactly the same (at the end)? 
  [maybe R dependence is not correct, index difference?]


Theoretical:

- Transversality conditions for dmp
- Linearization around a dmp
- Derive Kalman filter from minimum variance and derive seperation principle.
- Read some more ILC papers
- Read robotics book up to control chapters
- Read policy search review
- Read the Barrett WAM inertial specifications
- Read maximum principle chapter
- Read system identification survey


Extensions:
- Using quasi-Newton and system identification/ (recursive) least squares idea.
- Can we compute one K for all smooth robot trajectories?
- Can we couple the DMPs
- Evolving DMPs based on current position, as an oracle. Effects of error coupling on DMPs?
- Minimizing expectation with variance added for ILC to come up with a new update rule
- Check DDP and explore fully the connection with ILC (and regression?)
- Recursive pseudoinverse feasible? Connection to IDM and to adjoint system?
- Prove ILC convergence : keep Fk bounded and show that the cost fnc is convex
- Fast ways to parameterize LQR matrix K or ILC matrix F for different trajectories?
- Investigate LQR differences for different trajectories.
- Variational Bayes for estimating noise on positions and velocities?