ROADMAP

- Rereading 2 RL reviews
- Reading ILC review + paper
- Read control chapter of both books
- SL code debugging framework in Eclipse
- ATLAS/BLAS library integration with SL (wrapper?)
- Look at last ACC papers on ILC
- Start writing ACC paper
- Merge Guilherme's PMP idea with this folder.
- Implementing multithreading with SL
- Develop robust control framework for table tennis
- Read system identification survey
- Clean-up: Check out Yanlong's and Katharina's code.
- Derive Kalman filter from minimum variance and derive seperation principle.

Notes:
- compare kinematics with SL
- test performance of ILC with linearization at each iteration
- why doesn't ball and robot trajectory coincide in simulations? check kinematics in MATLAB
- check total least squares on the ball prediction dataset
- Why does ILC learning with feedback generate zero Finv?
- Load constraints umin and umax, joint limits from sensor offsets file and include in lift_constraints
- Try IDM as fb in MATLAB
- Test LQG and iLQR, DDP on MATLAB
- Implement REPS, PI2 algorithms on MATLAB, variational Bayes
- Articulate inverse dynamics in MATLAB does not match to SL! [in SL it is as good as NE]
implement the test function in MATLAB that checks for differences
- Why are two tracking LQRs not exactly the same (at the end)? 
  [maybe R dependence is not correct, index difference?]

Robot list:
- Test error feedback based weighting on the Q matrix
- Correct initial conditions of zero-phase filtering and test in MATLAB with ILC
- Implement/Learn RDMP first in MATLAB and then in SL
- Wrist problem in robot computer?

Extensions:
- Effects of error-inverse-covariance estimation on Q matrix and ILC convergence?
- Transversality conditions for a dmp? Linearization around a dmp?
- Kernel LQR for robot control? 
- Metronome Learning: how to define cost function? Semi markov decision process? Optimal stopping?
  relation to path following?
- Using quasi-Newton and system identification/ (recursive) least squares idea.
- Can we compute one (stable) K for all smooth robot trajectories? 
  Investigate LQR differences for different trajectories.
- Can we couple the DMPs
- Evolving DMPs based on current position, as an oracle. Effects of error coupling on DMPs?
- Minimizing expectation with variance added for ILC to come up with a new update rule
- Check DDP and explore fully the connection with ILC (and regression?)
- Recursive pseudoinverse feasible? Connection to IDM and to adjoint system?
- Prove ILC convergence : keep Fk bounded and show that the cost fnc is convex 
  (relation to convexity of \lambda_max?)
- Fast ways to parameterize LQR matrix K or ILC matrix F for different trajectories?
- Riemannian statistic as a way to estimate f(x,u) based on learning different Fs (linearization)?
- Chaos control to induce bifurcation to a stable orbit during learning? Does it apply only to human 
  motor control?
- GPUCB to estimate the value function instead of stage costs?

NOT UNDERSTOOD:
- Duality in Dynamic programming
- Relation of Kalman filtering to GP 
- MCMC methods