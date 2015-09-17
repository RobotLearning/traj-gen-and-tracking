ROADMAP

- Reading ILC review + papers
- Read control chapter of both books
- Look at last ACC papers on ILC
- Clean-up: Check out Yanlong's and Katharina's code.
- Read system identification survey
- Derive Kalman filter from minimum variance and derive seperation principle.
- SL code debugging framework in Eclipse
- ATLAS/BLAS library integration with SL (wrapper?)
- Write down robust trajectory generation framework for table tennis
- (Guilherme) Extrapolating ILC. Boosting ILC with weak learners?

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
- Total Gaussian Process?
- Relaxation Learning Control: how to define cost function? Semi markov decision process? Optimal stopping?
  relation to path following? Metronome? Apply monotonic convergence criteria and/or repetitive systems theory
  (Rogers and Owens,1992)?
- Kernel LQR for robot control? Related to representer theorem and/or scenario approach? 
  constructing a lyapunov function based on data (that upper bounds P?)
- Effects of error-inverse-covariance estimation on Q matrix and ILC convergence?
- Transversality conditions for a dmp? Linearization around a dmp?
- Using quasi-Newton and system identification / (recursive) least squares idea.
- Can we compute one (stable) K for all smooth robot trajectories? 
  Investigate LQR differences for different trajectories.
- Can we couple the DMPs
- Evolving DMPs based on current position, as an oracle. Effects of error coupling on DMPs?
- Minimizing expectation with variance added for ILC to come up with a new update rule
- Check DDP and explore fully the connection with ILC (and regression?)
- Recursive pseudoinverse feasible? Connection to IDM and to adjoint system? 
  Perturbation of lower triangular block Toeplitz matrix structure with nonlinear models.
- Prove ILC convergence : keep Fk bounded and show that the cost fnc is convex 
  (relation to convexity of \lambda_max?)
- Fast ways to parameterize LQR matrix K or ILC matrix F for different trajectories?
- Riemannian statistic as a way to estimate f(x,u) based on learning different Fs (linearization)?
- Chaos control to induce bifurcation to a stable orbit during learning? Does it apply only to human 
  motor control?
- Reduced rank approximation with Dirichlet process on Gaussian Process
- GPUCB to estimate the value function instead of stage costs?