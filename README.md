ROADMAP

- Train Bounce and contact models
- Train contact model based on step 3 and step 2
- Train EKF on ball data with EM to find C and g
- Get viable velocities based on Miki's code
- Test stochastic maximum principle in SL
- Write down draft for Table Tennis Performance Criterion

MATLAB/SL list:
- ATLAS/BLAS library integration with SL (wrapper?)
- How does aILC perform based on updated Kalman filter code
- test performance of ILC with linearization at each iteration
- Why does ILC learning with feedback generate zero Finv?
- Try IDM as fb in MATLAB
- Test LQG and iLQR, DDP on MATLAB
- Implement REPS, PI2 algorithms on MATLAB, variational Bayes
- Articulate inverse dynamics in MATLAB does not match to SL! [in SL it is as good as NE]
implement the test function in MATLAB that checks for differences
- Why are two tracking LQRs not exactly the same (at the end)? 
  [maybe R dependence is not correct, index difference?]

Robot list:
- Experiment with Ballgun fixed.
- Correct initial conditions of zero-phase filtering and test in MATLAB with ILC 
  [filtfilt is available in SL]
- Implement/Learn RDMP first in MATLAB and then in SL

Extensions:
- (Guilherme) Extrapolating ILC. Boosting ILC with weak learners?
- Relaxation Learning Control: how to define cost function? Semi markov decision process? Optimal stopping?
  relation to path following? Metronome? Apply monotonic convergence criteria and/or repetitive systems theory
  (Rogers and Owens,1992)?
- Kernel LQR for robot control? Related to representer theorem and/or scenario approach? 
  constructing a lyapunov function based on data (that upper bounds P?)
- Effects of error-inverse-covariance estimation on Q matrix and ILC convergence?
- Using quasi-Newton and system identification / (recursive) least squares idea.
- Can we compute one (stable) K for all smooth robot trajectories? 
  Investigate LQR differences for different trajectories.
- Transversality conditions for a dmp? Linearization around a dmp? 
  Evolving DMPs based on current position, as an oracle. Effects of error coupling on DMPs? 
  Can we couple the DMPs?
- Minimizing expectation with variance added for ILC to come up with a new update rule
- Check DDP and explore fully the connection with ILC (and regression?)
- Recursive pseudoinverse feasible? Connection to IDM and to adjoint system? 
- Perturbation of lower triangular block Toeplitz matrix structure with nonlinear models.
- Fast ways to parameterize LQR matrix K or ILC matrix F for different trajectories?
- Riemannian statistics as a way to estimate f(x,u) based on learning different Fs (linearization)?
  riemannian singular values?
- Chaos control to induce bifurcation to a stable orbit during learning? Does it apply only to human 
  motor control?
- Reduced rank approximation with Dirichlet process on Gaussian Process
