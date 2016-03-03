Robot experiments (MATLAB/SL/REAL):

- Test SL connection and correct frequency
- Implementing NLOPT in C ? [what to do with Kalman filter?]
- Can we incorporate table constraints?

ILC/control notes:
- How does aILC perform based on updated Kalman filter code
- test performance of ILC with linearization at each iteration
- Correct initial conditions of zero-phase filtering and test in MATLAB with ILC 
  [filtfilt is available in SL]
- Why does ILC learning with feedback generate zero Finv?
- Test LQG and iLQR, DDP on MATLAB
- Implement REPS, PI2 algorithms on MATLAB, variational Bayes
- Articulate inverse dynamics in MATLAB does not match to SL! [in SL it is as good as NE]
implement the test function in MATLAB that checks for differences
- Why are two tracking LQRs not exactly the same (at the end)? 
  [maybe R dependence is not correct, index difference?]

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
