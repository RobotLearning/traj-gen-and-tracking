TODO List:

Major TODO:

- Make MATLAB experiments for generalization: regression, convex hull learning. 
  Are the weights necessary for generalization? 
- Fix indexing issue on uffs vs. states
- How to filter correctly (test filter fncs on MATLAB)
- Try IDM as fb in MATLAB
- Check the different functions in SL for table tennis.
- Correct paper with new results, change methodology

- Check DDP and explore fully the connection between ILC and regression
- Recursive pseudoinverse feasible? Connection to IDM?
- Prove ILC convergence : keep Fk bounded and show that the cost fnc is convex
- Optimize DMPs with minimum jerk criterion
- Incorporate extending horizon as worst case for experiments
- Implement REPS, PI2 algorithms on MATLAB


Minor TODO:

- Articulate inverse dynamics in MATLAB does not match to SL! [in SL it is as good as NE]
implement the test function in MATLAB that checks for differences
- How to take inverse in end-goal learning in Mayer form
- Why does previous SL version show different results? [bringing back to q0]
- Symplectic Euler causes problems for convergence in RR
- Why does the width (h) of the basis functions matter?
- Why are two tracking LQRs not exactly the same (at the end)? 
  [maybe R dependence is not correct, index difference?]
- Effects of error coupling on DMPs?
- Why does ILC learning with feedback not improve?
- How to compute recursive pseudoinverse of F from A,B,Ainv,Binv fast? Can we use inverse dynamics lin.?
- Does aILC work on robot classes?
- Fast ways to construct, parameterize LQR matrix K or ILC matrix F online?
- Investigate LQR differences for different trajectories.
- Variational Bayes for estimating noise on positions and velocities?