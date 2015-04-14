TODO List:

SL steps

0. Fix indexing issue on uffs vs. states
3. How to filter correctly (test filter fncs on MATLAB)
4. Test articulate IDM in SL and check again in MATLAB
5. Test ILC in MATLAB with nominal robot parameters
   [what is the resulting LQR like?]
6. Cut down DMP first half ? [use kinematics in SL to visualize]
7. Test LQR for bringing back to q0 in SL
8. How to initialize to q0 with filtering turned on in SL?
9. Incorporate extending horizon as worst case for experiments
10. Try IDM as fb in MATLAB
11. Test adding fb to next iteration one step before
12. Test end-goal learning in SL
13. Check DDP and explore fully the connection between ILC and regression
14. Recursive pseudoinverse feasible? Connection to IDM?
15. Prove ILC convergence : keep Fk bounded and show that the cost fnc is convex
16. Optimize DMPs with minimum jerk criterion

Major TODO:

- Try simple ILC on the robot with a whole trajectory
- Are the weights necessary for generalization? 
- Check ddmp code and other functions in SL
- Adapt ILC fully to Yanlong's table tennis task. Check the different functions in SL for table tennis.
- Correct paper with new results, change methodology
- Implement REPS, PI2 algorithms on MATLAB
- Make MATLAB experiments for generalization: DDP, regression, convex hull learning. 

Minor TODO:

- Symplectic Euler causes problems for convergence in RR
- Why does the width (h) of the basis functions matter?
- Reference DMP (e.g. loaded dmp.txt) should start with zero velocity
- How can we learn with a smaller (e.g. 20) number of weights?
- Why are two tracking LQRs not exactly the same (at the end)? 
  [maybe R dependence is not correct, index difference?]
- Effects of error coupling on DMPs?
- Avoid large accelerations/large control inputs at start/end
  of DMP?
- Why does ILC learning with feedback not improve?
- How to compute recursive pseudoinverse of F from A,B,Ainv,Binv fast? Can we use inverse dynamics lin.?
- Does aILC work on robot classes?
- Fast ways to construct, parameterize LQR matrix K or ILC matrix F online?
- Investigate LQR differences for different trajectories.
- Why does only dmpILC show improvement on generalization?
- Test inverse dynamics in feedback mode
- Fix the indexing problem when up/downsampling with unom,K,F

Note:

- Adapting DMP for goal velocities can be much easier! Safe learning with DMPs?