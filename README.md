======== INSTRUCTIONS FOR IROS 2016 ========

A New Trajectory Generation Framework in Robotic Table Tennis

1. Run tableTennisPractice.m
2. OPT structure encodes the following options:

DRAW - draw simulation of table tennis if true.
LOOKUP - use lookup table instead of generating trajectories online.
TRAIN - train a lookup table using successfully returned balls.
VHP - use the Virtual Hitting Plane based method if true.
RECORD - record the simulation inside MATLAB.

3. STD structure encodes the following options:

pos - std of the initial state position
vel - std of the initial state velocity
camera - std of the observation noise from the 'cameras'.

4. tt is a Table Tennis class that can be used for table tennis practice:

tt.practice(q0,N)

runs the solo table tennis practice N times from the initial robot posture q0.

Okan Koc, Jan Peters, Guilherme Maeda, IROS 2016 Submission, March 2016.