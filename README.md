===== INSTRUCTIONS FOR IROS 2016 =====

A New Trajectory Generation Framework in Robotic Table Tennis

1. Run tableTennisPractice.m
2. OPT structure encodes the following options:

i) DRAW - draw simulation of table tennis if true.
ii) LOOKUP - use lookup table instead of generating trajectories online.
iii) TRAIN - train a lookup table using successfully returned balls.
iv) VHP - use the Virtual Hitting Plane based method if true.
v) RECORD - record the simulation inside MATLAB.

3. STD structure encodes the following options:

i) pos - std of the initial state position
ii) vel - std of the initial state velocity
iii) camera - std of the observation noise from the 'cameras'.

4. tt is a Table Tennis class that can be used for table tennis practice:

tt.practice(q0,N)

runs the solo table tennis practice N times from the initial robot posture q0.

Okan Koc, Jan Peters, Guilherme Maeda, IROS 2016 Submission, March 2016.