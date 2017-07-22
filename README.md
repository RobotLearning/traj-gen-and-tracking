======== INSTRUCTIONS  ========

Video containing example trajectories can be found in
video/ex_trajectories.mp4

To run a simulation containing solo trials with a ballgun, 
please do the following:

1. Run tableTennisPractice.m
2. OPT structure includes the following flags:

    DRAW - draw simulation of table tennis if true.
    PLAN - includes algorithms for generating trajectories and their options.
    LOOKUP - use lookup table instead of generating trajectories online.
    TRAIN - train a lookup table using successfully returned balls.
    VHP - use the Virtual Hitting Plane based method if true.
    RECORD - record the simulation inside MATLAB.
    DISTR - initial ball state related options
    CAMERA - std of the observation noise from the 'cameras'.

4. tt is a Table Tennis class that can be used for table tennis practice:

    tt.practice(q0,N)
    runs the solo table tennis practice N times from the initial robot posture q0.

Okan Koc, Jan Peters, Guilherme Maeda, Online optimal trajectory generation
for robot table tennis, 2017.
