 

The main idea of this tutorial is to use the trajectories explored by ILC
when learning how to reach positions rA, rB, rC and rD to quickly 
converge to the reaching of rE.
The basic idea is to learn correlations between trajectories and motor
commands in the form of Probabilistic Motor Primitives (promp).

Because there are several steps, the code was broken in the following
scripts

    step0_plot_save_all_refs.m
    step1_train_from_zero.m
    step2_conditioned_guess.m
    step3_ilc_conditioned_guess.m

Step0:
Run this script to create and store different reaching tasks.
./main_reaching_p2p/step0_plot_save_all_refs.m
There are five different reaching positions, named rA, rB, rC, rD, rE.

Step1:
For each rA, rB, rC, rD, rE the script
./main_reaching_p2p/step1_train_from_zero.m
will find the optimal torques of the bioRob via ILC. 
This script uses a naive approach where each task is learned independent of each other. We will store their solutions as 
  ilcsol_rA_30.mat
  ilcsol_rB_30.mat
  ilcsol_rC_30.mat
  ilcsol_rD_30.mat
  ilcsol_rE_30.mat
"30" means the number of ILC updates.
The 30 different trajectories of each task are then used to create a ProMP in the next step.

Step2:
Given the ILC solutions computed in step2, the script
./main_reaching_p2p/step2_conditioned_guess.m
will create ProMPs and condition each of them on the trajectory to reach position rE, such that the prediction of the torques can be used as the initial guess for ILC to try reach rE without much iterations.

Step3:
The script
./main_reaching_p2p/step3_ilc_conditioned_guess.m
will load the initial guesses from the previous step, and initialize ILC to reach rE.
