%% Load table stand values

S1 = [table_center + stand_x - stand_radius;
      dist_to_table - table_y + stand_y + stand_radius; 
      table_z - table_thickness];
S2 = [table_center + stand_x + stand_radius;
      dist_to_table - table_y + stand_y + stand_radius; 
      table_z - table_thickness];
S3 = [table_center + stand_x + stand_radius;
      dist_to_table - table_y + stand_y - stand_radius; 
      table_z - table_thickness];
S4 = [table_center + stand_x - stand_radius;
      dist_to_table - table_y + stand_y - stand_radius; 
      table_z - table_thickness];
S5 = [table_center + stand_x - stand_radius;
      dist_to_table - table_y + stand_y + stand_radius; 
      floor_level];
S6 = [table_center + stand_x + stand_radius;
      dist_to_table - table_y + stand_y + stand_radius; 
      floor_level];
S7 = [table_center + stand_x + stand_radius;
      dist_to_table - table_y + stand_y - stand_radius; 
      floor_level];
S8 = [table_center + stand_x - stand_radius;
      dist_to_table - table_y + stand_y - stand_radius; 
      floor_level];
SA = [S1,S2,S3,S4,S5,S6,S7,S8]';

S1 = [table_center - stand_x - stand_radius;
      dist_to_table - table_y + stand_y + stand_radius; 
      table_z - table_thickness];
S2 = [table_center - stand_x + stand_radius;
      dist_to_table - table_y + stand_y + stand_radius; 
      table_z - table_thickness];
S3 = [table_center - stand_x + stand_radius;
      dist_to_table - table_y + stand_y - stand_radius; 
      table_z - table_thickness];
S4 = [table_center - stand_x - stand_radius;
      dist_to_table - table_y + stand_y - stand_radius; 
      table_z - table_thickness];
S5 = [table_center - stand_x - stand_radius;
      dist_to_table - table_y + stand_y + stand_radius; 
      floor_level];
S6 = [table_center - stand_x + stand_radius;
      dist_to_table - table_y + stand_y + stand_radius; 
      floor_level];
S7 = [table_center - stand_x + stand_radius;
      dist_to_table - table_y + stand_y - stand_radius; 
      floor_level];
S8 = [table_center - stand_x - stand_radius;
      dist_to_table - table_y + stand_y - stand_radius; 
      floor_level];
SB = [S1,S2,S3,S4,S5,S6,S7,S8]';

S1 = [table_center - stand_x - stand_radius;
      dist_to_table - table_y - stand_y + stand_radius; 
      table_z - table_thickness];
S2 = [table_center - stand_x + stand_radius;
      dist_to_table - table_y - stand_y + stand_radius; 
      table_z - table_thickness];
S3 = [table_center - stand_x + stand_radius;
      dist_to_table - table_y - stand_y - stand_radius; 
      table_z - table_thickness];
S4 = [table_center - stand_x - stand_radius;
      dist_to_table - table_y - stand_y - stand_radius; 
      table_z - table_thickness];
S5 = [table_center - stand_x - stand_radius;
      dist_to_table - table_y - stand_y + stand_radius; 
      floor_level];
S6 = [table_center - stand_x + stand_radius;
      dist_to_table - table_y - stand_y + stand_radius; 
      floor_level];
S7 = [table_center - stand_x + stand_radius;
      dist_to_table - table_y - stand_y - stand_radius; 
      floor_level];
S8 = [table_center - stand_x - stand_radius;
      dist_to_table - table_y - stand_y - stand_radius; 
      floor_level];
SC = [S1,S2,S3,S4,S5,S6,S7,S8]';

S1 = [table_center + stand_x - stand_radius;
      dist_to_table - table_y - stand_y + stand_radius; 
      table_z - table_thickness];
S2 = [table_center + stand_x + stand_radius;
      dist_to_table - table_y - stand_y + stand_radius; 
      table_z - table_thickness];
S3 = [table_center + stand_x + stand_radius;
      dist_to_table - table_y - stand_y - stand_radius; 
      table_z - table_thickness];
S4 = [table_center + stand_x - stand_radius;
      dist_to_table - table_y - stand_y - stand_radius; 
      table_z - table_thickness];
S5 = [table_center + stand_x - stand_radius;
      dist_to_table - table_y - stand_y + stand_radius; 
      floor_level];
S6 = [table_center + stand_x + stand_radius;
      dist_to_table - table_y - stand_y + stand_radius; 
      floor_level];
S7 = [table_center + stand_x + stand_radius;
      dist_to_table - table_y - stand_y - stand_radius; 
      floor_level];
S8 = [table_center + stand_x - stand_radius;
      dist_to_table - table_y - stand_y - stand_radius; 
      floor_level];

SD = [S1,S2,S3,S4,S5,S6,S7,S8]';

%% Load robot stand

% TODO:
% robot base is at 0 z-level normally
% but bringing it down to arm level for simplicity for now
% later when adding joint1 graphics this will be modified
robot_base = -0.4;
robot_base_thickness = 0.1;
robot_base_radius = 0.2;

S1 = [robot_base_radius;
      robot_base_radius; 
      robot_base];
S2 = [robot_base_radius;
      -robot_base_radius; 
      robot_base];
S3 = [-robot_base_radius;
      -robot_base_radius;
      robot_base];
S4 = [-robot_base_radius;
      robot_base_radius;
      robot_base];
S5 = [robot_base_radius;
      robot_base_radius; 
      robot_base + robot_base_thickness];
S6 = [robot_base_radius;
      -robot_base_radius;  
      robot_base + robot_base_thickness];
S7 = [-robot_base_radius;
      -robot_base_radius; 
      robot_base + robot_base_thickness];
S8 = [-robot_base_radius;
      robot_base_radius; 
      robot_base + robot_base_thickness];
SR = [S1,S2,S3,S4,S5,S6,S7,S8]';
