function [data, data_aux] = load_bioRobData_incremental( filename, data, prm)
% type_u_ilc: defines what command is going to be predicted. When the code is
%             running properly it should be 'ilc' type. But to debug the code you
%             should use 'u'.
% type_q_qdot: decide if you should use the position or the velocity as the
%               variable q
%

    dt = prm.data.dt;
    type_u_ilc  = prm.data.u_ilc;
    type_q_qdot = prm.data.q_qdot;


    % load the file here
    path_data = ['../common_data/' filename];
    load_tmp = load(path_data);
    
    [data, data_aux] = format_data(load_tmp, dt, data, type_q_qdot, type_u_ilc);
    
    
    

end


















