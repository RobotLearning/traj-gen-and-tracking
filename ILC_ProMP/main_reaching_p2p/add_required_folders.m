function [] = add_required_folders()


    
    % defining paths
    thisPath = pwd;
    % add some necesary stuff into the path
    fprintf('Loading part of Peter Corkes toolbox... \n\n');

    if 0 % not used
        try
            addpath('~/Dropbox/myMatlabFunctions/');
        end
        try
            addpath('/localdata/Dropbox/myMatlabFunctions/');
        end
    end
    
    %cd ..

    % adding common stuff: PC toolbox
    addpath([pwd '/common_PCorkeRobotics/']);
    addpath([pwd '/common_PCorkeRobotics/robot/']);
    addpath([pwd '/common_PCorkeRobotics/common/']);
    
    % adding common stuff: ILC and bioRob
    addpath([pwd '/common_function/']);
    addpath([pwd '/common_function/trajectory/']);
    addpath([pwd '/common_function_plot']);
    
    addpath([pwd '/func_aux/']);
    
     
    cd(thisPath);

end