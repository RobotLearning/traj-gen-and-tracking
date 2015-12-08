function [extraFolder1] = add_required_folders_promp()

    
    
    % defining paths
    thisPath = pwd;
    % add some necesary stuff into the path
    fprintf('Loading part of Peter Corkes toolbox... \n\n');

    if 0 % make code self-contained
        try
            addpath('~/Dropbox/myMatlabFunctions/');
        end
        try
            addpath('/localdata/Dropbox/myMatlabFunctions/');
        end    
    end
    cd ..

    % adding common stuff: PC toolbox
    addpath([pwd '/common_PCorkeRobotics/']);
    addpath([pwd '/common_PCorkeRobotics/robot/']);
    addpath([pwd '/common_PCorkeRobotics/common/']);
    
    % adding common stuff: ILC and bioRob
    addpath([pwd '/common_function/']);
    addpath([pwd '/common_function/trajectory/']);
    addpath([pwd '/common_function_plot']);
    
    % adding the propm folder here so that the main function can be called
    addpath([pwd '/func_promp']);
    extraFolder1 = [pwd '/func_promp'];

    addpath([pwd '/func_aux/']);
    
    cd(thisPath);

end