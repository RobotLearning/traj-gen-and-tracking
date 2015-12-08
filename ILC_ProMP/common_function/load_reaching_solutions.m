function ilc_sol = load_reaching_solutions(data_folder)

    ilc_sol=[];
    
    % load the independently learned ILC solutions for each of the reaching
    % positions.
    files ={ 'ilcsol_rA_30.mat',  
             'ilcsol_rB_30.mat',
             'ilcsol_rC_30.mat',
             'ilcsol_rD_30.mat',
             'ilcsol_rE_30.mat' };    
    
    nSol = numel(files);
    for k=1:nSol
        ilc_sol{k} = load([data_folder files{k}]);
        ilc_sol{k}.loaded_file = files{k};
    end
         
         
         
end