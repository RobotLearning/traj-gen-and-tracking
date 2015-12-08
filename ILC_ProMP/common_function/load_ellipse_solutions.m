function ilc_sol = load_ellipse_solutions(data_folder)

    ilc_sol=[];

    files ={ 'ilcsol_eA_30.mat',  
             'ilcsol_eB_30.mat',
             'ilcsol_eC_30.mat',
             'ilcsol_eD_30.mat',
             'ilcsol_eE_30.mat' };
    
    
    nSol = numel(files);
    for k=1:nSol
        ilc_sol{k} = load([data_folder files{k}]);
        ilc_sol{k}.loaded_file = files{k};
    end
         
         
         
end