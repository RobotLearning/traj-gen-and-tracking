function ilc_sol = load_reaching_solutions_bootstrap(full_path)

%     ilc_sol=[];
% 
%     files ={ 'ilcsol_rA_30.mat',  
%              'ilcsol_rB_30.mat',
%              'ilcsol_rC_30.mat',
%              'ilcsol_rD_30.mat',
%              'ilcsol_rE_30.mat' };
    
    load([full_path])
    
    for k=1:1
        ilc_sol{k} = final_ilc{k}.final;
        ilc_sol{k}.loaded_file = ['previous trained ref:' final_ilc{k}.ref.name];
    end
%     nSol = numel(files);
%     for k=1:nSol
%         ilc_sol{k} = load([full_path]);
%         ilc_sol{k}.loaded_file = files{k};
%     end
         
         
         
end