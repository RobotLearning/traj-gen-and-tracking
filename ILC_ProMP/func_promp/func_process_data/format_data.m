function [data, data_aux] = format_data(load_tmp, dt, data, type_q_qdot, type_u_ilc)


    r_data = load_tmp.iter;
    data_aux = load_tmp.prm;
    

    ref = load_tmp.prm.ref;

    clear load_tmp;
    
    % if no dt is given, use the time of the original reference trajectory
    if isempty(dt)
        time = ref.t;
    else
        time = 0:dt:ref.t(end)-dt;
    end
    dt = time(2)-time(1);
    
    nDemoTrain = size(r_data,2);
    
    nStates = 3+3; % here I am treating the control signal u as the 
                  % 4th, 5th and 6th joints

                  
    if isempty(data) % create data structure here
        for j = 1:nStates
            data(j).q       = [];
            data(j).qdot    = [];
        end
    end
    
    for d = 2:nDemoTrain % always cut the first iteration because the ILC there is zero
        for j=1:nStates
            if j <= 3
                switch type_q_qdot
                    case 'q'
                    resampled_q = interp1(ref.t, r_data(d).q(:,j+1)',   time);
                    case 'qdot'
                    resampled_q = interp1(ref.t, r_data(d).qdot(:,j+1)',   time);
                    otherwise
                    
                end                
            else
                switch type_u_ilc
                    case 'ilc'
                    resampled_q = interp1(ref.t, r_data(d).uilc(:,j+1-3)', time);
                    case 'u'
                    resampled_q = interp1(ref.t, r_data(d).u(:,j+1-3)', time);
                end                
            end
            % augmenting the data.q.
            data(j).q = [data(j).q;   resampled_q];
        end
    end
    
    % compute velocity based on dt
    for j=1:nStates
        
        qdot_ = diff(data(j).q, 1, 2)/dt;
        qdot_ = [qdot_ qdot_(:,end)];
        data(j).qdot = qdot_;
         
    end
    
    if 0
        % check if resample was successfull
        figurew('q');
        for j=1:nStates-3
            disp(j)
            subplot(nStates,1,j); hold on;
            plot(ref.t, r_data(1).q(:,j+1), SBLUEBALLW(10));
            plot(time,  data(j).q(1,:), SREDBALLW(5));
        end
        for j=nStates-2:nStates
            disp(j)
            subplot(nStates,1,j); hold on;
            plot(ref.t, r_data(1).u(:,j+1-3), SBLUEBALLW(10));
            plot(time,  data(j).q(1,:), SREDBALLW(5));
        end    
    end    
    
     
    data_aux.ref = ref;
    data_aux.nTrainData = nDemoTrain;
    data_aux.nTrajSize = length(time);
    
    
    % resample the reference so that it is useful as a test data
    for j=1:6
         
        if j<=3
            % resample q
            refq_r = interp1(ref.t, ref.q(:,j+1)',   time);
            
            switch type_q_qdot
                case 'q'
                    refq_r = refq_r; % do nothing
                case 'qdot'
                    % use qdot as q, but just understand them as "states"
                    temp   = diff(refq_r)/dt;
                    refq_r = [temp temp(end)];
                otherwise
            end              
            % find qdot
            refqdot_r = diff(refq_r)/dt;
            refqdot_r = [refqdot_r refqdot_r(end)];
        else
            % resample u
            refq_r = interp1(ref.t, ref.u(:,j+1-3)',   time);

            % find udot (udot is not useful but keep for consistency)
            refqdot_r = diff(refq_r)/dt;
            refqdot_r = [refqdot_r refqdot_r(end)];            
        end
        % fill this additional structure with resampled values for the
        % reference
        ref_resampled(j).t    = time ;
        ref_resampled(j).q    = refq_r ;
        ref_resampled(j).qdot = refqdot_r ;
    end

    data_aux.ref_resampled = ref_resampled;
    
    
end