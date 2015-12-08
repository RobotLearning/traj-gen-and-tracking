function [mu_x, mu_xd, Sigma_t, Sigma_t1, Sigma_t_t1, Sigma_td_half] = ...
                     get_distributions_at_time_T1(w_mu, w_cov, basis, dt, timePoint)
% INPUT
%    p: promp
%    timePoint =  index*dt (usually timePoint belongs to [0 1]);

%     w_mu = p.w.mean;
%    
%     
%     % this is the covariance of the conditioned PMP
%     w_cov = p.weight.cov;
    
    % get index instead of normalized time
    timePointIndex = round( timePoint / dt);
    
    
    % This should get the value of all n bases at time t and t+1
    Psi_t   = basis.Gn(timePointIndex,:);
    Psi_td  = basis.Gndot(timePointIndex,:);
    Psi_tdd = basis.Gnddot(timePointIndex,:);

    Psi_t1  = basis.Gn(timePointIndex+1,:);
    Psi_t1d = basis.Gndot(timePointIndex+1,:);


    % Put basis and basis_d together for t and t+1
    Phi_t =  [ Psi_t' ,  Psi_td' ];
    Phi_t1 = [ Psi_t1' , Psi_t1d' ];

    % Put basis_d and basis_dd together for t
    Phi_td =  [ Psi_td' , Psi_tdd' ];

    
    mu_x  = Phi_t'  * w_mu;
    mu_xd = Phi_td' * w_mu;

    Sigma_t       = Phi_t'  * w_cov * Phi_t;  % current
    Sigma_t1      = Phi_t1' * w_cov * Phi_t1; % next
    Sigma_t_t1    = Phi_t'  * w_cov * Phi_t1; % cross
    Sigma_td_half = Phi_td' * w_cov * Phi_t;  % derivative
    
    
end