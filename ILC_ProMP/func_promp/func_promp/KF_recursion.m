
function [x_new, P_new] = KF_recursion(x_old, P_old, H, z, R_obs)

    
    tmp   =   H*P_old*H' + R_obs ;
    K     =   (P_old*H') / ( tmp  );

    % Covariance matrix is independent of measurement
    P_new = P_old - K*H*P_old;

    % x_old: first guess
    % z: observed values
    x_new = x_old + K*(z - H*x_old);

    
end