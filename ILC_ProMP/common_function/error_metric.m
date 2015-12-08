function rmse_metric    = error_metric( ref, iter )


    %dsquare = ( ref.xyz(:,1) - iter.xyz(:,1) ).^2 + ( ref.xyz(:,3) - iter.xyz(:,3) ).^2;
    
    dsquare = ( ref(:,1) - iter(:,1) ).^2 + ( ref(:,3) - iter(:,3) ).^2;

    rmse_metric = sqrt(  sum(dsquare) /size(dsquare,1)  );
    

end