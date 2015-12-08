function [basis_n, basisD_n, basisDD_n] = generateGaussianBasis( phase, mu, sigma)

    dbg = 0;

    basisCenter = mu;

    z   = phase.z;
    zd  = phase.zd;
    zdd = phase.zdd;


    %> (z - basisCenter)
    z_minus_center = bsxfun(@minus, z, basisCenter );
    
    if dbg % original
    % ====================================
    % original implementation
    % ====================================

        at = bsxfun(@times, z_minus_center, 1./sigma);

        % computing and normalizing basis (order 0)
        basis     = bsxfun(@times, exp( -0.5 * at.^2 ), 1./sigma/sqrt(2*pi) );
        basis_sum = sum(basis,2);
        basis_n   = bsxfun(@times, basis, 1 ./ basis_sum);
    end
        
    % ====================================
    % my simple implementation
    % ====================================
    clear basiss
    for k = 1:numel(mu)
        basiss(k,:) = pdf('Normal', z, mu(k), sigma(k));
    end
    basiss = basiss';

    normSum = sum(basiss,2); 
    if 1 % normalize
        basiss_n   = bsxfun(@times, basiss, 1./ normSum);
    else % skip normalization
        basiss_n   = basiss;
    end

    if dbg
        debugBasis(z, basis_n, basiss_n);
    end
    
    % ====================================
    % use my simple implementation
    % ====================================
    basis = basiss;
    basis_n = basiss_n;
    basis_sum = normSum;

    
    z_minus_center_sigma = bsxfun(@times, -z_minus_center, 1./(sigma.^2) );

    %> derivative of the basis
    basisD     =  z_minus_center_sigma .* basis; % closed form

    % normalizing basisD
    basisD_sum = sum(basisD,2);
    basisD_n_a = bsxfun(@times, basisD, basis_sum);
    basisD_n_b = bsxfun(@times, basis, basisD_sum);
    basisD_n   = bsxfun(@times, basisD_n_a - basisD_n_b, 1 ./(basis_sum.^2) );


    %> second derivative of the basis
    tmp =  bsxfun(@times,basis, -1./(sigma.^2) );
    basisDD = tmp + z_minus_center_sigma .* basisD;
    basisDD_sum = sum(basisDD,2);

    % normalizing basisDD
    basisDD_n_a  = bsxfun(@times, basisDD, basis_sum.^2);
    basisDD_n_b1 = bsxfun(@times, basisD, basis_sum);
    basisDD_n_b  = bsxfun(@times, basisDD_n_b1, basisD_sum);
    basisDD_n_c1 =  2 * basisD_sum.^2 - basis_sum .* basisDD_sum;
    basisDD_n_c  = bsxfun(@times, basis,  basisDD_n_c1);
    basisDD_n_d  = basisDD_n_a - 2 .* basisDD_n_b + basisDD_n_c;                                              
    basisDD_n    = bsxfun(@times, basisDD_n_d, 1 ./ basis_sum.^3);


    % This is extra code that can be found in 
    % +TrajectoryGenerators/NormalizedGaussiansBasisGenerator
    % I do not know why it is used, and it seems Gdd_n is not
    % correct
    basisDD_n = bsxfun(@times, basisDD_n, zd.^2) + bsxfun(@times, basisD_n, zdd );
    basisD_n =  bsxfun(@times, basisD_n,  zd);


    obj.basis.Gn    = basis_n;
    obj.basis.Gn_d  = basisD_n;
    obj.basis.Gn_dd = basisDD_n;


end





function [] = debugBasis(z, basis_n, basiss_n)

sb = struct('Color', 'b', 'LineStyle', '-', 'LineWidth', 1, 'Marker' , 'o',...
                   'MarkerFaceColor','b', 'MarkerSize',  5);
sr = struct('Color', 'r', 'LineStyle', '-', 'LineWidth', 1, 'Marker' , 'o',...
                   'MarkerFaceColor','r', 'MarkerSize',  2.5);

    figurew 'asdf'
    plot( z,  basis_n,  sb)
    plot( z,  basiss_n, sr)
    
% z: [100x1]
% z: [100x30]

    keyboard
    figurew 'asdf'
    plot( z_minus_center)
    
    figurew 'asdf'
    plot(basis_sum)
    
end


















