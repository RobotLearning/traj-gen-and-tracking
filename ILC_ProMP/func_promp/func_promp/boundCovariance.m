function [Sigma] = boundCovariance(Sigma, minCov, minCorr)

    if (length(minCov) > 0)
        for d = 1 : size(Sigma,1)

            if( isnan(Sigma(d,d) ) )
                Sigma(d,d) = minCov(d);
            end

            scaling = max(minCov(d) / Sigma(d,d), 1);

            if (~isinf(scaling))
                Sigma(:,d) = Sigma(:,d) .* sqrt(scaling);
                Sigma(d,:) = Sigma(d,:) .* sqrt(scaling);

                if (minCorr < 1)
                    for r = 1 : size(Sigma,1)
                        if(d ~= r)
                            %limit corrCoeff
                            corrCoeff = Sigma(d,r) / sqrt( Sigma(r,r) * Sigma(d,d) );
                            %assert(~isnan(corrCoeff));
                            if(~isnan(corrCoeff) && abs(corrCoeff) > minCorr)
                                Sigma(d,r) = sign(corrCoeff) * minCorr * sqrt( Sigma(r,r) * Sigma(d,d) );
                                Sigma(r,d) = Sigma(d,r);
                            end
                        end
                    end
                end
            else
                Sigma(d,d) = minCov(d);
            end

            if(Sigma(d,d) < minCov(d) )
                Sigma(d,d) = Sigma(d,d) + minCov(d);
            end

        end
    end

    Sigma(:,:) = (Sigma(:,:) + Sigma(:,:)' ) ./2; %Enforce symmetry

    [V,D] = eig(Sigma);
    Sigma = V * max(D,0) * V';
    
    %Fiddle with eigenvalues

%    eigenVal = eig(Sigma(:,:));
%    while (abs(min(eigenVal)/max(eigenVal)) <= 1e-12)
%        Sigma(:,:) = Sigma(:,:) + eye(size(Sigma(:,:))) * max(eigenVal) * 10^-10;
%        eigenVal = eig(Sigma(:,:));
%    end %while
% end %options

%    end %while
end
