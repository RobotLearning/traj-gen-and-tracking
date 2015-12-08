function [Sigma, cholSigma] = regularizeCovariance(Sigma, priorCov, numEffectiveSamples, priorCovWeight)


    count = 1;
    while (count < 100)
        Sigma_temp = (Sigma * numEffectiveSamples + priorCov * priorCovWeight) / (numEffectiveSamples + priorCovWeight);
        
        try 
            cholSigma = chol(Sigma_temp);
            Sigma = Sigma_temp;
            return
        catch E
            priorCovWeight = priorCovWeight * 2;
        end
        count = count + 1;
    end
    disp(Sigma); disp(eig(Sigma))
    error('Could not find decomposition for covariance matrix... HELP!!!\n');


end