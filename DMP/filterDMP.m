%% Function that filters q,qd,qdd with filtfilt to construct DMP

function [qFil,qdFil,qddFil] = filterDMP(t,q,qd,cutoff)

    dof = size(q,1);
    assert(dof == 7,'Degree of freedom is not 7!');
    
    if ~exist('filtfilt')
        
        % Fourier transform/inverse transform filter approach
        %{
        L = 2048; %500 * ceil(length(t)/500);
        NFFT = 2^nextpow2(L);
        Fs = L/(t(end)-t(1));
        f = Fs/2*linspace(0,1,NFFT/2+1);
        tLinStrike = linspace(t(1),t(end),L);
        qLin = interp1(t,q',tLinStrike);
        qdLin = interp1(t,qd',tLinStrike);
        qSpecs = fft(qLin,NFFT)/L;
        qdSpecs = fft(qdLin,NFFT)/L;

        % plot single-sided amplitude spectrum.
        %{
        figure;
        for j = 1:dof
            subplot(7,2,2*j-1);
            plot(f,2*abs(qSpecs(1:NFFT/2+1,j)));
            legend(['spectra\_',joints{j}]);
            subplot(7,2,2*j);
            plot(f,2*abs(qdSpecs(1:NFFT/2+1,j)));
            legend(['spectra\_',vel{j}]);
        end
        %}
        qSpecsZeroed = [qSpecs(f < cutoff,:); zeros(sum(f >= cutoff),dof)];
        qdSpecsZeroed = [qdSpecs(f < cutoff,:); zeros(sum(f >= cutoff),dof)];
        qFilLin = real(ifft(qSpecsZeroed,NFFT)*2*L);
        qdFilLin = real(ifft(qdSpecsZeroed,NFFT)*2*L);
        % interpolate back
        qFil = interp1(tLinStrike,qFilLin,t);
        qdFil = interp1(tLinStrike,qdFilLin,t);
        qFil = qFil'; 
        qdFil = qdFil';
        %}
        w_nyquist = floor(length(t)/2);
        w_cut = cutoff/w_nyquist;
        % filter, reverse data, filter again, and reverse data again
        qFil = q;
        qdFil = qd;
        for i = 1:dof
            q2 = filterButter2nd(q(i,:),w_cut);
            q2 = q2(end:-1:1);
            q2 = filterButter2nd(q2,w_cut);
            qFil(i,:) = q2(end:-1:1);
            qd2 = filterButter2nd(qd(i,:),w_cut);
            qd2 = qd2(end:-1:1);
            qd2 = filterButter2nd(qd2,w_cut);
            qdFil(i,:) = qd2(end:-1:1);
        end
    else
        w_nyquist = floor(length(t)/2);
        w_cut = cutoff/w_nyquist;
        [B,A] = butter(2,w_cut);
        qFil = filtfilt(B,A,q'); 
        qdFil = filtfilt(B,A,qd'); 
        qFil = qFil'; 
        qdFil = qdFil';
    end    

    qdd = diff(qdFil')./(repmat(diff(t),1,dof));
    qdd = [qdd; qdd(end,:)]';
    
    if ~exist('filtfilt')
        %{
        qddSpecs = fft(qdd',NFFT)/L;
        qddSpecsZeroed = [qddSpecs(f < cutoff,:); zeros(sum(f >= cutoff),dof)];
        qddFilLin = real(ifft(qddSpecsZeroed,NFFT)*2*L);
        % interpolate back
        qddFil = interp1(tLinStrike,qddFilLin,t);
        qddFil = qddFil';
        %}
        qddFil = qdd;
        for i = 1:dof
            qdd2 = filterButter2nd(qdd(i,:),w_cut);
            qdd2 = qdd2(end:-1:1);
            qdd2 = filterButter2nd(qdd2,w_cut);
            qddFil(i,:) = qdd2(end:-1:1);
        end
    else
        qddFil = filtfilt(B,A,qdd');
        qddFil = qddFil';
    end

end