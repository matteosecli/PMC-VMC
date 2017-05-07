function [x0, x0_error] = meanposition( x, bx, simlen )
%Evsp plots the GS energy vs the PMC simulation length p.
%   Detailed explanation goes here


%% Initialization

    % Define p and number of bins
    p = simlen-1;
    Nbins = floor(length(x)/simlen);
    
    % Rescale bx so to avoid overflow
    bx = bx/mean(bx);
    
    % Treat x as doubles
    x = double(x);

    x0_binned = zeros(1,Nbins);
    Gnp = zeros(1,Nbins);
    for BinIdx = 1:Nbins
        Gnp(BinIdx) = prod(bx(1+simlen*(BinIdx-1):simlen*BinIdx-1));
        x0_binned(BinIdx) = Gnp(BinIdx)*x(floor(simlen*(BinIdx-0.5)));
    end
    x0 = sum(x0_binned)/sum(Gnp);
    x0_error = x0*sqrt( (binanalysis(x0_binned,'No','No')/mean(x0_binned))^2 + (binanalysis(Gnp,'No','No')/mean(Gnp))^2 );

end
