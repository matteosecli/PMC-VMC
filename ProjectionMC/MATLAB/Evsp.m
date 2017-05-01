function [E0, E0_error, p_interval ] = Evsp( EL, bx )
%Evsp plots the GS energy vs the PMC simulation length p.
%   Detailed explanation goes here


%% Initialization

    % Define the binning lengths
    %Nbin_interval = [30:1:1000 1050:50:10000 10500:500:20000] ;
    Nbin_interval = [25000:1000:500000] ;
    %Nbin_interval = [25000:1000:100000];
    Length_interval = floor(length(EL)./Nbin_interval);
    p_interval = Length_interval-1;

    % Do a preallocation of the vector for speed
    E0 = zeros(1,length(Nbin_interval));
    E0_error = zeros(1,length(Nbin_interval));
    
    % Rescale bx so to avoid overflow
    bx = bx/mean(bx);
    
    for NumBinIdx = 1:length(Nbin_interval)
        fprintf( 'Binning technique: doing the calculations for %d bins.\n', Nbin_interval(NumBinIdx));
        E0_binned = zeros(1,Nbin_interval(NumBinIdx));
        Gnp = zeros(1,Nbin_interval(NumBinIdx));
        for BinIdx = 1:Nbin_interval(NumBinIdx)
            Gnp(BinIdx) = prod(bx(1+Length_interval(NumBinIdx)*(BinIdx-1):Length_interval(NumBinIdx)*BinIdx-1));
            E0_binned(BinIdx) = Gnp(BinIdx)*EL(Length_interval(NumBinIdx)*BinIdx);
        end
        E0(NumBinIdx) = sum(E0_binned)/sum(Gnp);
        %E0_error(NumBinIdx) = E0(NumBinIdx)*sqrt( (std(E0_binned)/mean(E0_binned))^2 + (std(Gnp)/mean(Gnp))^2 );
        E0_error(NumBinIdx) = E0(NumBinIdx)*sqrt( (binanalysis(E0_binned,'No','No')/mean(E0_binned))^2 + (binanalysis(Gnp,'No','No')/mean(Gnp))^2 );
    end

    
%% Plot
    
    % Create figure
    BinFigure = figure('PaperOrientation','landscape','PaperType','A3');

    % Create axes
    BinAxes = axes('Parent',BinFigure);
    hold(BinAxes,'on');

    % Create plot
    errorbar(p_interval, E0, E0_error,...
        'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
        'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
        'MarkerSize',4,...
        'Marker','square',...
        'LineWidth',2,...
        'LineStyle','none');

    % Helper function for labels
%        varname=@(var) inputname(1);

    % Create xlabel
    xlabel('$p$','FontSize',22,'Interpreter','latex');

    % Create title
    title('$E_{0}$ vs $p$','FontSize',24,'Interpreter','latex');

    % Create ylabel
    ylabel('$E_{0}$','FontSize',22,'Interpreter','latex');

    box(BinAxes,'on');
    % Set the remaining axes properties
    set(BinAxes,'FontSize',18);%,'YMinorTick','on','YScale','log');
    

end
