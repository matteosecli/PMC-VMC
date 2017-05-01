function Error = binanalysis( data, doplot, printmsg )
%BINANALYSIS gives the error on 'data' via the binning technique.
%   Detailed explanation goes here


%% Initialization

    % Decide whether to print msgs or not
    printmsgBOOL = true(1);
    printmsgBOOL = 0;
    if ( strcmp(printmsg, 'Yes') )
        printmsgBOOL = 1;
    end

    % Define the binning lengths
    %Nbin_interval = [30:1:1000 1050:50:10000 10500:500:20000] ;
    Nbin_interval = [10:1:100];
    Length_interval = floor(length(data)./Nbin_interval);

    % Do a preallocation of the vector for speed
    x_errors = zeros(1,length(Nbin_interval));
    
    for NumBinIdx = 1:length(Nbin_interval)
        if printmsgBOOL
            fprintf( 'Binning technique: doing the calculations for %d bins.\n', Nbin_interval(NumBinIdx));
        end
        x_binned_mean = zeros(1,Nbin_interval(NumBinIdx));
        for BinIdx = 1:Nbin_interval(NumBinIdx)
            x_binned_mean(BinIdx) = sum(data(1+Length_interval(NumBinIdx)*(BinIdx-1):Length_interval(NumBinIdx)*BinIdx))/Length_interval(NumBinIdx);
        end
        bin_std = std(x_binned_mean);
        x_errors(NumBinIdx) = bin_std/sqrt(Nbin_interval(NumBinIdx));
    end

    % Calculate the error as the max of the errors
    Error = max(x_errors);

    
%% Plot if asked so
    
    if ( strcmp(doplot, 'Yes') )
        % Create figure
        BinFigure = figure('PaperOrientation','landscape','PaperType','A3');

        % Create axes
        BinAxes = axes('Parent',BinFigure);
        hold(BinAxes,'on');

        % Create plot
        plot(Length_interval, x_errors,'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
            'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
            'MarkerSize',4,...
            'Marker','square',...
            'LineWidth',2,...
            'LineStyle','none');

        % Helper function for labels
%        varname=@(var) inputname(1);
        
        % Create xlabel
        xlabel('Bin length','FontSize',24,'Interpreter','latex');

        % Create title
        title(strcat('$\sigma_{',inputname(1),'}$ vs bin length'),'FontSize',24,'Interpreter','latex');

        % Create ylabel
        ylabel(strcat('$\sigma_{',inputname(1),'}$'),'FontSize',24,'Interpreter','latex');

        box(BinAxes,'on');
        % Set the remaining axes properties
        set(BinAxes,'FontSize',18);%,'YMinorTick','on','YScale','log');
    end
    

end