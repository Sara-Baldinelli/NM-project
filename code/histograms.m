% histograms for RSSA and RTC algorithms comparison
% x axis: voltage bins
% y axis: proportion of open calcium + potassium channels in a specific
% voltage bin over the total number of open channels over all bins
Ctot_values = [2, 5, 40];
num_bins = 100;  % define number of bins of the histogram

voltage_min = -90;  % minimum voltage (just for binning)
voltage_max = 90;   % maximum voltage 
bin_edges = linspace(voltage_min, voltage_max, num_bins + 1);

tmax = 200000;            % duration of each simulation
num_simulations = 1;  % number of simulations for each algorithm

for j = 1:length(Ctot_values)
    Ctot = Ctot_values(j);

    bin_totals_RSSA = zeros(num_bins, 1);   % to sum open channels for each bin
    total_open_channels_RSSA = 0;           % to sum open channels over all the bins

    bin_totals_2 = zeros(num_bins, 1);      % same for paper (RTC) simulation
    total_open_channels_2 = 0;
    
    % RSSA simulation
    for sim = 1:num_simulations
        [T, Dynamics] = MOESM4_RSSA_2(tmax, Ctot, Ctot, 100, 0.1, 0.001);
    
        voltage = Dynamics(:, 5);                         % voltage values
        open_channels = Dynamics(:, 3) + Dynamics(:, 4);  % total open channels
                                                                
        [bin_counts, ~, bin_indices] = histcounts(voltage, bin_edges);  % bin voltage values
        % bin_counts is a vector that contains the count of elements from voltage that fall into each of the bins
        % bin_edges is a vector that contains the edges of each bin (the total number of edges is 101)
        % bin_indices is a vector that tells you which bin each voltage value belongs to
        
        for i = 1:num_bins
            bin_totals_RSSA(i) = bin_totals_RSSA(i) + sum(open_channels(bin_indices == i));  % counts open channels for each bin
        end
        total_open_channels_RSSA = total_open_channels_RSSA + sum(open_channels);            % total number of open channels
    end


    % paper (RTC) simulation
    for sim = 1:num_simulations
        [V, M, N, t, Mtot, Ntot] = mlexactboth(tmax, Ctot, Ctot); 
    
        voltage_2 = V;                  
        open_channels_2 = M + N;        
       
        [bin_counts_2, ~, bin_indices_2] = histcounts(voltage_2, bin_edges);
        
        for i = 1:num_bins
            bin_totals_2(i) = bin_totals_2(i) + sum(open_channels_2(bin_indices_2 == i)); 
        end
        total_open_channels_2 = total_open_channels_2 + sum(open_channels_2);
    end

    % proportions of open channels for each bin
    bin_proportions_RSSA = bin_totals_RSSA / total_open_channels_RSSA;
    bin_proportions_2 = bin_totals_2 / total_open_channels_2;

    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;  % (left edges + right edges)/2 = midpoints

    figure;
    subplot(2, 1, 1); 
    bar(bin_centers, bin_proportions_RSSA, 'hist');
    xlabel('Voltage (mV)');
    ylabel('Proportion of Total Open Channels');
    title(['Total Open Channels vs Voltage (RSSA, Mtot = Ntot = ', num2str(Ctot), ')']);
    xlim([-100 100]);
    
    subplot(2, 1, 2); 
    bar(bin_centers, bin_proportions_2, 'hist');
    xlabel('Voltage (mV)');
    ylabel('Proportion of Total Open Channels');
    title(['Total Open Channels vs Voltage (RTC, Mtot = Ntot = ', num2str(Ctot), ')']);
    xlim([-100 100]);
end