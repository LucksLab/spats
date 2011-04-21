close all
clear all

n = 172; % number of nucleotides in the data file
counts_data = zeros(n, 4); 

basedir = '';

% ## READ COUNTS ##
% File format:   Index     Base    (+)-Channel_Count     (-)-Channel_Count
% Nucleotides are ordered from 5' to 3', top to bottom. 

% Since the signal decays from 3' to 5', this means that we read the 
% signal from end to start, that is, the data start at position k=n, 
% and end at position k=1.
fid = fopen(strcat(basedir, 'target_WT.adducts'), 'r');     % n=198

% # SEQUENCE COMPOSITIONS #
% RNase P: molecule starts at ID=3 (following an added GG at 5') and ends at ID=156.
% pT181 112: molecule starts at ID=15 and ends at ID=112. 
%            1-14 - 5' Structure Cassette (14 nt's),
%            15-126 - RNA molecule (112 nt's),
%            127-173 - 3' Structure Cassete (46 nt's).
% In signal representation (the opposite direction): the molecule ranges 
% from nucleotides 47 to 158.

% get the first line to advance the file pointer to the data section
first_line = fgetl(fid);

for i=1:n
    % there are 4 columns per line, where the last two contain the counts
    counts_data(i,:) = fscanf(fid, '%d %s\t %d %d', [1 4]);  
end
status = fclose(fid);

% Organize data by direction of signal decay (i.e., from position 1 to 
% position n).
signal_data = flipud(counts_data);

% ## ILLUSTRATE RAW DATA (as a sanity check) ##
figure;
signal_freq = signal_data(:,3)/sum(signal_data(1:end, 3));
noise_freq = signal_data(:,4)/sum(signal_data(1:end,4));
diff_signal = (signal_data(:,3) - signal_data(:,4));
diff_signal(find(diff_signal < 0)) = 0;
diff_freq = diff_signal/sum(diff_signal);
bar([signal_freq noise_freq diff_freq], 'group');
set(gca,'xlim',[1 n]);
title('Raw fragment frequencies from both channels + their diff. distribution');
legend('(+) Channel', '(-) Channel', 'Subtracted counts');
colormap Jet

% ## PRE-PROCESS DATA ##
% Here, we omit from analysis all nucleotides that had zero counts in both
% channels (X_k = Y_k = 0). Their reactivities are automatically set to 0. 
% This way, we make sure each position k in the resulting sequence has at
% least one non-zero count. This is required for the feasibility and 
% stability of the numerical optimization procedure. Later, we re-assemble 
% the entire sequence. 
analyzed_ind = sort(union(find(signal_data(:, 3) > 0), find(signal_data(:, 4) > 0)));
analyzed_data = signal_data(analyzed_ind, :);

% Locate all positions for which X_k=0 or Y_k=0, and set the counts to 1.
% This correction is applied for numerical stability, but can be omitted
% later!
zero_Ys = find(analyzed_data(:,4) == 0);
analyzed_data(zero_Ys, 4) = 1;
zero_Xs = find(analyzed_data(:,3) == 0);
analyzed_data(zero_Xs, 3) = 1;

% This is the length of the analyzed sequence, including the complete
% fragments counts, which are at the end.
N = length(analyzed_ind); 

% ## COMPUTE MLE ##
% ## OPTION 1: APPLY ALGORITHM 1 ##
% We first transform the counts into frequencies.
% Note that position "N" (last line) pertains to counts of complete
% fragments, and should be treated differently than all other counts!
% These are the w_k's.
Plus_Channel_Freq_Data = analyzed_data(1:N,3)/sum(analyzed_data(1:N,3));
X_0 = analyzed_data(N,3);
% These are the p_k's.
Minus_Channel_Freq_Data = analyzed_data(1:(N-1),4)/sum(analyzed_data(1:N,4));
p_0_hat = analyzed_data(N,4)/sum(analyzed_data(1:N,4));

% # ILLUSTRATE SIGNAL + NOISE DECAY
figure;
subplot(2,1,1);
embedded_Plus_Channel_Freq = zeros(n,1);
embedded_Plus_Channel_Freq(analyzed_ind(1:N)) = Plus_Channel_Freq_Data;
bar(embedded_Plus_Channel_Freq);
set(gca, 'xlim', [1 n]);
title('(+) Channel Input (X length distribution)');
colormap Winter
subplot(2,1,2);
embedded_Minus_Channel_Freq = zeros(n,1);
embedded_Minus_Channel_Freq(analyzed_ind(1:N-1)) = Minus_Channel_Freq_Data;
embedded_Minus_Channel_Freq(n) = p_0_hat;
bar(embedded_Minus_Channel_Freq);
set(gca, 'xlim', [1 n]);
title('(-) Channel Input (Y length distribution)');
colormap Winter

% Compute the gammas.
current_dropoff_rates = zeros(N-1,1);
for j=1:N-1
    current_dropoff_rates(j) = analyzed_data(j,4)/sum(analyzed_data(j:N,4));
end

% consistency check
if ((X_0 == 0) || (p_0_hat == 0))
    disp('Zero complete fragments in at least one of the channels!');
    disp('To continue analysis, type the command return and then press Enter.');
    %waitforbuttonpress
    keyboard
end

% We now recover the Poisson rate estimate (c_hat).
c_estimate = log(p_0_hat) - log(X_0/sum(analyzed_data(:,3)));

% consistency check
if (c_estimate < 0)
    disp('Estimated Poisson rate is negative!');
    disp('To continue analysis, type the command return and then press Enter.');
    %waitforbuttonpress
    keyboard
end

% We now reconstruct the "initial distribution" Theta.
signal_portion = zeros(N-1,1);
for m=1:(N-1)
    signal_portion(m) = log(1 + (Plus_Channel_Freq_Data(m)/(sum(Plus_Channel_Freq_Data((m+1):end)))));
end
noise_portion = log(1 - current_dropoff_rates);
recovered_Theta = (1/c_estimate)*(signal_portion + noise_portion);
% NORMALIZE after setting all negatives to zero.
recovered_Theta(find(recovered_Theta < 0)) = 0;
recovered_Theta = recovered_Theta/sum(recovered_Theta);

% This is the noise decay.
noise_decay_factors = -Minus_Channel_Freq_Data./noise_portion;

figure;
subplot(2, 1, 1);
bar(Minus_Channel_Freq_Data, 'hist');
colormap Winter
set(gca,'xlim',[1 N-1]);
title('Empirical Noise Length Distribution');

subplot(2, 1, 2);
plot((1:1:N-1), noise_decay_factors, 'm*-');
xlim([1 N-1]);
title('Noise Decay Curve');

% This is the signal decay.
signal_decay_factors = Plus_Channel_Freq_Data(1:N-1)./signal_portion;

figure;
subplot(2, 1, 1);
bar(Plus_Channel_Freq_Data(1:N-1), 'hist');
colormap Winter
set(gca,'xlim',[1 N-1]);
title('Empirical Signal Length Distribution');

subplot(2, 1, 2);
plot((1:1:N-1), signal_decay_factors, 'm*-');
xlim([1 N-1]);
title('Signal Decay Curve');

% Illustrate both corrections together
figure;
plot((1:1:N-1), signal_decay_factors, 'm*-', (1:1:N-1), noise_decay_factors, 'c*-');
xlim([1 N-1]);

% Illustrate the ideal case
ideal_theta_vec = (1/(N-1))*ones(1,N-1);
k_fragment_prob = zeros(1,N-1);  

% Compute the empirical length distribution.
c = 1;
elongation_const_rate = 0.995;
for m=1:(N-1)
    k_fragment_prob(m) = exp(-c)*(exp(c*sum(ideal_theta_vec(1,m:N-1))*((elongation_const_rate)^(m-1))) - exp(c*sum(ideal_theta_vec(1,m+1:N-1))*((elongation_const_rate)^(m))));
end
prob_no_modification = 1-sum(k_fragment_prob);
k_fragment_prob = k_fragment_prob/sum(k_fragment_prob);
ideal_decay_factors = k_fragment_prob./ideal_theta_vec;

figure;
plot((1:1:N-1), ideal_decay_factors/max(ideal_decay_factors), 'm*-');
xlim([1 N-1]);
title('Model-Based Noisy Signal Decay Curve for Ideal Uniform Distribution');

figure;
subplot(2, 1, 1);
bar(ideal_theta_vec, 'hist');
colormap Winter
set(gca,'xlim',[1 N-1]);
title('Uniform Distribution');

subplot(2, 1, 2);
plot((1:1:N-1), ideal_decay_factors/max(ideal_decay_factors), 'm*-');
xlim([1 N-1]);
title('Signal Decay Curve for a Uniform Distribution');

