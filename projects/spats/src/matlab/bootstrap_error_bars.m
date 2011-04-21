% This program bootstraps the data observed for a given RNA 
% molecule, by drawing samples from the empirical length 
% distibution observed for that molecule. The program then 
% computes the sample-based variance to generate error bars. 
close all
clear all

RandStream.setDefaultStream(RandStream('mt19937ar', 'seed', sum(100*clock)));

% ## ALGORITHM INPUT PARAMS ##
% The optimization flag options are: 
% 1 = set negatives to zero and normalize, 
% 2 = apply numerical optimization.
optim_flag = 2;  
% This is the number of iterations to be used by the numerical optimization
% procedure.
num_iter = 20;  % this is the max number of iterations allowed.
accuracy = 0.0000001;  % this is the desired optimization accuracy.
num_samples = 500; % this is the number of bootstrap samples
n = 172; % number of nucleotides to be read from the data file

basedir = '';

empirical_counts_data = zeros(n, 4); 
empirical_signal_data = zeros(n, 4); 
signal_data =  zeros(n, 4); 
ML_Theta = zeros(n-1, num_samples);
ML_Gamma =  zeros(n-1, num_samples);
ML_c = zeros(1, num_samples);
intermediate_c = zeros(1, num_samples);
p_0_hats = zeros(1, num_samples);
status_flags = zeros(1, num_samples);
sum_negatives = zeros(1, num_samples);
num_negatives = zeros(1, num_samples);

% ## RETRIEVE OBSERVED COUNTS ## 
fid = fopen(strcat(basedir, 'target_WT.adducts'), 'r');     % n=198

% Skip the contents of the first line since it pertains to complete 
% fragments (we retrieve to advance the file pointer). 
first_line = fgetl(fid);
for i=1:n
    % there are 4 columns per line, where the last two contain the counts
    empirical_counts_data(i,:) = fscanf(fid, '%d %s\t %d %d', [1 4]);  
end
status = fclose(fid);

% Determine sample sizes for both channels.
experiment_sample_size = sum(empirical_counts_data(:, 3));  
control_sample_size = sum(empirical_counts_data(:, 4));  

% Organize data by direction of signal decay (i.e., from position 1 to 
% position n) and compute the empirical frequencies.
empirical_signal_data = flipud(empirical_counts_data);

% Construct the empirical length distribution
experiment_length_freqs = empirical_signal_data(1:n,3)/sum(empirical_signal_data(:,3));
control_length_freqs = empirical_signal_data(1:n,4)/sum(empirical_signal_data(:,4));

% ## LOOP OVER NUMEROUS SAMPLES ##
for counter=1:num_samples
     
    signal_data(:,1:2) = flipud(empirical_counts_data(:,1:2));
    
    % ## DRAW RANDOM COUNTS FROM EMPIRICAL DISTRIBUTIONS ##
    experiment_random_indices = randp(experiment_length_freqs, experiment_sample_size, 1);
    control_random_indices = randp(control_length_freqs, control_sample_size, 1);
    % NOTE: this data is already given in the direction of signal decay 
    % (i.e., it corresponds to positions 1 to n).
    experiment_random_counts = accumarray(experiment_random_indices, 1, [n 1]);
    signal_data(:, 3) = experiment_random_counts;
    control_random_counts = accumarray(control_random_indices, 1, [n 1]);
    signal_data(:, 4) = control_random_counts;
    
    % ## PRE-PROCESS DATA ##
    % Here, we omit from analysis all nucleotides that had zero counts in both
    % channels (X_k = Y_k = 0). Their reactivities are automatically set to 0. 
    % This way, we make sure each position k in the resulting sequence has at
    % least one non-zero count. This is required for the feasibility and 
    % stability of the numerical optimization procedure. Later, we re-assemble 
    % the entire sequence. 
    analyzed_ind = sort(union(find(signal_data(:, 3) > 0), find(signal_data(:, 4) > 0)));
    analyzed_data = zeros(length(analyzed_ind), 4);
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
    opt_Theta = zeros(N-1, 1);
    opt_Gamma = zeros(N-1, 1);
                                
    % ## COMPUTE MLE ##
    [opt_Theta, opt_Gamma, ML_c(counter), intermediate_c(counter), sum_negatives(counter), num_negatives(counter), p_0_hats(counter), status_flags(counter)] = MLE_one_molecule(analyzed_data, optim_flag, num_iter, accuracy);
                 
    % ## RE-ASSEMBLE ORIGINAL SEQUENCE ##
    ML_Theta(analyzed_ind(1:N-1),counter) = opt_Theta;
    ML_Gamma(analyzed_ind(1:N-1),counter) = opt_Gamma;
      
end

% ## COMPUTE SAMPLE-BASED VARIANCE ## 
% Screen out all samples where data quality was too bad
valid_final_Theta = ML_Theta(:, find(status_flags > 0));
size(valid_final_Theta)

valid_final_c = ML_c(find(status_flags > 0));

bootstrap_Theta_std = std(valid_final_Theta, 0, 2);
bootstrap_final_c_std = std(valid_final_c);
bootstrap_intermediate_c_std = std(intermediate_c);

save workspace


