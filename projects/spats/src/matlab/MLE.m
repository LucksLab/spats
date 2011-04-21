% This program reads a file with fragment counts and constructs 
% the estimated SHAPE relative reactivities from the counts. 
close all
clear all

global X_counts
global X_0
global Y_counts
global current_num_weights_vec
global current_den_weights_vec
global current_elongation_rates
global current_contributing_ind
global current_c 
global current_Theta
global current_c_weight 

max_iter = 10;

% ## READ COUNTS ##
% File format:   Index     Base    (+)-Channel_Count     (-)-Channel_Count
% Nucleotides are ordered from 5' to 3', top to bottom. 

% Since the signal decays from 3' to 5', this means that we read the 
% signal from end to start, that is, the data start at position k=n, 
% and end at position k=1.
% n=172 / n=192 (for pT181_Sense_122/132.adducts)
% n=198 for RNase P specificity domain
n = 172; % number of nucleotides in the data file
counts_data = zeros(n, 4); 

basedir = '';

fid = fopen(strcat(basedir, 'pT181_Sense_112.adducts'), 'r');  

% # SEQUENCE COMPOSITIONS #
% RNase P: molecule starts at ID=3 (following an added GG at 5') and ends at ID=156.
% pT181 112: molecule starts at ID=15 and ends at ID=126. 
%            1-14 - 5' Structure Cassette (14 nt's),
%            15-126 - RNA molecule (112 nt's),
%            127-172 - 3' Structure Cassete (46 nt's).
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

% ## ILLUSTRATE RAW DATA ##
figure;
signal_freq = signal_data(:,3)/sum(signal_data(1:end, 3));
noise_freq = signal_data(:,4)/sum(signal_data(1:end,4));
diff_signal = (signal_data(:,3) - signal_data(:,4));
diff_signal(find(diff_signal < 0)) = 0;
diff_freq = diff_signal/sum(diff_signal);
bar([signal_freq noise_freq diff_freq], 'group');
set(gca, 'xlim', [1 n]);
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
zero_Ys = find(analyzed_data(:,4) == 0);
analyzed_data(zero_Ys, 4) = 1;
zero_Xs = find(analyzed_data(:,3) == 0);
analyzed_data(zero_Xs, 3) = 1;

% This is the length of the analyzed sequence, including the complete
% fragments counts, which are at the end.
N = length(analyzed_ind); 

% Storage for the algorithm's progress reports
Theta_estimates = zeros(N-1, max_iter+1);
Gamma_estimates =  zeros(N-1, max_iter+1);
c_estimates = zeros(1, max_iter+1);
log_lik_vals = zeros(3, max_iter+1);

% ## COMPUTE MLE ##
% ## APPLY ALGORITHM 1 ##
% We first transform the counts into frequencies.
% Note that position "N" (last line) pertains to counts of complete
% fragments, and should be treated differently than all other counts!
Plus_Channel_Freq_Data = analyzed_data(1:N,3)/sum(analyzed_data(1:N,3));
X_0 = analyzed_data(N,3);
Minus_Channel_Freq_Data = analyzed_data(1:(N-1),4)/sum(analyzed_data(:,4));
p_0_hat = analyzed_data(N,4)/sum(analyzed_data(:,4));

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

% ## RESOLVE INCONSISTENCIES ## 
% ## APPLY ALGORITHM 2 ##
% First, zero all negative entries, and normalize the rest.
prelim_Theta = recovered_Theta;  % save for illustration purposes
negative_entries = find(recovered_Theta < 0);
delta = 1 - sum(recovered_Theta(negative_entries));
recovered_Theta(negative_entries) = 0;
recovered_Theta = recovered_Theta/sum(recovered_Theta);

% Numerical optimization. 
% Split into 3 local convex optimization problems, and iteratively solve 
% them until convergence is observed or a maximum number of iterations 
% is reached.

% Step 0: Initialize.
% Save initial values.
Theta_estimates(:,1) = recovered_Theta;
Gamma_estimates(:,1) = current_dropoff_rates;
c_estimates(1) = c_estimate;
WF_Theta = zeros(N-1,1);

% These quantities remain fixed throughout all iterations.
X_counts = analyzed_data(1:N-1, 3)';   % The X array does NOT include X_0.
Y_counts = analyzed_data(1:N, 4)';   % The Y array includes Y_0 as its last entry.
% This array is used by the water filling procedure.
partial_X_count_sums = zeros(N-1,1);
% This array is used by the Gamma optimization procedure.
S_vec = zeros(N,1);  
total_combined_sum = sum(Y_counts) + sum(X_counts) + analyzed_data(N, 3);
S_vec(1) = total_combined_sum;  % This is S_0
for k=2:(N-1)
    partial_X_count_sums(k) = sum(X_counts(1:k-1));
    S_vec(k) = total_combined_sum - sum(Y_counts(1:k-1)) - partial_X_count_sums(k);
end
S_vec(N) = total_combined_sum - sum(Y_counts(1:N-1)) - sum(X_counts);  % This is S_n = X_0 + Y_0

% These matrices are used by the water filling routine.
WF_num_matrix = ones(N-1,N-1)*diag(partial_X_count_sums, 0);
WF_den_matrix = [WF_num_matrix(:,2:end) sum(X_counts(1:N-1))*ones(N-1,1)];

log_lik_vals(:,1) = ones(3,1)*log_likelihood_fun(recovered_Theta, current_dropoff_rates, c_estimate);

% Set the values for iteration 1.
current_c = c_estimate;
current_Theta = recovered_Theta;
current_T_vec = exp(current_c*current_Theta) - 1;
current_elongation_rates = 1 - current_dropoff_rates;
iter = 1;
relative_log_change = 1;

% LOOP
tic;
while ((iter <= max_iter) && (relative_log_change > 0.0000001))
    
    disp(sprintf('Starting iteration no. %d', iter));
    
    % Step 1: Solve P1 = maximize H(c) over R+. 
    % Here, Gamma and Theta are kept fixed.
    current_c_weight = sum(current_Theta.*partial_X_count_sums) - sum(X_counts) - X_0;

    disp('Solving P1...');
    
    % minimize -H(c) without the non-negativity constraint
    disp('Debugging message: starting fminunc');
    options = optimset('LargeScale', 'on', 'GradObj','on');
    argmin = fminunc(@H_fun, current_c, options);
    
    % We showed that H(c) is strictly concave in c over H(c)'s domain, 
    % which is D = {c | c > "c_thres"}, where "c_thres" can be negative!
    % For this reason, if the unconstrained-maximum is negative, then 
    % we know that the function must decrease when further increasing c, 
    % and hence c=0 is the optimal c value under the constraint c>=0.
    % NOTE: for a valid experiment, we do expect the unconstrained 
    % optimium to be positive!
    if (argmin < 0)
        argmin = 0;
        disp('P1: optimal c is negative --> c was set to 0.');
    end
    current_c = argmin;
    log_lik_vals(1, iter+1) = log_likelihood_fun(current_Theta, current_dropoff_rates, current_c);
    
    % Step 2: Solve P2 = maximize Lg(Gamma) over all Gammas in the unit hypercube. 
    % Here, c and Theta are kept fixed.
    
    disp('Solving P2...');
    
    b_vec = X_counts' + Y_counts(1:N-1)' - current_T_vec.*(S_vec(2:N) + Y_counts(1:N-1)');
    discriminant_vec = b_vec.*b_vec + 4.*current_T_vec.*Y_counts(1:N-1)'.*S_vec(1:N-1);
    current_dropoff_rates = (b_vec + sqrt(discriminant_vec))./(2*S_vec(1:N-1));
    current_dropoff_rates(find(current_dropoff_rates) < 0) = 0;
    current_elongation_rates = 1 - current_dropoff_rates;
    log_lik_vals(2, iter+1) = log_likelihood_fun(current_Theta, current_dropoff_rates, current_c);

    % Step 3: Solve P3 = maximize F(Theta) over all non-negative Theta's whose elements 
    % sum to 1. Here, c and Gamma are kept fixed. 
    disp('Solving P3...');
    
    % This matrix contains the elongation rates (1-gamma_1),...,(1-gamma_n), 
    % where each row reads like (1-gamma_1),...,(1-gamma_n). 
    current_elongation_matrix = ones(N-1,N-1)*diag(current_elongation_rates,0);
               
    % These are the break points of the water-fill function, where we
    % evaluate it.
    current_alphas = current_c*(partial_X_count_sums + (transpose(X_counts)./current_dropoff_rates));
    % These are the weights to be subtracted from the break-point 
    % values in the numerator and in the denominator. 
    current_WF_num_matrix = current_c*WF_num_matrix;
    current_WF_den_matrix = current_c*WF_den_matrix;
        
    [WaterFill_vals ordered_ind] = Compute_WF(current_alphas, current_elongation_matrix, current_WF_num_matrix, current_WF_den_matrix);
    % Find neighborhood of intersection.    
    [cutoff_ind] = find(WaterFill_vals(:,1) > exp(current_c), 1, 'first');        
        
    ordered_nu_values = WaterFill_vals(:,2);
    % All the indices from the cutoff point (and below) are set to zero.
    zeroed_ind = ordered_ind(cutoff_ind+1:end);
    current_contributing_ind = ordered_ind(1:cutoff_ind);
    WF_Theta(zeroed_ind) = 0;
        
    % Next, we find the \nu* for which the total product equals e^c. 
    % These are the constants to be subtracted from \nu when evaulating
    % the total product (i.e., the water-fill function).
    current_num_weights_vec = current_c*partial_X_count_sums;
    current_den_weights_vec = current_c*[partial_X_count_sums(2:end); sum(X_counts(1:end))];
  
    % The search should not go below c*sum(X_k's) since in this 
    % regime some of the terms might yield negative numbers 
    % (since \nu - c*sum(X_k's) at the denominator can be negative, 
    % while the expression in the numerator is positive). Some of 
    % the terms also yield values less than 1, and so the function 
    % stops increasing. 
    % We also limit the search to be between the pertaining two 
    % adjacent alphas since our implementation of "WF_fun" relies 
    % on using a *fixed* and known set of contributing and of 
    % zeroed indices.
    search_lower_end = max(ordered_nu_values(cutoff_ind), current_c*sum(X_counts) + 0.01);
    search_interval = [search_lower_end, ordered_nu_values(cutoff_ind-1)];
    [nu_star] = fzero(@WF_fun, search_interval);
           
    % Finally, we determine the non-zero thetas.
    opt_quotients = (nu_star - current_num_weights_vec(current_contributing_ind))./(nu_star - current_den_weights_vec(current_contributing_ind));
    WF_Theta(current_contributing_ind) = (log(opt_quotients) + log(current_elongation_rates(current_contributing_ind)))/current_c;
        
    current_Theta = WF_Theta;
    % Update the T_k's for usage in the next iteration
    current_T_vec = exp(current_c*current_Theta) - 1;
    log_lik_vals(3, iter+1) = log_likelihood_fun(current_Theta, current_dropoff_rates, current_c);
        
    % Keep progress records
    c_estimates(iter+1) = current_c;
    Gamma_estimates(:,iter+1) = current_dropoff_rates;
    Theta_estimates(:,iter+1) = current_Theta;
    relative_log_change = -(log_lik_vals(3, iter+1) - log_lik_vals(3, iter))/log_lik_vals(3, iter);
    iter = iter + 1;
end
toc

% ## RE-ASSEMBLE ORIGINAL SEQUENCE ##
all_Thetas = zeros(n-1, iter);
all_Thetas(analyzed_ind(1:N-1),:) = Theta_estimates(:, 1:iter);
final_Theta = all_Thetas(:, iter);
all_Gammas = zeros(n-1, iter);
all_Gammmas(analyzed_ind(1:N-1),:) = Gamma_estimates(:, 1:iter);

% ## ILLUSTRATE FINAL ANALYSIS RESULTS ##
Channels_Data = zeros(n-1, 2);
Channels_Data(analyzed_ind(1:N-1), 1) = Plus_Channel_Freq_Data(1:N-1)/sum(Plus_Channel_Freq_Data(1:N-1));
Channels_Data(:, 2) = final_Theta;
figure;
bar(Channels_Data, 'group');
set(gca, 'xlim', [1 n-1]);
legend('Input Distribution', 'Optimized Distribution');
colormap Winter
title('Input and Output Distributions');

figure;
with_negatives_Theta = zeros(n-1,1);
with_negatives_Theta(analyzed_ind(1:N-1)) = prelim_Theta;
bar([Channels_Data(:, 1) with_negatives_Theta final_Theta], 'group');
set(gca, 'xlim', [1 n-1]);
legend('Input Distribution', 'Algorithm 1''s Output', 'Final Output');
title('Input, Intermediate, and Final Distributions');
colormap Jet

figure;
normalized_Theta = zeros(n-1,1);
normalized_Theta(analyzed_ind(1:N-1)) = recovered_Theta;
bar([with_negatives_Theta normalized_Theta final_Theta], 'group');
set(gca, 'xlim', [1 n-1]);
legend('Initial Results (with negatives)', 'Naively Corrected Distribution', 'Numerically Corrected Distribution');
colormap Jet

% ## OUTPUT TO EXCEL ## 
% Flip back to 5' --> 3' order
Nucleotide_IDs = counts_data(1:n, 1);
Thetas_for_Excel = flipud(final_Theta);
Output_for_Excel = [Nucleotide_IDs [-1; Thetas_for_Excel]];
all_c = [-111, c_estimates];
padding = -1*ones(1, iter);
Thetas_for_Log = [Nucleotide_IDs [padding; flipud(all_Thetas)]];
Gammas_for_Log = [Nucleotide_IDs [padding; flipud(all_Gammas)]];

save(strcat(basedir, 'Results.txt'), 'Output_for_Excel', '-ascii', '-double', '-tabs');
save(strcat(basedir, 'LogFile.txt'), 'Thetas_for_Log', 'Gammas_for_Log', 'all_c', '-ascii', '-double', '-tabs');


