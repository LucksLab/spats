% This function computes the MLE of the Theta vector when given raw 
% fragment counts from both the (-) and the (+) channels ("analyzed_data"). 
% The optimization flag ("optim_flag") options are: 
% 1 = set negatives to zero and normalize, 
% 2 = apply numerical optimization.
% The "max_iter" parameter specifies the maximum number of iterations to 
% be used for numerical optimization. 
function [Theta, Gamma, c, initial_c, delta, num_neg, p_0_hat, return_flag] = MLE_one_molecule(analyzed_data, optim_flag, max_iter, precision)

global X_counts
global X_0
global Y_counts
global current_num_weights_vec
global current_den_weights_vec 
global current_contributing_ind
global current_c 
global current_Theta
global current_elongation_rates
global current_c_weight 

N = size(analyzed_data, 1);

% ## COMPUTE MLE ##
% ## APPLY ALGORITHM 1 ##
% We first transform the counts into frequencies.
% Note that position "N" (last line) pertains to counts of complete
% fragments, and should be treated differently than all other counts!
Plus_Channel_Freq_Data = analyzed_data(1:N,3)/sum(analyzed_data(1:N,3));
X_0 = analyzed_data(N,3);
Minus_Channel_Freq_Data = analyzed_data(1:(N-1),4)/sum(analyzed_data(:,4));
p_0_hat = analyzed_data(N,4)/sum(analyzed_data(:,4));

% Compute the gammas.
current_dropoff_rates = zeros(N-1,1);
for j=1:N-1
    current_dropoff_rates(j) = analyzed_data(j,4)/sum(analyzed_data(j:N,4));
end

% consistency check
if ((X_0 == 0) || (p_0_hat == 0))
    disp('Zero complete fragments in at least one of the channels!');
    return_flag = -1;
    Theta = zeros(N-1,1);
    Gamma = zeros(N-1,1);
    c = 0;
    initial_c = 0;
    num_neg = 0;
    delta = 0;
    return
end

% We now recover the Poisson rate estimate (c_hat).
c_estimate = log(p_0_hat) - log(X_0/sum(analyzed_data(:,3)));

% consistency check
if (c_estimate < 0)
    disp('Estimated Poisson rate is negative!');
    return_flag = -2;
    Theta = zeros(N-1,1);
    Gamma = zeros(N-1,1);
    c = c_estimate;
    initial_c = c_estimate;
    num_neg = 0;
    delta = 0;
    return
end

recovered_Theta = zeros(N-1,1);

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
negative_entries = find(recovered_Theta < 0);
delta = 1 - sum(recovered_Theta(negative_entries));
num_neg = length(negative_entries);
recovered_Theta(negative_entries) = 0;
recovered_Theta = recovered_Theta/sum(recovered_Theta);
% takes care of numerical precision issues
if (sum(recovered_Theta) > 1)
    recovered_Theta = recovered_Theta/sum(recovered_Theta);
end

if (optim_flag == 1)
    Theta = recovered_Theta;
    Gamma = current_dropoff_rates;
    c = c_estimate;
    initial_c = c_estimate;
    return_flag = 1;
    return
else % optim_flag == 2
    % Numerical optimization. 
    % Split into 3 local convex optimization problems, and iteratively solve 
    % them until convergence is observed or a maximum number of iterations 
    % is reached.

    % Storage for the algorithm's progress reports
    Theta_estimates = zeros(N-1, max_iter+1);
    Gamma_estimates =  zeros(N-1, max_iter+1);
    c_estimates = zeros(1, max_iter+1);
    log_lik_vals = zeros(3, max_iter+1);

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

    % Set the values for interation 1.
    current_c = c_estimate;
    current_Theta = recovered_Theta;
    current_T_vec = exp(current_c*current_Theta) - 1;
    current_elongation_rates = 1 - current_dropoff_rates;
    iter = 1;
    relative_log_change = 1;

    % LOOP
    while ((iter <= max_iter) && (relative_log_change > precision))
        
        disp(sprintf('Starting iteration no. %d', iter));

        % Step 1: Solve P1 = maximize H(c) over R+. 
        % Here, P and Theta are kept fixed.
        current_c_weight = sum(current_Theta.*partial_X_count_sums) - sum(X_counts) - X_0;
        disp('Solving P1...');
        % minimize -H(c) without the non-negativity constraint
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

        disp(sprintf('ITERATION %d', iter));
        disp(sprintf('c estimate is %f', current_c));

        % Keep progress records
        c_estimates(iter+1) = current_c;
        Gamma_estimates(:,iter+1) = current_dropoff_rates;
        Theta_estimates(:,iter+1) = current_Theta;
        relative_log_change = -(log_lik_vals(3, iter+1) - log_lik_vals(3, iter))/log_lik_vals(3, iter);
        iter = iter + 1;
    end  % LOOP

        Theta = Theta_estimates(:,iter);
        Gamma = Gamma_estimates(:,iter);
        c = c_estimates(iter);
        initial_c = c_estimate;
        return_flag = 1;
end
