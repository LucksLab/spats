% This function computes the log-likelihood value for the joint 
% likelihood function, which integrates data from both channels.
function fval = log_likelihood_fun(Theta_vals, Gamma_vals, c_val)

global X_counts 
global X_0
global Y_counts

theta_size = length(Theta_vals);

S_coeff = zeros(1, theta_size);
% Compute noise-based coefficients.
for k=1:theta_size
    S_coeff(k) = Y_counts(theta_size + 1) + X_0 + sum(X_counts(k+1:theta_size)) + sum(Y_counts(k+1:theta_size));
end

% Compute modification-based coefficients.
partial_X_count_sums = zeros(theta_size,1);
for k=2:theta_size
    partial_X_count_sums(k) = sum(X_counts(1:k-1));
end

c_coeff = sum(Theta_vals.*partial_X_count_sums) - sum(X_counts) - X_0;

exp_deltas = exp(c_val*Theta_vals) + Gamma_vals - 1;

positive_Y_ind = find(Y_counts(1:theta_size) > 0);  % we don't want Y_0 to be included in this index list
positive_X_ind = find(X_counts > 0);

% We do not include zero X and Y counts in the summation due to
% numerical issues.
fval = Y_counts(positive_Y_ind)*log(Gamma_vals(positive_Y_ind)) +  S_coeff*log(1-Gamma_vals) + c_val*c_coeff + X_counts(positive_X_ind)*log(exp_deltas(positive_X_ind));
