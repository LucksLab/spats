% This function computes the water-filling function values for a given 
% vector of break-point values.
function [fval, sorted_ind] = Compute_WF(break_points, elong_rates, num_weights, den_weights)

global current_c
global X_counts
global X_0

% We sort all break points in descending order, and then evaluate
% them top-to-bottom, that is, we gradually decrease \nu*, and at each 
% break point, an additional factor is added to the total product.  
[sorted_bp, sorted_ind] = sort(break_points, 'descend');
% We are interested in evaluating the function only at the 
% breakpoints that lie above c*sum(X_k's).
[lowest_ind] = find(sorted_bp < current_c*(sum(X_counts)+X_0), 1, 'first');
if isempty(lowest_ind)
    % When all breakpoints are above the total sum of counts 
    % (including X_0), then the optimal solution is obtained exactly 
    % at this total sum. Hence, it is redundant to perform water 
    % filling, but for convenience, we currently allow the software 
    % to look for this point. 
    % TO BE IMPROVED IN FUTURE VERSIONS.
    ordered_eval_points = [sorted_bp(2:end); current_c*(sum(X_counts))];
else
    % We added the total sum as a last point to avoid numerical precision
    % issues.
    ordered_eval_points = [sorted_bp(2:lowest_ind); current_c*(sum(X_counts))];
end
eval_size = length(ordered_eval_points);
sorted_eval_ind = sorted_ind(1:eval_size);
% This matrix contains the relevant ordered break points column 
% vector in each of its columns.
sorted_bp_matrix = diag(ordered_eval_points, 0)*ones(eval_size,eval_size);

num_matrix = (sorted_bp_matrix - num_weights(1:eval_size,sorted_eval_ind)).*elong_rates(1:eval_size,sorted_eval_ind);
num_matrix = tril(num_matrix,0.0) + triu(ones(eval_size,eval_size),1.0);

den_matrix = sorted_bp_matrix - den_weights(1:eval_size,sorted_eval_ind);
den_matrix = tril(den_matrix,0.0) + triu(ones(eval_size,eval_size),1.0);

quotient_matrix = num_matrix./den_matrix;

% NOTE: all the matrix entries, except for those in the last line, 
% should be equal to or greater than 1. This is because all evaluation 
% points except for the bottom one are not below c*sum(X_k's), which 
% implies that each of the individual terms we compute must correspond 
% to e^(c*theta) > 1. In other words, we are guaranteed not to "cross" 
% the "limit point" on the theta axis towards the "undefined" domain. 
% However, this is not true for the bottom point, where some terms 
% might correspond to the "undefined" domain of their corresponding 
% fuctions. In such cases, we will have a value smaller than 1 
% (it can even be negative), and in order not to decrease the total 
% product, we assign e^c to exp(c*theta_k), since this is the maximum 
% value that this expression can assume under the constraints 0<=theta_k<=1. 
quotient_matrix(eval_size, find(quotient_matrix(eval_size,:) < 1)) = exp(current_c);

% Due to numerical precision issues, some entries are slightly 
% smaller than 1, although they pertain to valid expressions 
% (i.e., inside the function's domain). This happens when the 
% values at the numerator and denominator are extremely close, 
% but still, the numerator is slightly smaller. This issue 
% arises, for example, when consecutive alpha's are 
% theoretically identical, but numerically they are a little different. 
% In such cases, the numerator might be smaller than the denominator. 
% To such cases, we set those entries to 1.  
quotient_matrix(find(quotient_matrix < 1.0)) = 1.0;

fval = [prod(quotient_matrix, 2) ordered_eval_points];
