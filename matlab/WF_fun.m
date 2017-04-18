% This function computes the water-filling function's value at a given \nu.
function fval = WF_fun(nu)

global current_c 
global current_contributing_ind
global current_elongation_rates
global current_num_weights_vec
global current_den_weights_vec

% We order the relevant constants and then subtract them from nu in both 
% the numerator and the denominator. 
quotient_vec = (nu - current_num_weights_vec(current_contributing_ind))./(nu - current_den_weights_vec(current_contributing_ind));
factor_vec = quotient_vec.*current_elongation_rates(current_contributing_ind);

fval = prod(factor_vec) - exp(current_c);
