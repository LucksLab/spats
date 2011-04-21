% This function copmutes H(c) for a given c
function [fval, Gradval] = H_fun(c)

global X_counts
global current_c_weight  
global current_elongation_rates 
global current_Theta

weighted_exp = exp(c*current_Theta);
deltas = weighted_exp - current_elongation_rates;

fval = -current_c_weight*c - X_counts*log(deltas);

numerators = current_Theta.*weighted_exp;

Gradval = -current_c_weight - X_counts*(numerators./deltas);