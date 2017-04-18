% This function computes the entropy of a given probability distribution D.
function fval = entropy(D)

% ************************************************************************
% *** NOTE: add security checks here - input must represent a prob.
% dist.!!!
% ************************************************************************
included_ind = D > 0;
sub_D = D(included_ind);
fval = -sum(sub_D.*log(sub_D));