function S = effective_path_count(energy, threshold)
% =========================================================================
% Program : effective_path_count.m
% Description :
%   "有效路徑數 S" (chap4_list.md Sim-15/Sim-16): the minimum number of
%   LARGEST entries of an energy vector whose cumulative sum reaches a
%   given fraction of the total energy. Used two ways in Stage 3:
%     - Sim-15: energy = sinc(l - l_i)^2 across candidate delay taps l
%       -> quantifies how many delay bins share one fractional path's energy.
%     - Sim-16: energy = |H(row,:)|^2 across one row of the DD-domain
%       channel matrix -> quantifies per-symbol detection complexity.
%
% Input  : energy    - nonnegative vector (any shape)
%          threshold - fraction of total energy to capture (default 0.99)
% Output : S         - effective count (integer, >= 0; 0 if energy is all zero)
% =========================================================================

if nargin < 2 || isempty(threshold), threshold = 0.99; end

e = sort(energy(:), 'descend');
total = sum(e);
if total <= 0
    S = 0;
    return;
end

cs = cumsum(e) / total;
idx = find(cs >= threshold, 1, 'first');
if isempty(idx)
    S = numel(e);
else
    S = idx;
end

end
