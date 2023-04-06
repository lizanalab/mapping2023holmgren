function [B,twom] = fractalglobulepoly_f(A, D, gamma, alpha)
%MODULARITY returns fractal globule polymer modularity matrix for network given
% by adjacency matrix A and distance matrix D.
%
% Version: 2.2.0
% Date: Mon 26 Aug 2019
%
%Works for directed and undirected networks
%
%   Input: A:  NxN adjacency matrices of a undirected network
%          D:  NxN distance matrix defined for indicies i, j as |i-j|,
%           explicitly defined as certain regions have uncoded sequences
%          gamma: resolution parameter
%          alpha: decay exponent for fractal globule polymer null-model
%
%   Output: B: function handle where B(i) returns the ith column of
%          [N]x[N] modularity matrix of the monolayer network
%           with adjacency matrix A
%           twom: normalisation constant
%
%   Example of usage: [B, twom]=modularity(A, .7, 1);
%          [S,Q]= genlouvain(B);
%          Q=Q/twom;
%
%   Notes:
%     The matrix A is assumed to be square. This assumption is not checked
%     here.
%
%     By using this code, the user implicitly acknowledges that the authors
%     accept no liability associated with that use.  (What are you doing
%     with it anyway that might cause there to be a potential liability?!?)
%
%   References:
%     Newman, Mark E. J. and Michelle Girvan. "Finding and Evaluating
%     Community Structure in Networks", Physical Review E 69, 026113 (2004).
%
%     Lee, Sang Hoon et. al. "Mapping the spectrum of 3D communities in human
%     chromosome conformation capture data", Scientific Reports 9, 6859 (2019).


if nargin < 2 || isempty(gamma)
    gamma = 1;
    alpha = 1;
elseif nargin < 3 || isempty(alpha)
    alpha = 1;
end

k = full(sum(A));
twom = sum(k);

Ds = ((double(D)+eye(size(D))).^(-alpha)-eye(size(D))).*(k'*k);
kD = full(sum(Ds));
sumD = sum(kD);

B = @(i) A(:,i) - gamma*Ds(:,i)/(sumD/twom);

end
