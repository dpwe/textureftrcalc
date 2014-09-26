function F = specflatness(X,W)
% F = specflatness(X,W)
%   Calculate spectral flatness of columns of X, using bands
%   defined by rows of W.
% 2010-08-13 Dan Ellis dpwe@ee.columbia.edu

amean = W*X;
gmean = exp(W*log(X));

F = gmean./amean;
