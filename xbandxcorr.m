function Y = xbandxcorr(X,D)
% Y = xbandxcorr(X,D)
%     Calculate cross-band correlations of signal envelopes.
%     X is a set of N signals, one per row.
%     D is a list of "steps"
%     Y returns N x len(D) - sum(D) values, as the normalized correlations
%     between each band, and the band D(i) steps higher.
%     D defaults to 1:(N-1)
% 2010-10-06 Dan Ellis dpwe@ee.columbia.edu

[N, nc] = size(X);

if nargin < 2;  D = 1:(N-1); end

lD = length(D);

nop = lD*N - sum(D);

Y = zeros(1,nop);
yp = 0;

% prenormalize rows so correlations are just inner products
% remove means
MX = mean(X,2);
NX = X - repmat(MX,1,nc);
% normalize energies
VX = sqrt(sum(NX.^2,2));
NVX = diag(1./VX)*NX;

for d = D
  for i = 1:N-d
    yp = yp + 1;
    Y(yp) = correlation(NVX(i,:),NVX(i+d,:));
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = correlation(X,Y)

C = sum(X.*Y);