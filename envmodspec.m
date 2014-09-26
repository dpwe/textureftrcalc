function Y = envmodspec(X)
%  X is an envelope.  Take its FFT, then integrate energy in octave
%  bands. 
%  Works for X as an array of envelopes (sgram-style i.e. process
%  by row)
%  First output column correspnds to dc, then first 2 bins, then
%  2..4, right up to 2nd half, thus num output cols = 1+log_2(length(X)/2).
% 2010-08-13 Dan Ellis dpwe@ee.columbia.edu

nfft = size(X,2);

W = hann(nfft)';
WX = (X - repmat(mean(X,2),1,nfft))*diag(W);
F = abs(fft(WX,nfft,2));
F = F(:,1:nfft/2);
% sum in octaves, starting at 0.5 Hz, normalized by total energy
ly = 1;
y = 1;
E = sum(F,2);
F = diag(1./E)*F;
nocts = 1+round(log(nfft/2)/log(2));

for o = 1:nocts
  
%  [ly,y]
  Y(:,o) = sum(F(:,ly:y),2);
  ly = y;
  y = 2*y;

end
