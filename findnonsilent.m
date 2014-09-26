function [startsamp,endsamp] = findnonsilent(d,sr,win_ms,thresh_db)
% [startsamp,endsamp] = findnonsilent(d,sr,win_ms,thresh_db)
%   Takes waveform in d (at sr) and breaks it into nonoverlapping 
%   blocks of win_ms (32 ms), finds the largest block, then returns 
%   the sample indices of the first and last time blocks exceed 
%   thresh_db (40 dB) below the peak.  Thus d(startsamp:endsamp)
%   truncates d to exclude initial and final near-silent regions.
% 2010-08-10 Dan Ellis dpwe@ee.columbia.edu 

if nargin < 3;  win_ms = 32;  end
if nargin < 4;  thresh_db = 40; end

% Figure energy in 32 ms blocks (nonoverlapping, rect win)
eblklen = round(win_ms/1000 * sr);
nblk = floor(length(d)/eblklen);
ed = 10*log10(sum(reshape(d(1:(nblk*eblklen)).^2,eblklen,nblk)));
maxed = max(ed);
firstnz = min(find(ed > (maxed - thresh_db)));
lastnz = max(find(ed > (maxed - thresh_db)));

startsamp = (firstnz-1)*eblklen + 1;
endsamp = lastnz*eblklen;

% special case: if nothing is trimmed from end, return every last
% sample (don't truncate on eblklen blocks)
if lastnz == length(ed)
  endsamp = length(d);
end

