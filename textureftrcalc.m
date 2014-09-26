function [F,E,X,frqs] = textureftrcalc(d,sr,thresh_db)
% F = textureftrcalc(d,sr,thresh_db)
%   Calculate "texture features" for input soundfile.
%   Input data is waveform d (at samplerate sr)
%   OR file named by string d
%   OR list of files in cell array d
%   Output F is a 200xN array of features (or cell array of arrays
%   if input is cell array).
%   thresh_db specifies threshold in dB below peak for silence
%   trimming (default 40.0).
% 2010-08-10 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2;  sr = 0;  end  % just to avoid problems passing down
if nargin < 3;  thresh_db = 40; end

if iscell(d)
  nd = length(d);
  for i = 1:nd
    F{i} = textureftrcalc(d{i},sr);
  end
  return
end

if isstr(d)
  if strcmp(d(end-3:end),'.mp3')
    [d,sr] = mp3read(d,0,1,2);  % always downsample by 2?
  else
    [d,sr] = wavread(d);
    if size(d,2) > 1
      d = mean(d,2);
    end
  end
end

% d is now mono waveform data sampled at SR
% resample to 8 kHz
sro = 8000;
if (sr ~= sro)
  gg = gcd(sro,sr);
  d = resample(d,sro/gg,sr/gg);
  sr = sro;
end

%%%% PLAN:
% 1. trim off any silence at start/end
% 2. Apply AGC
% 3. take auditory spectrogram
% 4. break into 10s windows with 5s hop
% 5. Calculate mean, variance, skew, kurtosis in each band
% 6. Calculate spectral flatness mean & variance in each band
% 7. Calculate modulation energy proportions in octaves for each band

% 1. trim off any silence at start/end
win_ms = 32.0;
[initialsil,finalsil] = findnonsilent(d,sr, win_ms, thresh_db);
d = d(initialsil:finalsil);

% 2. Apply AGC
[d,D,E] = tf_agc(d,sr);
% D is the STFT (on 32ms time base) of the original signal, 
% and E is the smoothed energy envelope.  Thus the STFT of 
% the AGC'd signal is (approximately)
D = D./E;

% 3. Take auditory spectrogram
nfft = 2*(size(D,1)-1);  % fft from tf_agc
[f2b,frqs] = fft2melmx(nfft,sr);
f2b = f2b(:,1:size(D,1));
DA = f2b*abs(D);

frametime = nfft/2/sr; % hop time from tf_agc

% 4. Break into ~10s windows with 5s hop
envhoptime = 4.0;
envwindurtarg = 8.0;
envhop = round(envhoptime/frametime);
envwin = 2^ceil(log(envwindurtarg/frametime)/log(2)); % enclosing power of 2
envwindur = envwin * frametime;

% How many of these envwin frames can we fit into DA with envhop hop?
nenvfr = 1+floor((size(DA,2)-envwin)/envhop);

disp(['trimmed dur = ',num2str(length(d)/sr),' s, ',num2str(nenvfr),' frames']);

nbands = size(DA,1);

nftrs = 10*nbands;
F = zeros(nenvfr,nftrs);

for fr = 1:nenvfr
  
  envfr = 20+20*log10(DA(:,(fr-1)*envhop + [1:envwin]));
  
  % 5. Calculate mean, variance, skew, kurtosis in each band

  % Map to dB, clamp low end at 40 dB below peak
  eceil = percentile(envfr(:),0.95);
  efloor = eceil - 40.0;
  envfr(envfr(:) < efloor) = NaN;
  
  %F(fr,1:nbands) = mean(envfr,2);
  %F(fr,nbands + [1:nbands]) = std(envfr,0,2)./F(fr,1:nbands)';
  MVSK = rowmeanstdskewkurt(envfr); % ignores NaN values
  F(fr,[1:nbands]) = MVSK(:,1);
  F(fr,nbands+[1:nbands]) = MVSK(:,2);
  F(fr,2*nbands+[1:nbands]) = MVSK(:,3);
  F(fr,3*nbands+[1:nbands]) = MVSK(:,4);
  
  % 6. Calculate spectral flatness mean & variance in each band

  
  % 7. Calculate modulation energy proportions in octaves for each band

  % replace NaNs with floor value
  envfr(isnan(envfr(:))) = efloor;
  MSE = envmodspec(envfr);
  % last col is 16-31 Hz band (for 16 ms = 62.5 Hz frame SR)
  % so take last 6 cols to get down to 0.5-1 Hz band
  % e.g. 0.5-1 1-2 2-4  4-8 8-16 16-31
  F(fr,4*nbands+[1:nbands]) = MSE(:,end-5);
  F(fr,5*nbands+[1:nbands]) = MSE(:,end-4);
  F(fr,6*nbands+[1:nbands]) = MSE(:,end-3);
  F(fr,7*nbands+[1:nbands]) = MSE(:,end-2);
  F(fr,8*nbands+[1:nbands]) = MSE(:,end-1);
  F(fr,9*nbands+[1:nbands]) = MSE(:,end-0);

  % 8. Cross-band correlations
%  lags = [3 6 9]; % 36 elements
  lags = [3 6 9 1 2 4 5 7 8 10 11 12];  % 138 elements
  F(fr,10*nbands+[1:(nbands*length(lags)-sum(lags))]) = xbandxcorr(envfr,lags);
  
  if fr == 1
    % debug output
    X = envfr;
  end
  
end

E = DA;

% debug output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = rowmeanstdskewkurt(X)
% Calculate row-wise means, stddevs, skews, kurts, but ignore NaN vals
[nr,nc] = size(X);
Y = zeros(nr,4);
for i = 1:nr
  dd = X(i,~isnan(X(i,:)));
  Y(i,1) = mean(dd,2);
  Y(i,2) = std(dd,0,2);
  Y(i,3) = skewness(dd,1,2);
  Y(i,4) = kurtosis(dd,1,2);
end  
