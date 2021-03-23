function bic2 = bicoh(data,fy,fx,pad);

% BIC2 The bicoherence is close to 0 when the coherence of Fourier
% frequency triplets (v, f, and v + f) behave randomly across S segments
% and it is close to 1 when the phase of the numerator (bispectrum) remains
% constant over segments; see Kim et al., (1979).
%
% Use as
%   bic2 = bicoh(dataw,fy,fx,pad);
%
% where
%   data      =  raw signal in a FieldTrip raw structure style. See ft_datatype_raw.m
%   fy        =  vector of frequencies in spectrum (y-axes comodulogram)
%   fx        =  [begin end], frequency band of interest (x-axes comodulogram)
%   pad       =  number, zero pad in seconds. Determines the frequency
%              resolution of the x-axes of the comodulogram.
%
% Output:
%   bic2.powspctrm   =  bicoherence spectrum
%
% References:
%
%   Kim, Y. C., and Powers, E. J. (1979). Digital bispectral analysis and its
%     applications to nonlinear wave interactions. IEEE Trans. Plasma Sci. 7,
%     120–131. doi: 10.1109/tps.1979.4317207
%
% Copyright (C) 2016, Diego Lozano-Soldevilla 04-Oct-2016 23:19:34
%
% License
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see http://www.gnu.org/licenses/.

if isfield(data,'fsample');
  fsample = data.fsample;
else
  fsample = 1./mean(diff(data.time{1}));
end

nchan = size(data.label,1);
nwind = unique(cellfun('length',data.trial));
ntrials = size(data.trial,2);
nfft = pad * fsample;

if size(nwind,2)>1;
  error('trials with different samples: not supported');
end

if nfft < nwind;
  error('the padding that you specified is shorter than the data');
end

% define frequency x-axes of the comodulogram based on the size of the
% sliding window 'nfft'. Minimum nfft will be same size as nwind
freqfft = (fsample/(nfft))*(0:((nfft)-1)/2);

% select only the frequencies of your interest
fxw =  nearest(freqfft,fx);
fxind = fxw(1):fxw(2);
fx = freqfft(fxind);
fyind = find(ismember(freqfft,fy)==1);

bic2 = [];
bic2.powspctrm = zeros(nchan,length(fyind),length(fxind));

dat = raw2dat(data);

for c=1:nchan;
  d = squeeze(dat(c,:,:))';
  for k=1:ntrials;
    signal = d(k,:);
    signal = signal - mean(signal);
    X = fft(hann(size(signal,2)).*signal',nfft);  % Take the FFT of segment with Hanning.
    
    for i=1:length(fxind);    % For each freq pair, compute numerator.
      for j=1:length(fyind);
        Bsp(j,i,k) = X(fxind(i))*X(fyind(j))*conj(X(fxind(i)+fyind(j)));
        de1(j,i,k) = abs(X(fxind(i))*X(fyind(j)))^2;
        de2(j,i,k) = abs(X(fxind(i)+fyind(j)))^2;
      end
    end
  end
  
  Bsp = mean(Bsp,3);
  de1 = mean(de1,3);
  de2 = mean(de2,3);
  bic2.powspctrm(c,:,:) = abs(Bsp).^2./(de1.*de2);
end

bic2.label     = data.label;
bic2.dimord    = 'chan_freq_freq';
bic2.fx        = fx;
bic2.fy        = fy;

if isfield(data,'grad');
  bic2.grad = data.grad;
end
if isfield(data,'elec');
  bic2.elec = data.elec;
end


function dat = raw2dat(data)
% RAW2DAT transforms the fieldtrip raw structure trial into dat matrix
%
% Diego Lozano-Soldevilla CERCO 21-Aug-2015 18:24:30

nchan = size(data.label,1);
nsmps = unique(cellfun(@length,data.trial));
if length(nsmps)>1;
  error('trials contain different amount of samples');
end
ntrls = size(data.trial,2);

dat = zeros(nchan,nsmps,ntrls);
for k=1:ntrls;
  dat(:,:,k) = data.trial{1,k};
end
