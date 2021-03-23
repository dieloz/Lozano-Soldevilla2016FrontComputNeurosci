function freq = cfcoh(data,fy,fx,width,pad)

% CFCOH performs coherence between a signal and the time-course of
% the power at higher frequencies across trials; see Osipova et al.,(2008).
%
% Use as
%   freq = cfcoh(data,fy,fx,width,pad)
%
% where
%   data      =  raw signal in a FieldTrip raw structure style. See ft_datatype_raw.m
%   fy        =  vector of frequencies in spectrum (y-axes comodulogram)
%   fx        =  [begin end], frequency band of interest (x-axes comodulogram)
%   width     =  number, width of the wavelet, determines the temporal and
%              spectral resolution of the x-axes of the comodulogram.
%   pad       =  number, zero pad in seconds. Determines the frequency
%              resolution of the x-axes of the comodulogram.
%
% Output:
%   freq.Sl   =  FFT of the raw signal
%   freq.SPl  =  FFT of high frequency power
%   freq.C    =  Cross-spectra between raw signal and high frequency power
%
% References: 
% 
%   Osipova, D., Hermes, D., and Jensen, O. (2008). Gamma power is
%     phase-locked to posterior alpha activity. PLoS.One. 3, e3990.
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

i = sqrt(-1);
nchan = size(data.label,1);
fsample = data.fsample;
ntrials = length(data.trial);
nwind = unique(cellfun('length',data.trial));
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
dfreq   = fsample/nfft;

% select only the frequencies of your interest
fxind =  nearest(freqfft,fx);
fx = freqfft(fxind(1):fxind(2));

% hanning tapers applied for raw signal when calculating the phase consistency
tapdat = hanning(nwind);

C   = zeros(ntrials,nchan,length(fy),length(fx));
Sl  = zeros(ntrials,nchan,length(fy),length(fx));
SPl = zeros(ntrials,nchan,length(fy),length(fx));

%% COHERENCE COMPUTATION
% Calculate the temporal development of the high-frequency power using a
% sliding time-window N points long after applying a Hanning taper. Do this
% when looping over the frequencies on the y-axis.
for c = 1:nchan;
  for l = 1:ntrials;
    for k = 1:length(fy);
      
      f = fy(k);
      n = floor(width*fsample/f);
      tapy = hanning(n)';
      
      % amplitude envelope of high frequencies: formula (1) in Lozano-Soldevilla et al 2016 Front Comput Neurosci
      sP = abs(conv(data.trial{1,l}(c,:),tapy.*exp(i*2*pi*f.*(1:n)/fsample),'same'));
      sl  = fft(tapdat'.*data.trial{1,l}(c,:),nfft);  % formula (2)
      spl = fft(tapdat'.*sP,nfft);                    % formula (3)
      ccc = sl.*conj(spl);                            % cross spectral density matrix
      Sl(l,c,k,:)  = sl(fxind(1):fxind(2));  clear sl;
      SPl(l,c,k,:) = spl(fxind(1):fxind(2)); clear spl;
      C(l,c,k,:)   = ccc(fxind(1):fxind(2)); clear ccc; 
    end
  end
end

freq.label = data.label;
freq.Sl  = Sl;  clear Sl;
freq.SPl = SPl; clear SPl;
freq.C   = C;   clear C;
freq.c_dimord = 'rpt_chan_freq_freq';
freq.COH = squeeze(abs(mean(freq.C,1)).^2 ./(mean(abs(freq.Sl).^2,1) .* mean(abs(freq.SPl).^2,1))); % formula (4)
freq.coh_dimord = 'chan_freq_freq';
% coherence = abs(mean(freq.C,3))./sqrt( (mean(abs(freq.Sl).^2,3) .* mean(abs(freq.SPl).^2,3)) );

freq.fx      = fx;
freq.fy      = fy;
freq.fsample = fsample;
freq.dfreq   = dfreq;
