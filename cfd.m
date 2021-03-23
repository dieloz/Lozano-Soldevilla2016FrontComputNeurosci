function CFD = cfd(freq,fx,fwidth,jack)
% CFD computes the phase-slope index between the phase of slower
% oscillations and the power envelope of faster oscillations across trials;
% see Jiang et al., (2015).
%
% Use as
%   CFD = cfd(freq,fx,fwidth,jack)
%
% where
%   freq        =  frequency data type in a FieldTrip raw structure style. See ft_datatype_freq.m
%   fx          =  vector of frequencies in spectrum x-axes comodulogram).
%   fwidth      =  frequency bandwidth to estimate the phase slope [f-fwidth/2,f+fwidth/2]
%   jack        =  jacknife estimate of the phase slope.
%
% Output:
%   CFD.PS      = cross-frequency directionality index
%   CFD.PSstd   = if jack==1, standard error of the cross-frequency phase-slope
%   CFD.time    = fx;
%   CFD.freq    = fy;

%
% References:
%
%   Jiang, H., Bahramisharif, A., van Gerven, J.M.A., and Jensen, O.
%     (2015). Measuring directionality between neuronal oscillations of
%     different frequencies. Neuroimage 118, 359-367.
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

% cross-spectral density matrix between signal spectrum (Sl) and high frequency
% envelope spectrum (SPl)

C = freq.Sl.*conj(freq.SPl);

ntrials = size(C,1);
nchan   = size(C,2);

fy     = freq.fy;
dealtf = freq.dfreq;

if jack ==0    
  CS  = mean(C,1)./sqrt((mean(abs(freq.Sl).^2,1) .* mean(abs(freq.SPl).^2,1)));
  dim = size(CS);
  CS  = reshape(CS,dim(2:4));
  PS  = phaseslope(CS,fy,fx,fwidth,dealtf);

else
  jackid = 1:ntrials;
  PSjack = zeros(nchan,numel(fy),numel(fx),ntrials);
  for c = 1:nchan;
    for i = 1:ntrials;
      jacktmp = jackid;
      jacktmp(i) = [];
      Sltmp   = freq.Sl(jacktmp,c,:,:);
      SPltmp  = freq.SPl(jacktmp,c,:,:);
      Ctmp    = C(jacktmp,c,:,:);
      CSjack   = mean(Ctmp,1)./sqrt((mean(abs(Sltmp).^2,1) .* mean(abs(SPltmp).^2,1)));
      dim = size(CSjack);
      CSjack = reshape(CSjack,dim(2:4));
      PSjack(c,:,:,i) = phaseslope(CSjack,fy,fx,fwidth,dealtf);
    end
  end
  
  PS    = squeeze(mean(PSjack,4));
  PSstd = squeeze(std(PSjack,0,4))*sqrt(ntrials);
end

CFD.PS = PS;
if (jack ==1)
  CFD.PSstd = PSstd;
end
CFD.label = freq.label;
CFD.time = fx;
CFD.freq = fy;
CFD.dimord = 'chan_freq_time';


function PS = phaseslope(CS,fy,fx,fwidth,dealtf)

nchan = size(CS,1);
PS = zeros(nchan,numel(fy),numel(fx));

for c = 1:nchan;
  for k = 1:numel(fx);
    maxfreqbin = floor(fx(k)/dealtf);
    step       = round(0.5*fwidth/dealtf);
    for j =1:numel(fy);
      CStmp = squeeze(CS(c,j,:));
      PStmp = imag(conj(CStmp(1:end-1)).*CStmp(2:end));
      PS(c,j,k)= mean(PStmp(maxfreqbin-step:maxfreqbin+step));
    end
  end
end
