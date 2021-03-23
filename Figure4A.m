%% Figure 4A
% Simulated data with spurious phase-amplitude coupling using 3 harmonics
% with no frequency fluctuations with phase difference = (j-1) * pi/2
%
% Lozano-Soldevilla, D., ter Huurne, N., and Oostenveld, R. (2016).
% Neuronal Oscillations with Non-sinusoidal Morphology Produce Spurious
% Phase-to-Amplitude Coupling and Directionality. Front. Comput. Neurosci.
% 10.
%
% The code below uses the following wave equation to generate a low
% frequency signal with harmonics:
%
% signal = (1/j^2)*cos(j*f*2*pi*t + (j-1) * phi; See formula 10 in Lozano-Soldevilla et al., 2016
%
% where j index the n-harmonic of frequency f weighted a specific phase
% phi. (j-1) * phi determines the wave morphology:
%   - phi of (j-1) * pi/2 over the sum of  j-th harmonics generates a
%   left-sided sawtooth wave 
%   - phi of (j-1) * 3*pi/2 over the sum of j-th generates a 
%   right-sided sawtooth wave
%   - phi of (j-1) * 2*pi over the sum of j-th generates an
%   amplitude asymmetric waveform
%
% The simulations below used the fundamental and the firsts two harmonics
% (i.e. j = 3)
%
% The normalization factor 1/j^2 gives more weight the lower frequencies
% yielding the harmonic time series more curved shape that decrease the
% spurious cross-frequency spurious coupling due to sharp edge artifacts
% (Kramer et al 2008 J Neurosci Methods).
%
% For more details on harmonic wave generation see sections 3.3 and 3.4 in
% Maccarone, T. J. (2013). The biphase explained: understanding the
% asymmetries in coupled Fourier components of astronomical time series.
% Monthly Notices of the Royal Astronomical Society, 435(4), 3547-3558.
% doi: 10.1093/mnras/stt1546
%

dt = 0.001;
fsample = 1/dt;
t0 = 0.1;
s=[];

for i = 1:1800;
  
  f = 1/(t0);     % create low frequency signal with no period fluctuation (i.e. constant frequecy)
  t = 0:dt:1/f-dt;
  
  tmp=zeros(1,size(t,2));
  har=[];
  for j=1:3;
    har(j,:) = ((1/j^2)*cos(j*f*2*pi*t + (j-1) * pi/2)); % Maccarone 2013 MNRAS equation, see formula 10 in Lozano-Soldevilla et al., 2016
    tmp = tmp + har(j,:);
  end
  s = [s tmp];
end

t1 = (1:size(har,2))./fsample;
figure;
subplot(221);plot(t1,har');xlim([0 0.1]);hold all;
subplot(221);plot(t1,sum(har,1)','Color','k');xlim([0 0.1]);
legend('F0','H1','H2','total');
title('base function');
xlabel('Time (s)');
ylabel('Amplitude (a.u.)');
xlim([0 0.12]);

N = length(s);
t = (1:N)./fsample;
s = s + 0.1*randn(1,N); % add noise to the signal

% compute spectral power of the signal
df = 1.0 / (N*dt);
faxis = (0:N/2)*df;
pow = abs(fft(s.*hann(N)'-mean(s.*hann(N)')));
pow = pow(1:N/2+1);
pow = 10.0*log10(pow / max(pow));

subplot(223);plot(faxis,pow);xlim([2 120]);
title('power spectrum');
xlabel('Frequency (Hz)');
ylabel('dB');


%% create FieldTrip raw structure
data = [];
data.label{1}   = 'hp';
data.fsample    = fsample;
data.trial{1}   = s; clear s;
data.time{1}    = (1:size(data.trial{1},2))./data.fsample;
data.trialinfo  = [1];

cfg         = [];
cfg.length  = 2;
cfg.overlap = 0.5;
dataw = ft_redefinetrial(cfg,data); % cut data into epochs


%%
fy = 5:1:120;
fx = [0 30];
width = 5;
pad = 2;
freq = cfcoh(dataw,fy,fx,width,pad);

figure;
subplot(121);
imagesc(freq.fx,freq.fy,freq.COH,[0 1]);
title('Coherence');colorbar;axis xy;
set(gcf, 'Renderer', 'painters');
xlabel('Frequency phase (Hz)');
ylabel('Frequency amplitude (Hz)');

%%
fwidth = 4;
jack   = 1;
fx     = freq.fx(6:55);
psi    = cfd(freq,fx,fwidth,jack);

figure;
subplot(121);
imagesc(psi.time,psi.freq,squeeze(psi.PS./psi.PSstd),[-2.5 2.5]);
title('CFD');colorbar;axis xy;colormap('bluewhitered');
set(gcf, 'Renderer', 'painters');
xlabel('Frequency phase (Hz)');
ylabel('Frequency amplitude (Hz)');

%% bicoherence  
fy = 5:1:120;
fx = [5 30];
pad = 2;
bic2 = bicoh(dataw,fy,fx,pad);

figure;
subplot(121);imagesc(bic2.fx, bic2.fy, squeeze(bic2.powspctrm),[0 1]);
xlim([5 30]);
title('Bicoherence');colorbar;axis xy;
xlabel('Frequency (Hz)');
ylabel('Frequency (Hz)');
set(gcf, 'Renderer', 'painters');

