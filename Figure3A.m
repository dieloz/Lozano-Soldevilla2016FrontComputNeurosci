%% Simulation Figure 3A
% Simulated data with true phase coupling with sinusoids
%
% Lozano-Soldevilla, D., ter Huurne, N., and Oostenveld, R. (2016).
% Neuronal Oscillations with Non-sinusoidal Morphology Produce Spurious
% Phase-to-Amplitude Coupling and Directionality. Front. Comput. Neurosci.
% 10.
%
% The code below is an adpatation from (exam3): Kramer, M. A., Tort, A. B.
% L., & Kopell, N. J. (2008). Sharp edge artifacts and spurious coupling in
% EEG frequency comodulation measures. Journal of Neuroscience Methods,
% 170(2), 352-357. doi: 10.1016/j.jneumeth.2008.01.020

dt = 0.001;
fsample = 1/dt;
t0 = 0.1; % mean period of alpha oscillations -> 1/t0 = 10 Hz

s = [];

for k=1:1730;
    f = 1/(t0+rand()*0.015);              % Create alpha oscillations fluctuating in period
    s1 = cos(2.0*pi*(0:dt*f:1-dt*f));
    
    good = find(s1 > 0.999);              % detect specific low frequency phase
    s2 = zeros(1,length(s1));
 
    stemp = rand(1,3000);
    stemp = ft_preproc_bandpassfilter(stemp, 1000, [60, 80]); % Bandpass filter the gamma band
    stemp = stemp(2000:2049);             % Duration 50 ms.
    stemp = 5.0*hanning(50)'.*stemp;      % Hanning tapered to increase gamma amplitude
    
    rindex = ceil(rand()*10);             
    s2(rindex+good(1):rindex+good(1)+50-1)=stemp; % Add the gamma at specific alpha phase
    s = [s, s1+s2];
end

t1 = (1:size(s1,2))./fsample;

N = length(s);
t = (1:N)./fsample;
s = s + 0.1*randn(1,N); % add noise to the signal

figure;
subplot(221);plot(t,s);xlim([5.9 7.1]);
title('raw signal');
xlabel('Time (s)');
ylabel('Amplitude (a.u.)');

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
