%% initialization

impulseName = 'wav/semMidCenterImp.wav';
respName = 'wav/semMidCenterResp.wav';

scaleAmplitude = 2;

% mode analysis
nbins = 2048;   % stft analysis half bandwidth, bins
nskip = 32;    % stft hop size, samples
beta = 3;   % stft peak half width, bins
offset = 0.0;   % response start time delay, frames
rmin = 0.05;    % response end time level, ratio
tmin = 0.0; % minimum response time, seconds

%% compute STFTs and plot the response spectrogram

% load impulse and response
[impulse, fs] = wavread(impulseName);
[response, fs] = wavread(respName);

impulseSTFT = stft(sum(impulse,2),nbins,nskip);
figure(1);
responseSTFT = ftgram(sum(response,2), fs, 'music', 'nbins', nbins, 'nskip', nskip);

%% computer impulse response and plot it

impulseSpectrum = mean(abs(impulseSTFT),2)/max(mean(abs(impulseSTFT),2));
responseSpectrum = mean(abs(responseSTFT),2)/max(mean(abs(responseSTFT),2));
irSpectrum = (responseSpectrum./impulseSpectrum).*scaleAmplitude;

% deconvolution
irSTFT = responseSTFT./impulseSTFT;
% compute the spectrum
%irSpectrum = mean(abs(irSTFT),2)/max(mean(abs(irSTFT),2));

% define time, frequency axes
[~, nframes] = size(irSTFT);
f = fs/2*[0:nbins]'/nbins;
t = [0.5:nframes-0.5]*nskip/fs;

% figure(2);
% plot(f, 20*log10(irSpectrum), '-'); grid;
% title('Spectrum');
% xlabel('frequency, Hz'); ylabel('power, dB');
% xlim([20 10000]);

%% find mode frequencies

[gammam, im] = localmax(20*log10(abs(max(irSpectrum, 2))));
fm = (im-1)/nbins*fs/2;

nmode = length(fm);

index = [1 2 4 5 6 7 8 9 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 28 30 31 32 33 34 36 37 39 40 41 42 ];
nmode = length(index);
fm = fm(index);
gammam = gammam(index);

figure(3);
plot(f, 20*log10(abs(irSpectrum)), '-', fm, gammam, 'o'); grid;
title('bowl spectrum (-), mode frequencies (o)');
xlabel('frequency, Hz'); ylabel('power, dB');
xlim([20 10000]);

%% estimate T60 and amplitudes

delta = round(offset*fs/nskip);
gammam = zeros(nmode,1);
rt60m = zeros(nmode,1);
gm = zeros(nmode,1);
for m = [1:nmode],
    % extract mode response
    psdm = mean(abs(responseSTFT(round((nbins-1)*fm(m)*2/fs)+[-beta:beta],:)),1)';

    % estimate decay time, initial amplitude
    [gmax, fmax] = max(psdm);
    fstart = fmax + delta;
    fstop = max(find(psdm > rmin*gmax, 1, 'last'), fstart+round(tmin*fs/nskip));
    index = [fstart+1:fstop];

    basis = [ones(fstop-fstart,1) t(index)'];
    theta = basis\(20*log10(psdm(index)));

    gm(m) = 10.^((theta(1)+theta(2)*t(fmax+1))/20);
    rt60m(m) = -60/theta(2);

    % plot mode response
    % figure(4);
    % plot(t, 20*log10(psdm), '-', t(index), basis*theta, '-', 'LineWidth', 2); grid;
    % title(['response energy profile, mode ', int2str(m)]);
    % xlabel('time, seconds'); ylabel('power, dB');

    % pause(0.5);

end;

% normalize mode amplitude
gm = gm/max(gm);

% plot mode amplitudes, decay times
figure(4);
subplot(2,1,1); plot(fm, 20*log10(gm), '-o'); grid; xlim([500 11000]);
title('Mode amplitudes');
xlabel('mode frequency, Hz'); ylabel('amplitude, dB');

subplot(2,1,2); plot(fm, rt60m, '-o'); grid; xlim([500 11000]);
title('Mode T60s');
xlabel('mode frequency, Hz'); ylabel('60 dB decay time, seconds');








