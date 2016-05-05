%% Notes
% - Currently, the deconvolved IR is used to detect the frequency and the
%   amplitude of the modes and the response IR is used to detect the T60.

%% Global Parameters
% - whichTest (0-6): which test to run (see next section to see options) 
% - modeMatching (0-1): tests if modes are matching both ways
% - nmodes: the maximum number of modes to be considered
% - nmodes: the maximum number of modes to be considered for "the proof"

whichTest = 0;
modeMatching = 1;
nmodes = 30;
nmodesProof = 30;

%% STFT Parameters

nbins = 2048;   % stft analysis half bandwidth, bins
nskip = 32;     % stft hop size, samples
beta = 3;       % stft peak half width, bins
offset = 0.0;   % response start time delay, frames
rmin = 0.05;    % response end time level, ratio
tmin = 0.0;     % minimum response time, seconds

%% Init with squareBig -> squareSmall

if whichTest == 0 % squareBig -> squareSmall
    impulseName0 = 'wav/squareBigCenterImp.wav';
    respName0 = 'wav/squareBigCenterResp.wav';
    impulseName1 = 'wav/squareSmallCenterImp.wav';
    respName1 = 'wav/squareSmallCenterResp.wav';
    impulseNameProof = 'wav/squareMidCenterImp.wav';
    respNameProof = 'wav/squareMidCenterResp.wav';
    interp = 0.55; % TODO: not tested
elseif whichTest == 1 % semBig -> semSmall
    impulseName0 = 'wav/semBigCenterImp.wav';
    respName0 = 'wav/semBigCenterResp.wav';
    impulseName1 = 'wav/semSmallCenterImp.wav';
    respName1 = 'wav/semSmallCenterResp.wav';
    impulseNameProof = 'wav/semMidCenterImp.wav';
    respNameProof = 'wav/semMidCenterResp.wav';
    interp = 0.47; % TODO: not tested
elseif whichTest == 2 % roundBig -> roundSmall
    impulseName0 = 'wav/roundBigCenterImp.wav';
    respName0 = 'wav/roundBigCenterResp.wav';
    impulseName1 = 'wav/roundSmallCenterImp.wav';
    respName1 = 'wav/roundSmallCenterResp.wav';
    impulseNameProof = 'wav/roundMidCenterImp.wav';
    respNameProof = 'wav/roundMidCenterResp.wav';
    interp = 0.42; % TODO: not tested
elseif whichTest == 3 % squareBig -> roundBig
    impulseName0 = 'wav/squareBigCenterImp.wav';
    respName0 = 'wav/squareBigCenterResp.wav';
    impulseName1 = 'wav/roundBigCenterImp.wav';
    respName1 = 'wav/roundBigCenterResp.wav';
    impulseNameProof = 'wav/semBigCenterImp.wav';
    respNameProof = 'wav/semBigCenterResp.wav';
    interp = 0.42; % TODO: not tested
elseif whichTest == 4 % squareSmall -> roundSmall
    impulseName0 = 'wav/squareSmallCenterImp.wav';
    respName0 = 'wav/squareSmallCenterResp.wav';
    impulseName1 = 'wav/roundSmallCenterImp.wav';
    respName1 = 'wav/roundSmallCenterResp.wav';
    impulseNameProof = 'wav/semSmallCenterImp.wav';
    respNameProof = 'wav/semSmallCenterResp.wav';
    interp = 0.42; % TODO: not tested
elseif whichTest == 5 % squareBig -> roundSmall
    impulseName0 = 'wav/squareBigCenterImp.wav';
    respName0 = 'wav/squareBigCenterResp.wav';
    impulseName1 = 'wav/roundSmallCenterImp.wav';
    respName1 = 'wav/roundSmallCenterResp.wav';
    impulseNameProof = 'wav/semMidCenterImp.wav';
    respNameProof = 'wav/semMidCenterResp.wav';
    interp = 0.42; % TODO: not tested
elseif whichTest == 6 % roundBig -> squareSmall
    impulseName0 = 'wav/roundBigCenterImp.wav';
    respName0 = 'wav/roundBigCenterResp.wav';
    impulseName1 = 'wav/squareSmallCenterImp.wav';
    respName1 = 'wav/squareSmallCenterResp.wav';
    impulseNameProof = 'wav/semMidCenterImp.wav';
    respNameProof = 'wav/semMidCenterResp.wav';
    interp = 0.42; % TODO: not tested
end

%% compute STFTs and plot the response spectrogram

% load impulse and response
[impulse0, fs] = wavread(impulseName0);
[response0, fs] = wavread(respName0);
[impulse1, fs] = wavread(impulseName1);
[response1, fs] = wavread(respName1);
[impulseProof, fs] = wavread(impulseNameProof);
[responseProof, fs] = wavread(respNameProof);

impulseSTFT0 = stft(sum(impulse0,2),nbins,nskip);
% Uncomment to plot spectrogram
%figure(1);
%responseSTFT0 = ftgram(sum(response0,2), fs, 'music', 'nbins', nbins, 'nskip', nskip);
responseSTFT0 = stft(sum(response0,2),nbins, nskip); % Comment if plotting spectrogram

impulseSTFT1 = stft(sum(impulse1,2),nbins,nskip);
% Uncomment to plot spectrogram
%figure(2);
%responseSTFT1 = ftgram(sum(response1,2), fs, 'music', 'nbins', nbins, 'nskip', nskip);
responseSTFT1 = stft(sum(response1,2),nbins, nskip); % Comment if plotting spectrogram

impulseSTFTProof = stft(sum(impulseProof,2),nbins,nskip);
% Uncomment to plot spectrogram
%figure(2);
%responseSTFT1 = ftgram(sum(response1,2), fs, 'music', 'nbins', nbins, 'nskip', nskip);
responseSTFTProof = stft(sum(responseProof,2),nbins, nskip); % Comment if plotting spectrogram

%% compute impulse response and plot it
% TODO from here should do steps for 1 (with only have 0 here)

impulseSpectrum0 = mean(abs(impulseSTFT0),2)/max(mean(abs(impulseSTFT0),2));
responseSpectrum0 = mean(abs(responseSTFT0),2)/max(mean(abs(responseSTFT0),2));
irSpectrum0 = (responseSpectrum0./impulseSpectrum0); % deconvolution
irSpectrumDB0 = 20*log10(irSpectrum0/max(irSpectrum0)); % converting to dB

impulseSpectrum1 = mean(abs(impulseSTFT1),2)/max(mean(abs(impulseSTFT1),2));
responseSpectrum1 = mean(abs(responseSTFT1),2)/max(mean(abs(responseSTFT1),2));
irSpectrum1 = (responseSpectrum1./impulseSpectrum1); % deconvolution
irSpectrumDB1 = 20*log10(irSpectrum1/max(irSpectrum1)); % converting to dB

impulseSpectrumProof = mean(abs(impulseSTFTProof),2)/max(mean(abs(impulseSTFTProof),2));
responseSpectrumProof = mean(abs(responseSTFTProof),2)/max(mean(abs(responseSTFTProof),2));
irSpectrumProof = (responseSpectrumProof./impulseSpectrumProof); % deconvolution
irSpectrumDBProof = 20*log10(irSpectrumProof/max(irSpectrumProof)); % converting to dB

% define time, frequency axes
[~, nframes0] = size(responseSTFT0);
f0 = fs/2*[0:nbins]'/nbins;
t0 = [0.5:nframes0-0.5]*nskip/fs;

[~, nframes1] = size(responseSTFT1);
f1 = fs/2*[0:nbins]'/nbins;
t1 = [0.5:nframes1-0.5]*nskip/fs;

[~, nframesProof] = size(responseSTFTProof);
fProof = fs/2*[0:nbins]'/nbins;
tProof = [0.5:nframesProof-0.5]*nskip/fs;

%% find mode frequencies and amplitudes

[gammam0, im0] = localmax(irSpectrumDB0);
gammam0 = gammam0-max(gammam0);
fm0 = (im0-1)/nbins*fs/2;
nmode0 = length(fm0);

[gammam1, im1] = localmax(irSpectrumDB1);
gammam1 = gammam1-max(gammam1);
fm1 = (im1-1)/nbins*fs/2;
nmode1 = length(fm1);

[gammamProof, imProof] = localmax(irSpectrumDBProof);
gammamProof = gammamProof-max(gammamProof);
fmProof = (imProof-1)/nbins*fs/2;
nmodeProof = length(fmProof);

%% estimate T60s

delta0 = round(offset*fs/nskip);
rt60m0 = zeros(nmode0,1);
for m = [1:nmode0],
    % extract mode response
    psdm0 = mean(abs(responseSTFT0(max(1,round((nbins-1)*fm0(m)*2/fs)+[-beta:beta]),:)),1)';

    % estimate decay time, initial amplitude
    [gmax0, fmax0] = max(psdm0);
    fstart0 = fmax0 + delta0;
    fstop0 = max(find(psdm0 > rmin*gmax0, 1, 'last'), fstart0+round(tmin*fs/nskip));
    index0 = [fstart0+1:fstop0];

    basis0 = [ones(fstop0-fstart0,1) t0(index0)'];
    theta0 = basis0\(20*log10(psdm0(index0)));

    rt60m0(m) = -60/theta0(2);
end;

delta1 = round(offset*fs/nskip);
rt60m1 = zeros(nmode1,1);
for m = [1:nmode1],
    % extract mode response
    psdm1 = mean(abs(responseSTFT1(max(1,round((nbins-1)*fm1(m)*2/fs)+[-beta:beta]),:)),1)';

    % estimate decay time, initial amplitude
    [gmax1, fmax1] = max(psdm1);
    fstart1 = fmax1 + delta1;
    fstop1 = max(find(psdm1 > rmin*gmax1, 1, 'last'), fstart1+round(tmin*fs/nskip));
    index1 = [fstart1+1:fstop1];

    basis1 = [ones(fstop1-fstart1,1) t1(index1)'];
    theta1 = basis1\(20*log10(psdm1(index1)));

    rt60m1(m) = -60/theta1(2);
end;

deltaProof = round(offset*fs/nskip);
rt60mProof = zeros(nmodeProof,1);
for m = [1:nmodeProof],
    % extract mode response
    psdmProof = mean(abs(responseSTFTProof(max(1,round((nbins-1)*fmProof(m)*2/fs)+[-beta:beta]),:)),1)';

    % estimate decay time, initial amplitude
    [gmaxProof, fmaxProof] = max(psdmProof);
    fstartProof = fmaxProof + deltaProof;
    fstopProof = max(find(psdmProof > rmin*gmaxProof, 1, 'last'), fstartProof+round(tmin*fs/nskip));
    indexProof = [fstartProof+1:fstopProof];

    basisProof = [ones(fstopProof-fstartProof,1) tProof(indexProof)'];
    thetaProof = basisProof\(20*log10(psdmProof(indexProof)));

    rt60mProof(m) = -60/thetaProof(2);
end;

%% detecting the most powerful modes 

sorted_gammam0 = sort(gammam0,'descend');
maxIndex0 = find(ismember(gammam0,sorted_gammam0(1:nmodes)));

sorted_gammam1 = sort(gammam1,'descend');
maxIndex1 = find(ismember(gammam1,sorted_gammam1(1:nmodes)));

sorted_gammamProof = sort(gammamProof,'descend');
maxIndexProof = find(ismember(gammamProof,sorted_gammamProof(1:nmodesProof)));

%% Pairing Modes

% Measuring the evolution of the first mode
transpositionIndex = fm1(maxIndex1(1)) - fm0(maxIndex0(1));

scaledFm0 = fm0 + transpositionIndex;
pairedIndex1 = zeros(nmodes,1);
for m = [1:nmodes]
    [idx idx] = min(abs(fm1-scaledFm0(maxIndex0(m))));
    pairedIndex1(m) = idx;
end

scaledFm1 = fm1 - transpositionIndex;
pairedIndex0 = zeros(nmodes,1);
for m = [1:nmodes]
    [idx idx] = min(abs(fm0-scaledFm1(maxIndex1(m))));
    pairedIndex0(m) = idx;
end

clear matchedIndex0;
clear matchedIndex1;
cnt = 1;
for m = [1:nmodes]
    if modeMatching == 1
        if maxIndex0(m) == pairedIndex0(m) && maxIndex1(m) == pairedIndex1(m)
            matchedIndex0(cnt) = pairedIndex0(m);
            matchedIndex1(cnt) = pairedIndex1(m);
            cnt = cnt+1;
        end
    else 
        matchedIndex0(cnt) = pairedIndex0(m);
        matchedIndex1(cnt) = pairedIndex1(m);
        cnt = cnt+1;
    end
end

disp(sprintf('%i modes were matched',(cnt-1)));

%% Computing Theoretical Response

fmMorph = fm0(matchedIndex0) + (fm1(matchedIndex1)-fm0(matchedIndex0))*interp;
gammamMorph = gammam0(matchedIndex0) + (gammam1(matchedIndex1)-gammam0(matchedIndex0))*interp;
rt60mMorph = rt60m0(matchedIndex0) + (rt60m1(matchedIndex1)-rt60m0(matchedIndex0))*interp;

%% Plotting Frequeny Responses and Detected Modes 

% figure(1);
% plot(f0, irSpectrumDB0, '-', fm0(maxIndex0), gammam0(maxIndex0), 'o', fm0(matchedIndex0), gammam0(matchedIndex0), 'ro'); grid;
% title('Detected Modes for Object 0');
% xlabel('frequency, Hz'); ylabel('power, dB');
% xlim([20 5000]);
% 
% figure(2);
% plot(f1, irSpectrumDB1, '-', fm1(maxIndex1), gammam1(maxIndex1), 'o', fm1(matchedIndex1), gammam1(matchedIndex1), 'o'); grid;
% title('Detected Modes for Object 1');
% xlabel('frequency, Hz'); ylabel('power, dB');
% xlim([20 5000]);

%% Plotting measured modes T60s and amplitudes

% figure(3);
% subplot(2,1,1);
% stem(fm0(matchedIndex0), gammam0(matchedIndex0)); grid;
% title('Mode Amplitudes');
% xlabel('mode frequency, Hz'); ylabel('amplitude, dB');
% xlim([20 5000]);
% subplot(2,1,2); 
% stem(fm0(matchedIndex0), rt60m0(matchedIndex0)); grid;
% title('Mode T60s');
% xlabel('mode frequency, Hz'); ylabel('60 dB decay time, seconds');
% xlim([20 5000]);
% 
% figure(4);
% subplot(2,1,1);
% stem(fm1(matchedIndex1), gammam1(matchedIndex1)); grid;
% title('Mode Amplitudes');
% xlabel('mode frequency, Hz'); ylabel('amplitude, dB');
% xlim([20 5000]);
% subplot(2,1,2); 
% stem(fm1(matchedIndex1), rt60m1(matchedIndex1)); grid;
% title('Mode T60s');
% xlabel('mode frequency, Hz'); ylabel('60 dB decay time, seconds');
% xlim([20 5000]);

%% Plotting transposed modes against measured modes 

% figure(4);
% subplot(2,1,1);
% oneValue = repmat(1,1,length(matchedIndex1));
% stem(fm1(matchedIndex1), oneValue, 'blue'); % measured modes
% hold on;
% stem(scaledFm0(matchedIndex0), oneValue, 'red'); % transposed modes
% hold off;
% grid;
% set(gca,'YTickLabel',[])
% xlabel('mode frequency, Hz');
% title('Object 1 ');
% xlim([20 2000]);
% subplot(2,1,2);
% stem(fm0(matchedIndex0), oneValue, 'blue'); % measured modes
% hold on;
% stem(scaledFm1(matchedIndex1), oneValue, 'red'); % transposed modes
% hold off;
% grid;
% set(gca,'YTickLabel',[])
% xlabel('mode frequency, Hz');
% title('Object 0');
% xlim([20 2000]);

%% Plotting morphed modes against measured modes for middle object

figure(5);
subplot(2,1,1);
stem(fmProof(maxIndexProof), gammamProof(maxIndexProof), 'blue');
hold on;
stem(fmMorph, gammamMorph, 'red');
hold off;
grid;
title('Mode Amplitudes');
xlabel('mode frequency, Hz'); ylabel('amplitude, dB');
xlim([20 3000]);
subplot(2,1,2); 
stem(fmProof(maxIndexProof), rt60mProof(maxIndexProof), 'blue');
hold on;
stem(fmMorph, rt60mMorph, 'red');
hold off;
grid;
title('Mode T60s');
xlabel('mode frequency, Hz'); ylabel('60 dB decay time, seconds');
xlim([20 3000]);

