% TODO: should turn that into a function that plots an interactive graph

%% Parameters

fs = 96000;
%interp = 0.42; % round
interp = 0.47; % semi
%interp = 0.55; % square
nbins = 2048; % FFT window size for theoritical IR
nskip = 32; % FFT hop size for theoritical IR

% change the file names in function of your needs...
% file1Fm = fopen('modes/squareBigCenterFreq.txt','r');
% file1Gm = fopen('modes/squareBigCenterGain.txt','r');
% file1Rt60m = fopen('modes/squareBigCenterT60.txt','r');
% file2Fm = fopen('modes/squareSmallCenterFreq.txt','r');
% file2Gm = fopen('modes/squareSmallCenterGain.txt','r');
% file2Rt60m = fopen('modes/squareSmallCenterT60.txt','r');

file1Fm = fopen('modes/semBigCenterFreq.txt','r');
file1Gm = fopen('modes/semBigCenterGain.txt','r');
file1Rt60m = fopen('modes/semBigCenterT60.txt','r');
file2Fm = fopen('modes/semSmallCenterFreq.txt','r');
file2Gm = fopen('modes/semSmallCenterGain.txt','r');
file2Rt60m = fopen('modes/semSmallCenterT60.txt','r');

% file1Fm = fopen('modes/roundBigCenterFreq.txt','r');
% file1Gm = fopen('modes/roundBigCenterGain.txt','r');
% file1Rt60m = fopen('modes/roundBigCenterT60.txt','r');
% file2Fm = fopen('modes/roundSmallCenterFreq.txt','r');
% file2Gm = fopen('modes/roundSmallCenterGain.txt','r');
% file2Rt60m = fopen('modes/roundSmallCenterT60.txt','r');


%% Extracting the data & computing theoritical response

fm1 = fscanf(file1Fm,'%f'); % mode frequencies
gm1 = fscanf(file1Gm,'%f'); % mode gains
rt60m1 = fscanf(file1Rt60m,'%f'); % mode T60

fm2 = fscanf(file2Fm,'%f');
gm2 = fscanf(file2Gm,'%f');
rt60m2 = fscanf(file2Rt60m,'%f');

% Adjusting the size of 1 in function of 2 (in case...)
nmode = length(fm2);
fm1 = fm1(1:nmode);
gm1 = gm1(1:nmode);
rt60m1 = rt60m1(1:nmode);

% calculating the parameters of the theoritical response
modeFreqs = fm1 + (fm2-fm1)*interp;
modeGains = gm1 + (gm2-gm1)*interp;
modeT60 = rt60m1 + (rt60m2-rt60m1)*interp;

%% Ploting the theoritical modes

% converting the modes gains to normalized amplitudes
% gm1 = gm1/max(gm1);
% gmdB1 = 20*log10(gm1);
% gmZerodB1 = gmdB1-min(gmdB1);
% 
% gm2 = gm2/max(gm2);
% gmdB2 = 20*log10(gm2);
% gmZerodB2 = gmdB2-min(gmdB2);

modeGains = modeGains/max(modeGains);
gmdB = 20*log10(modeGains);
gmZerodB = gmdB-min(gmdB);

figure(20);
stem(modeFreqs, gmdB, 'blue');
hold on;
stem(fm, 20*log10(gm), 'red');
hold off;
grid;
title('Mode Amplitudes');
xlabel('mode frequency, Hz'); ylabel('amplitude, dB');

% figure(21);
% stem(modeFreqs, modeT60, 'blue');
% hold on;
% stem(fm, rt60m, 'red'); 
% hold off;
% grid;
% title('Mode T60s');
% xlabel('mode frequency, Hz'); ylabel('60 dB decay time, seconds');

%% Converting the modes into a matrix of biquad coeffs

nPoints = 4096; % Number of points for freqz
totalH = zeros(nPoints,1); % Frequency response of the modeled instrument

for m = [1:nmode],
    
    % freq, gain and T60s are converted to biquad coeffs 
    w = 2*pi*modeFreqs(m)/fs;
    r = 0.001.^(1/(modeT60(m)*fs));
    A = [1 -2*r*cos(w) r^2];
    B = [modeGains(m) 0 -modeGains(m)];
    %B = [gm(m) 0 -1];
    
    % computing the frequency response of the biquad and plotting it 
    [H,F] = freqz(B,A,nPoints);
    totalH = totalH + H; % adding the frequency responses together
    
    % creating the coeffs matrix
    if m == 1
        sos = [B A];
    else
        sos = [sos; B A];
    end

end;

fplot = fs/2*[0:(nPoints-1)]'/(nPoints-1);

totalH = totalH / max(totalH); % normalizing the frequency response
HFindB = 20*log10(totalH); % to dB

irSpectrum = irSpectrum ./ max(irSpectrum);

% figure(2);
% plot(fplot,real(HFindB),'-',f, 20*log10(abs(irSpectrum)),'-'); grid;
% title('Computed Frequency Response of the Model');
% xlabel('Frequency, Hz'); ylabel('Amplitude, dB');
% xlim([20 10000]);

%% Computing the impulse response of the system and saving it to audio

dur = 0.4; % the duration of the response in seconds

impulse = [1 zeros(1,fs*dur-1)]; % generating the impulse
response = zeros(1,fs*dur); % instantiating the response

% computing the impulse response
for m = [1:nmode],
    B = sos(m,1:3);
    A = sos(m,4:6);
    response = response + filter(B,A,impulse);
end;

% normalizing it and writing it to an audio file
response = response/max(response);
audiowrite('res.wav',response,fs);

% plotting the response
% figure(3);
% plot(response);
% title('Theoritical Impulse Response');
% xlabel('Time, Samples'); ylabel('Gain');

%% Computing the frequency response of the impulse response

responseSTFT = stft(sum(response',2),nbins,nskip);
responseSpectrum = mean(abs(responseSTFT),2)/max(mean(abs(responseSTFT),2));
f = fs/2*[0:nbins]'/nbins;

% figure(4);
% plot(f, 20*log10(responseSpectrum), '-',f, 20*log10(abs(irSpectrum)),'-'); grid;
% title('Measured Frequency Response of the Model');
% xlabel('frequency, Hz'); ylabel('power, dB');
% xlim([20 10000]);
% ylim([-100 0]);









