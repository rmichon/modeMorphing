% TODO: should turn that into a function that plots an interactive graph

%% Parameters

fs = 96000;
interp = 0; % a value between 0 and 1
nbins = 2048; % FFT window size for theoritical IR
nskip = 32; % FFT hop size for theoritical IR

% change the file names in function of your needs...
file1Fm = fopen('modes/squareBigCenterFreq.txt','r');
file1Gm = fopen('modes/squareBigCenterGain.txt','r');
file1Rt60m = fopen('modes/squareBigCenterT60.txt','r');
file2Fm = fopen('modes/squareSmallCenterFreq.txt','r');
file2Gm = fopen('modes/squareSmallCenterGain.txt','r');
file2Rt60m = fopen('modes/squareSmallCenterT60.txt','r');

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
fm = fm1 + (fm2-fm1)*interp;
gm = gm1 + (gm2-gm1)*interp;
rt60m = rt60m1 + (rt60m2-rt60m1)*interp;

%% Ploting the theoritical modes

% converting the modes gains to normalized amplitudes
% gm1 = gm1/max(gm1);
% gmdB1 = 20*log10(gm1);
% gmZerodB1 = gmdB1-min(gmdB1);
% 
% gm2 = gm2/max(gm2);
% gmdB2 = 20*log10(gm2);
% gmZerodB2 = gmdB2-min(gmdB2);

gm = gm/max(gm);
gmdB = 20*log10(gm);
gmZerodB = gmdB-min(gmdB);

figure(1);
subplot(2,1,1);
stem(fm, gmZerodB); grid;
%hold on;
%stem(fm1, gmZerodB1, 'blue');
%stem(fm2, gmZerodB2, 'red');
%hold off;
title('Mode Amplitudes');
xlabel('mode frequency, Hz'); ylabel('amplitude, dB');

subplot(2,1,2); 
stem(fm, rt60m); grid;
title('Mode T60s');
xlabel('mode frequency, Hz'); ylabel('60 dB decay time, seconds');

%% Converting the modes into a matrix of biquad coeffs

figure(2);
for m = [1:nmode],
    
    % freq, gain and T60s are converted to biquad coeffs 
    w = 2*pi*fm(m)/fs;
    r = 0.001.^(1/(rt60m(m)*fs));
    A = [1 -2*r*cos(w) r^2];
    B = [gm(m) 0 -gm(m)];
    
    % computing the frequency response of the biquad and plotting it 
    [H,F] = freqz(B,A,1024);
    HdB = 20*log10(H);
    fplot = [0:length(F)-1]/length(F)*fs/2;
    plot(fplot,real(HdB));
    hold on;
    
    % creating the coeffs matrix
    if m == 1
        sos = [B A];
    else
        sos = [sos; B A];
    end

end;

hold off;
grid;
title('Biquads Frequency Response');
xlabel('Frequency, Hz'); ylabel('Amplitude, dB');
xlim([20 10000]);

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
figure(3);
plot(response);
title('Theoritical Impulse Response');
xlabel('Time, Samples'); ylabel('Gain');

%% Computing the frequency response of the impulse response

responseSTFT = stft(sum(response',2),nbins,nskip);
responseSpectrum = mean(abs(responseSTFT),2)/max(mean(abs(responseSTFT),2));
f = fs/2*[0:nbins]'/nbins;

figure(4);
plot(f, 20*log10(responseSpectrum), '-'); grid;
title('Spectrum');
xlabel('frequency, Hz'); ylabel('power, dB');
xlim([20 10000]);









