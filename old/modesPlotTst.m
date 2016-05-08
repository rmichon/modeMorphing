% fs = 44100;
% 
% fm = [400 1000];
% rt60m = [0.1 1];
% gm = [1 1];
% nmodes = length(fm);

figure(1);
for m = [1:nmode],

    w = 2*pi*fm(m)/fs;
    r = 0.001.^(1/(rt60m(m)*fs));
    A = [1 -2*r*cos(w) r^2];
    B = [gm(m) 0 -gm(m)];

    [H,F] = freqz(B,A,1024);
    %HdB = dbn(H, -100);
    HdB = 20*log10(H);
    fplot = [0:length(F)-1]/length(F)*fs/2;
    plot(fplot,real(HdB));
    hold on;
    if m == 1
        sos = [B A];
    else
        sos = [sos; B A];
    end

end;

hold off;
xlim([20 10000]);

%% Trying the analytical approach to plot the frequency response

[B,A] = sos2tf(sos);
[H,F] = freqz(B,A,2048);
HdB = 20*log10(H);
fplot = [0:length(F)-1]/length(F)*fs/2;

% figure(2);
% plot(fplot,real(HdB));
% xlim([20 10000]);

%% Impulse approach

fs = 96000;
dur = 0.4;

impulse = [1 zeros(1,fs*dur-1)];
response = zeros(1,fs*dur);

for m = [1:nmode],
    B = sos(m,1:3);
    A = sos(m,4:6);
    response = response + filter(B,A,impulse);
end;

response = response/max(response);
audiowrite('res.wav',response,fs);

figure(6);
plot(response);

nbins = 2048;
nskip = 32; 

responseSTFT = stft(sum(response',2),nbins,nskip);
responseSpectrum = mean(abs(responseSTFT),2)/max(mean(abs(responseSTFT),2));

figure(7);
plot(f, 20*log10(responseSpectrum), '-'); grid;
title('Spectrum');
xlabel('frequency, Hz'); ylabel('power, dB');
xlim([20 10000]);







