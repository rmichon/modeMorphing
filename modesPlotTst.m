% fs = 44100;
% 
% fm = [400 1000];
% rt60m = [0.1 1];
% gm = [1 1];
% nmodes = length(fm);

%sos = 

figure(1);
for m = [1:nmode],

    w = 2*pi*fm(m)/fs;
    r = 0.001.^(1/(rt60m(m)*fs));
    A = [1 -2*r*cos(w) r^2];
    B = [1 0 -1];

    [H,F] = freqz(B,A,1024);
    %HdB = dbn(H, -100);
    HdB = 20*log10(H*gm(m));
    fplot = [0:length(F)-1]/length(F)*fs/2;
    plot(fplot,real(HdB));
    hold on;
    sos = [sos; B A];

end;

hold off;
xlim([20 10000]);

%% Trying the analytical approach to plot the frequency response

[B,A] = sos2tf(sos);
[H,F] = freqz(B,A,2048);
HdB = 20*log10(H);
fplot = [0:length(F)-1]/length(F)*fs/2;

figure(2);
plot(fplot,real(HdB));
xlim([20 10000]);

%% Impulse approach

%fs = 96000;
%dur = 0.4;

%impulse = [1 zeros(1,96000*0.4)];

%B = sos







