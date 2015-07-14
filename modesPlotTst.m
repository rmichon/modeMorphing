% fs = 44100;
% 
% fm = [400 1000];
% rt60m = [0.1 1];
% gm = [1 1];
% nmodes = length(fm);

for m = [1:nmode],

    w = 2*pi*fm(m)/fs;
    r = (1/rt60m(m)*fs).^(0.001);
    A = [1 -2*r*cos(w) r^2];
    B = [1 0 -1];

    [H,F] = freqz(B,A,1024);
    %HdB = dbn(H, -100);
    HdB = 20*log10(H*gm(m));
    fplot = [0:length(F)-1]/length(F)*fs/2;
    plot(fplot,real(HdB));
    hold on;
    %xlim([20 10000]);

end;

hold off;
xlim([20 10000]);