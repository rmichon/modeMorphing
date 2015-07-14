% TODO: should turn that into a function that plots an interactive graph

interp = 0.5;

file1Fm = fopen('modes/squareBigCenterFreq.txt','r');
fm1 = fscanf(file1Fm,'%f');
file1Gm = fopen('modes/squareBigCenterGain.txt','r');
gm1 = fscanf(file1Gm,'%f');
file1Rt60m = fopen('modes/squareBigCenterT60.txt','r');
rt60m1 = fscanf(file1Rt60m,'%f');

file2Fm = fopen('modes/squareSmallCenterFreq.txt','r');
fm2 = fscanf(file2Fm,'%f');
file2Gm = fopen('modes/squareSmallCenterGain.txt','r');
gm2 = fscanf(file2Gm,'%f');
file2Rt60m = fopen('modes/squareSmallCenterT60.txt','r');
rt60m2 = fscanf(file2Rt60m,'%f');

fm1 = fm1(1:length(fm2));
gm1 = gm1(1:length(gm2));
rt60m1 = rt60m1(1:length(rt60m2));

fm = fm1 + (fm2-fm1)*interp;
gm = gm1 + (gm2-gm1)*interp;
rt60m = rt60m1 + (rt60m2-rt60m1)*interp;

gm1 = gm1/max(gm1);
gmdB1 = 20*log10(gm1);
gmZerodB1 = gmdB1-min(gmdB1);

gm2 = gm2/max(gm2);
gmdB2 = 20*log10(gm2);
gmZerodB2 = gmdB2-min(gmdB2);

gm = gm/max(gm);
gmdB = 20*log10(gm);
gmZerodB = gmdB-min(gmdB);

figure(1);
subplot(2,1,1);
stem(fm, gmZerodB);
%hold on;
%stem(fm1, gmZerodB1, 'blue');
%stem(fm2, gmZerodB2, 'red');
grid;
%hold off;
title('Mode Amplitudes');
xlabel('mode frequency, Hz'); ylabel('amplitude, dB');

subplot(2,1,2); 
stem(fm, rt60m); grid;
title('Mode T60s');
xlabel('mode frequency, Hz'); ylabel('60 dB decay time, seconds');
