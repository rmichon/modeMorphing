% TODO: should turn that into a function that plots an interactive graph

file1Fm = fopen('modes/squareBigCenterFreq.txt','r');
fm1 = fscanf(file1Fm,'%f');
file1Gm = fopen('modes/squareBigCenterGain.txt','r');
gm1 = fscanf(file1Gm,'%f');
file1Rt60m = fopen('modes/squareBigCenterT60.txt','r');
rt60m1 = fscanf(file1Rt60m,'%f');

file2Fm = fopen('modes/squareBigCenterFreq.txt','r');
fm2 = fscanf(file2Fm,'%f');
file2Gm = fopen('modes/squareBigCenterGain.txt','r');
gm2 = fscanf(file2Gm,'%f');
file2Rt60m = fopen('modes/squareBigCenterT60.txt','r');
rt60m2 = fscanf(file2Rt60m,'%f');
%hi

gm1 = gm1/max(gm1);
gmdB1 = 20*log10(gm1);
gmZerodB1 = gmdB1-min(gmdB1);

gm2 = gm2/max(gm2);
gmdB2 = 20*log10(gm2);
gmZerodB2 = gmdB2-min(gmdB2);

figure(1);
subplot(2,1,1);
stem(fm1, gmZerodB1, 'blue');
hold on;
stem(fm2, gmZerodB2, 'red');
hold off
grid;
title('Mode Amplitudes');
xlabel('mode frequency, Hz'); ylabel('amplitude, dB');

subplot(2,1,2); 
stem(fm1, rt60m1); grid;
title('Mode T60s');
xlabel('mode frequency, Hz'); ylabel('60 dB decay time, seconds');
