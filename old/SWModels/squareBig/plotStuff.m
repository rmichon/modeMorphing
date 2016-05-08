format shortG
freqFile = fopen('freq.txt','r');
freq = fscanf(freqFile,'%f');
%freqTr = freq*0.92958;
freqTr = freq - 11;

massXFile = fopen('massX.txt','r');
massX = fscanf(massXFile,'%f');

massYFile = fopen('massY.txt','r');
massY = fscanf(massYFile,'%f');

massZFile = fopen('massZ.txt','r');
massZ = fscanf(massZFile,'%f');

mass = massX + massY + massZ;


figure(10);
plot(freq,mass);
%plot(freq,massX,'g',freq,massY,'r',freq,massZ,'b');
xlim([0 5000]);
legend('X','Y','Z');