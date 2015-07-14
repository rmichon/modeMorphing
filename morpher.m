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

