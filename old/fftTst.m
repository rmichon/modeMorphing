respName = 'wav/semBigCenterResp.wav';

[response, fs] = wavread(respName);

nb = 4096*4;
responseFFT = fft(response,nb);
responseFFT = responseFFT(1:round(length(responseFFT)/2));
irFFT = responseFFT;
irFFTdB = dbn(irFFT);
f = fs/2*[0:(length(responseFFT)-1)]'/(length(responseFFT)-1);
figure(10);
plot(f,irFFTdB); grid;
xlim([20 1000]);