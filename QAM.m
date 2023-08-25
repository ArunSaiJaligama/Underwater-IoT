M  = 16;
k = log2(M);
n = 10000;
sps = 1;
snr_dB = 12;
snr = 10^(snr_dB/10);
x = randi([0,1],n*k,1);
symbols = bit2int(x,k);
x_mod = qammod(symbols,M,'bin');
receivedSignal = awgn(x_mod,snr,'measured');
scatfig = scatterplot(receivedSignal,1,0,'b.');
hold on
scatterplot(x_mod,1,0,'r*',scatfig);
hold off