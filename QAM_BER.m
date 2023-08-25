snr_dB = -10:1:15;

for i = 1:1:length(snr_dB)
    snr(i) = 10^(snr_dB(i)/10); 
end

for j = 1:1:length(snr)
    EbNo = snr(j);
    M  = 64;
    k = log2(M);
    n = 50000;
    x = randi([0,1],n*k,1);
    in_symbols = bit2int(x,k);
    x_mod = qammod(in_symbols,M,'bin');
    receivedSignal = awgn(x_mod,EbNo,"measured");
    outsymbols = qamdemod(receivedSignal,M,'bin');
    y = int2bit(outsymbols,k);
    [numErrors,ber] = biterr(x,y);
    BER(j) = ber;
end

semilogy(snr_dB,BER,'r');
xlabel('snr [dB]');
ylabel('Bit Error Rate');
title('QAM modulation');