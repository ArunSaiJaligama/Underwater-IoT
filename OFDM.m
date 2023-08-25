snr_dB = -10:2:20;
for q = 1:1:length(snr_dB)
    snr(q) = 10^(snr_dB(q)/10); 
end

 % rx_grid = zeros(M,N);
 BER = zeros(1,length(snr));

for j = 1:1:length(snr)
    M  = 64;
    l = 5;
    k = log2(16);
    n = 6400;
    cp = l;
    x = randi([0,1],n*k,1);
    in_symbols = bit2int(x,k);
    x_mod = qammod(in_symbols,16,'bin');
    N = length(x_mod)/M;
    rx_grid = zeros(M,N);
    h = sqrt(1/2)*(randn(1,l)+1i*randn(1,l));
    H = fft(h,M);
    h = [h zeros(1,M-l)];
    grid = reshape(x_mod,[M,N]);
    noise = sqrt(1/(2*snr(j)))*(randn(1,M)+1i*randn(1,M));
    for n = 1:1:N
        tx_sym = grid(:,n);
        tx_sym_ifft = ifft(tx_sym,M);
        tx_sym_ifft = tx_sym_ifft.';
        tx_cp = [tx_sym_ifft(M-cp+1:M) tx_sym_ifft];
        rx = cconv(tx_cp,h,M)+noise;
        %rx_with_noise = awgn(rx,snr(j),"measured");
        % rx_no_cp = rx_with_noise(cp+1:cp+M);
        rx_final = fft(rx,M);%change
        rx_equ = rx_final./H;
        rx_grid(:,n) = rx_equ;
    end
    y_demod = reshape(rx_grid,[M*N,1]);
    outsymbols = qamdemod(y_demod,16,'bin');
    y = int2bit(outsymbols,k);
    [numErrors,ber] = biterr(x,y);
    BER(j) = ber;
end

semilogy(snr_dB,BER,'r');
xlabel('snr [dB]');
ylabel('Bit Error Rate');
title('OFDM');