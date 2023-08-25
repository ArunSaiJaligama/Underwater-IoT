clear all
clc
snr_dB = -20:1:20;
snr = zeros(1,length(snr_dB));
for q = 1:1:length(snr_dB)
    snr(q) = 10^(snr_dB(q)/10); 
end

pe_analytical = zeros(1,length(snr));
for ii = 1:1:length(snr)
    MM = 128;
    ll = 16;
    pe_analytical(ii) = (1/2)*(1-sqrt((ll/MM)*snr(ii)/(1+(ll/MM)*snr(ii))));
    %pe_analytical(ii) = qfunc(sqrt((ll/MM)*snr(ii)));
end
tic;
 % rx_grid = zeros(M,N);
 BER = zeros(1,length(snr));
    M  = 128;
    l = 16;
    index = 2;
    %k = log2(index);
    %n = 128;
    cp = l;
for j = 1:1:length(snr)
    avg_ber = zeros(1,1000);
    for loop= 1:1:1000
        x = randi([0,1],M,1); %data
        %in_symbols = bit2int(x,k);
        %x_mod = qammod(in_symbols,index,'bin');
        x_mod = zeros(1,length(x));
        for jj=1:1:length(x)
            if(x(jj)==1)
                x_mod(jj) = 1;
            else
                x_mod(jj) = -1;
            end
        end
        N = length(x_mod); % time instants
        %rx_grid = zeros(1,M); 
        h = sqrt(1/2)*(randn(1,l)+1i*randn(1,l));
        H = fft(h,M);
        %h = [h zeros(1,M-cp)];
        %grid = reshape(x_mod,[M,N]);
        noise = sqrt(1/(2*snr(j)))*(randn(1,M)+1i*randn(1,M));
            %tx_sym = grid(:,nn);
        tx_sym_ifft = ifft(x_mod,M);
            %tx_sym_ifft = tx_sym_ifft.';
        tx_cp = [tx_sym_ifft(M-cp+1:M) tx_sym_ifft];
        rx = conv(h,tx_cp);
            %rx_with_noise = awgn(rx,snr(j),"measured");
        rx_no_cp = rx(cp+1:cp+M)+noise;
        rx_final = fft(rx_no_cp,M);%change
        rx_equ = rx_final./H;
            %rx_grid(:,nn) = rx_equ;
        %end
    %y_demod = reshape(rx_grid,[M*N,1]);
    %outsymbols = qamdemod(y_demod,index,'bin');
    %y = int2bit(outsymbols,k);
        y = zeros(length(rx_equ),1);
    for kk=1:1:length(rx_equ)
        if(real(rx_equ(kk))>0)
            y(kk) = 1;
        else
            y(kk) = 0;
        end
    end
    [numErrors,ber] = biterr(x,y);
    avg_ber(loop) = ber; 
    end
    BER(j) = double(mean(avg_ber));
end
toc;
hold on;
semilogy(snr_dB,BER,'r');
semilogy(snr_dB,pe_analytical,'b');
xlabel('snr [dB]');
ylabel('Bit Error Rate');
%xlim([-5 20]);
% ylim([10^(-5)]);
legend('BER Simulation','BER Analytical');
title('OFDM');