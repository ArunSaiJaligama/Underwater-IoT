% number of Doppler bins (time slots)
N=16;
% number of delay bins (subcarriers)
M=64;
% normalized DFT matrix
Fn=dftmtx(N);
Fn=Fn/norm(Fn);
% subcarrier spacing
delta_f=15e3;
% block duration
T=1/delta_f;
% carrier frequency
fc=4e9;
% speed of light
c=299792458;
% OTFS grid delay and Doppler resolution
delay_resolution = 1/(M*delta_f);
Doppler_resolution = 1/(N*T);
%modulation size
mod_size=4;
% number of information symbols in one frame
N_syms_per_frame=N*M;
% number of information bits in one frame
N_bits_per_frame=N*M*log2(mod_size);
% generate random bits
tx_info_bits=randi([0,1],N_bits_per_frame,1);
% QAM modulation
tx_info_symbols=qammod(tx_info_bits,mod_size,'gray','InputType','bit');
% Generate the MxN OTFS delay-Doppler frame
X=reshape(tx_info_symbols,M,N);
X_tilda=X*Fn';
s=reshape(X_tilda,1,N*M);

Es = mean(abs(qammod(0:mod_size-1,mod_size).^2));
% SNR=Es/noise power
SNR_dB = 25;
SNR=10.^(SNR_dB/10);
% noise power
sigma_w_2=Es/SNR;
% generate Gaussian noise samples with variance=sigma_w_2
noise = sqrt(sigma_w_2/2)*(randn(N*M,1) + 1i*randn(N*M,1));
% add AWGN to the received signal
r=s+noise;

Y_tilda=reshape(r,[M,N]);
Y=Y_tilda*Fn;
% vectorize Y
y=reshape(Y.',N*M,1);
% QAM demodulation
x_hat=qamdemod(x_hat,mod_size,'gray');
% [numErrors,ber] = biterr(x_hat,y);

