clc;
clear;
close all;

%PARAMETERS
Fs = 44100;
fc = 10e3;
N  = 128;     % subcarriers
L  = 20;      % CP length
T  = 58e-3;

data = load('J:\PWC\Assignment1\Signals_task3/signal6.mat');
R = data.R(:);
t = data.t(:);

figure; plot(t,R); title('Received Signal');

%Downconversion and LPF
I = R .* cos(2*pi*fc*t);
Q = -R .* sin(2*pi*fc*t);

[B,A] = butter(8, 0.05);
[H,F] = freqz(B,A,1024,Fs);

%Filter amplitude response
figure;
plot(F, abs(H), 'LineWidth', 2);
grid on;
xlabel('Frequency (Hz)', 'FontSize',14);
ylabel('Amplitude','FontSize',14);
title('Amplitude Frequency Response of LPF','FontSize',14);

I = filter(B,A,I);
Q = filter(B,A,Q);

r_bb = I + 1j*Q;

%Downsampling
Ns = round(Fs*T/N); %Original sample rate
r_bb = r_bb(1:Ns:end);

figure; 
plot(abs(r_bb)); 
title('Baseband');

%Synchronization (CP correlation)
corr_window = N/2;
corr = zeros(1, length(r_bb)-2*corr_window);

for k = 1:length(corr)
    a = r_bb(k:k+corr_window-1);
    b = r_bb(k+corr_window:k+2*corr_window-1);

    corr(k) = abs(sum(a .* conj(b))) / sqrt(sum(abs(a).^2)*sum(abs(b).^2));
end

[~, start] = max(corr);

figure; 
plot(corr); 
fprintf('Start index = %d\n', start); %check start index

%OFDM symbol extraction
sym_len = N + L;
start = start - 5;   % small offset inside CP (trial and error)

num_sym = floor((length(r_bb)-start)/sym_len);

X = zeros(N, num_sym);

for i = 1:num_sym
    idx = start + (i-1)*sym_len;
    block = r_bb(idx:idx+N-1);   % take FFT window
    X(:,i) = fft(block);
end

fprintf('Number of OFDM symbols = %d\n', num_sym);

%Channel estimation
randn('state',100);   % keep SAME as transmitter
P = sign(randn(1,N/2));

Xpilot = zeros(N,1);
Xpilot(1:2:end) = 2*P;

Ypilot = X(:,1); 

idx = 1:2:N;
H_est = zeros(N,1);
H_est(idx) = Ypilot(idx) ./ Xpilot(idx);

H_est = interp1(idx, H_est(idx), 1:N, 'linear', 'extrap').';

%Equalization
Y = X(:,2:end);
Xeq = Y ./ H_est;

figure;
plot(real(Xeq(:)), imag(Xeq(:)), '.');
title('Constellation'); 
grid on;

%Length decoding
len_sym = Xeq(:,1);

bits_len = zeros(1, 2*N);

idx = 1;
for k = 1:N
    I = real(len_sym(k));
    Q = imag(len_sym(k));
    
    if (I >= 0 && Q >= 0)
        bits_len(idx:idx+1) = [0 0];
    elseif (I < 0 && Q > 0)
        bits_len(idx:idx+1) = [1 0];
    elseif (I < 0 && Q < 0)
        bits_len(idx:idx+1) = [1 1];
    elseif (I >= 0 && Q < 0)
        bits_len(idx:idx+1) = [0 1];
    end
    
    idx = idx + 2;
end

trellis = poly2trellis(6,[77 45]);
len_dec = vitdec(bits_len(1:30), trellis, 15, 'term','hard');

len_bits = len_dec(1:10);

%Trying both MSB and LSB bit order
lm1 = len_bits * (2.^(9:-1:0)).';
lm2 = len_bits * (2.^(0:9)).';

fprintf('Length1 = %d, Length2 = %d\n', lm1, lm2);

lm = lm2;   %try both to see which decodes

fprintf('Message length = %d chars\n', lm);

%Demodulation
payload = Xeq(:,2:end);
payload = payload(:);

Nsymbols = length(payload);
bits = zeros(1, 2*Nsymbols);

idx = 1;
for k = 1:Nsymbols
    I = real(payload(k));
    Q = imag(payload(k));
    
    if (I >= 0 && Q >= 0)
        bits(idx:idx+1) = [0 0];
    elseif (I < 0 && Q >= 0)
        bits(idx:idx+1) = [1 0];
    elseif (I < 0 && Q < 0)
        bits(idx:idx+1) = [1 1];
    elseif (I >= 0 && Q < 0)
        bits(idx:idx+1) = [0 1];
    end
    
    idx = idx + 2;
end

%Viterbi decoding
num_bits = 7*lm;
mem = 5;
total_in = num_bits + mem;
needed = 2*total_in;

bits = bits(1:needed);

decoded = vitdec(double(bits), trellis, 30, 'term','hard');
decoded = decoded(1:num_bits);

%ASCII Decoding
wh = 2.^(6:-1:0);
num_chars = floor(length(decoded)/7);
m = '';
for l = 1:num_chars
    bits7 = decoded(7*(l-1)+1 : 7*l);
    m = [m char(bits7 * wh')];
end

disp(m);