%tx1_audio_ofdm.m
%clear; clc; close all;

params = ofdm_params();

Fs = params.Fs;
Tsamp = params.Tsamp;
fc = params.fc;
N_sc = params.N_sc;
N_cp = params.N_cp;
T_ofdm = params.T_ofdm;
F_ofdm = params.F_ofdm;
T_cp = params.T_cp;
T_sym = params.T_sym;
L = params.L;
delta_f = params.delta_f;

decoding_active = params.decoding_active;

header_length = params.header_length;

fid = fopen('input.txt','r');
fileData = fread(fid,'*uint8');
fclose(fid);
file_len = length(fileData);
fprintf('File: %d bytes = %d bits\n', length(fileData), length(fileData)*8);

if file_len*8 < 20000
    warning('File too small! Need >= 20000 bits. Have %d bits', file_len*8);
end

%convert bytes to bits
data_bits = double(reshape(de2bi(fileData, 8, 'left-msb').', 1, []));

%32-bit header with file length for receiver to know how many bytes to recover
header = double(de2bi(file_len, 32, 'left-msb'));
header = header(:).';

%header + data from file
bits = [header  data_bits];
Nbits = length(bits);

if mod(Nbits, 2) ~= 0
    bits  = [bits 0];   % pad one bit to make even
    Nbits = length(bits);
end

fprintf('File: %d bytes = %d bits\n', file_len, Nbits);

%Channel Coding
[coded_bits, N_coded_bits] = channel_coding.encode(bits, decoding_active);

%QPSK Modulation
[a, N_mod_symbols] = qpsk_modulation.modulate(coded_bits, params);

%ofdm symbol structure
N_active = params.N_active; % data include continous pilot
N_dc = params.N_dc;
N_outer_zero =  N_sc-N_active-N_dc; % 24
N_zero_each_side = N_outer_zero/2; % 12
pilot_idx = params.pilot_idx;
pilot_value = params.pilot_value;
N_data = N_active-length(pilot_idx); % real data without pilot

left_guard_idx = params.left_guard_idx;
dc_idx = params.dc_idx;
right_guard_idx = params.right_guard_idx;

reserved_idx = [left_guard_idx dc_idx right_guard_idx pilot_idx];
reserved_idx = sort(reserved_idx);

data_idx = setdiff(1:N_sc, reserved_idx);
data_idx = sort(data_idx);
Nofdm_data_symbols = ceil(N_mod_symbols/N_data); %OFDM data sysmbols number
fprintf('Num of OFDM Symbol: %d\n', Nofdm_data_symbols);

%last with zero
Npad = Nofdm_data_symbols*N_data - length(a);
a_padded = [a zeros(1, Npad)];

%mapping to one OFDM symbols
ofdm_symbols = zeros(Nofdm_data_symbols, N_sc);
for i = 1:Nofdm_data_symbols
    one_ofdm_symbol = zeros(1,N_sc);

    start_idx = (i-1)*N_data + 1;
    end_idx   = i*N_data;

    data_for_this_ofdm = a_padded(start_idx:end_idx);

    one_ofdm_symbol(data_idx) = data_for_this_ofdm;
    one_ofdm_symbol(pilot_idx) = pilot_value;

    ofdm_symbols(i,:) = one_ofdm_symbol;
end

%Data IFFT and CP
x = ifft(ofdm_symbols, N_sc, 2);
x_cp = [x(:, end-N_cp+1:end) x];

%Paralel to serial
tx_data_bb = reshape(x_cp.', 1, []);

%Header Generation
header_bits = de2bi(Nbits, header_length, 'left-msb'); %left most big value

% Channel Coding
[coded_header_bits, N_coded_header_bits] = channel_coding.encode(header_bits, decoding_active);

% QPSK Modulation
[header_symbols, N_header_mod_symbols] = qpsk_modulation.modulate(coded_header_bits, params);

% padded the last remain bit with 0 one last symbol
N_header_pad = N_data - N_header_mod_symbols;
header_symbols_padded = [header_symbols zeros(1, N_header_pad)];

% make one OFDM symbol
header_ofdm_symbol = zeros(1, N_sc);
header_ofdm_symbol(data_idx) = header_symbols_padded;
header_ofdm_symbol(pilot_idx) = pilot_value;

% Header IFFT and CP
x_header = ifft(header_ofdm_symbol, N_sc);
x_header_cp = [x_header(end-N_cp+1:end) x_header];

%% Preamble Generation
X_preamble = zeros(1, N_sc);
rng(100, 'twister');
P = sign(randn(1, N_sc/2));
X_preamble(1:2:end) = 2 * P;

% Preamble IFFT and CP
x_preamble = ifft(X_preamble, N_sc);
x_preamble_cp = [x_preamble(end-N_cp+1:end) x_preamble];

% combine preamble, header, data, and tail guard
N_tail_guard = params.N_tail_guard;
ofdm_symbol_len = N_sc + N_cp;
tx_bb = [x_preamble_cp, x_header_cp, tx_data_bb, zeros(1, N_tail_guard*ofdm_symbol_len)];

% Upsampling
tx_bb_up = interp(tx_bb, L);

n = 0:length(tx_bb_up)-1;
t = n/Fs;

%Upconversion
I = real(tx_bb_up);
Q = imag(tx_bb_up);

tx_audio = I.*cos(2*pi*fc*t) - Q.*sin(2*pi*fc*t);

%Normalization for audio channel (0.8 or 0.9) 
tx_audio = tx_audio / max(abs(tx_audio)) * 0.8;

%Bit Calculation
N_total_ofdm_symbols = 1 + 1 + Nofdm_data_symbols + N_tail_guard; % preamble + header + data + tail
T_total_theory = N_total_ofdm_symbols * T_sym;
T_total_audio = length(tx_audio)/Fs;

bit_rate = Nbits/T_total_audio;

fprintf('Theoretical frame duration: %.3f s\n', T_total_theory);
fprintf('Audio frame duration: %.3f s\n', T_total_audio);
fprintf('bit rate: %.2f kbps\n', bit_rate/1000);

%PLOT 1 - TX Signal spectrum
figure;
t_plot = (0:length(tx_audio)-1)/Fs;
plot(t_plot, tx_audio, 'LineWidth', 0.5);
title('Transmitted Signal (Time Domain)', 'FontSize', 14);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Amplitude', 'FontSize', 12);
grid on;

figure;
[Pxx_tx, F_tx] = pwelch(tx_audio, hamming(1024), 512, 1024, Fs);
plot(F_tx/1000, 10*log10(Pxx_tx), 'b', 'LineWidth', 1.5);
title('Transmitted Signal Spectrum', 'FontSize', 14);
xlabel('Frequency (kHz)', 'FontSize', 12);
ylabel('PSD (dB/Hz)', 'FontSize', 12);
xline(fc/1000, 'r--', sprintf('fc = %d Hz', fc), 'FontSize', 11);
grid on;
xlim([0 Fs/2000]);

%PAPR
signal_power = abs(tx_audio).^2;
PAPR_dB      = 10*log10(max(signal_power)/mean(signal_power));
fprintf('PAPR = %.2f dB\n', PAPR_dB);
if PAPR_dB > 10
    fprintf('WARNING: High PAPR (%.1f dB) — consider scrambler\n', PAPR_dB);
end

%Save original data
save('tx_config.mat','tx_audio','bits', ...
     'N_coded_bits','N_mod_symbols','bit_rate','coded_bits');

audiowrite('ofdm_tx.wav', tx_audio, Fs);
fprintf('Saved ofdm_tx.wav\n');