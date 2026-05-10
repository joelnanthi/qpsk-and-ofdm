%rx1_audio_ofdm.m
clear; clc; close all;

%check if tx already run
if ~isfile('tx_config.mat')
    error('tx_config.mat not found');
end

if ~isfile('recorded_ofdm.wav')
    error('recorded_ofdm.wav not found');
end

%% Parameter
params = ofdm_params();

Fsamp = params.Fsamp;
fc = params.fc;
Nsc = params.Nsc;
Ncp = params.Ncp;
L = params.L;

pilot_idx = params.pilot_idx;
pilot_value = params.pilot_value;
data_idx = params.data_idx;

M = params.M;
bits_per_symbol = params.bits_per_symbol;
N_data = params.N_data;

decoding_active = params.decoding_active;

header_length = params.header_length;

%% Receiver side
load('tx_config.mat', ...
     'bits', ...
     'N_mod_symbols', ...
     'bit_rate');

[rx_audio, Fs_read] = audioread('recorded_ofdm.wav');
rx_audio = rx_audio.';

if Fs_read ~= Fsamp
    error('Sampling frequency mismatch');
end

n_rx = 0:length(rx_audio)-1;
t_rx = n_rx/Fsamp;

%PLOT 1 - RXD signal time domain
figure;
plot(t_rx, rx_audio, 'LineWidth', 0.3);
title('Received Signal (Time Domain)', 'FontSize', 14);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Amplitude', 'FontSize', 12);
grid on;

%% Down Conversion and LPF
y_I = rx_audio.*cos(2*pi*fc*t_rx);
y_Q = -rx_audio.*sin(2*pi*fc*t_rx);

% Low Pass Filter
[B,A] = butter(8,0.1);
[H,F] = freqz(B,A,1024,Fsamp);

r_I = filter(B,A,y_I);
r_Q = filter(B,A,y_Q);
r_bb = r_I + 1j*r_Q;

%PLOT 2 Baseband signal amplitude
figure;
plot(abs(r_bb), 'LineWidth', 0.5);
title('Received Baseband Signal Magnitude', 'FontSize', 14);
xlabel('Sample index (downsampled)', 'FontSize', 12);
ylabel('|r_{bb}|', 'FontSize', 12);
grid on;

% Downsample
r_bb = r_bb(1:L:end);

%% Synchronization and Removal of Cyclic Prefix
half_pilot = Nsc/2;
last_index = (length(r_bb)-(2*half_pilot))+1;
correlation = zeros(1,last_index);

% this loop is sliding the window to calculate correlation for
% different side by side half pilot
% the output is vector correlation that represent the index and its correlation
for idx = 1:last_index
    first_idx  = idx : idx+half_pilot-1;
    second_idx = idx+half_pilot : idx+2*half_pilot-1;

    first_pilot = r_bb(first_idx);
    second_pilot = r_bb(second_idx);

    num = abs(sum(first_pilot .* conj(second_pilot)));
    denum = sqrt(sum(abs(first_pilot).^2) * sum(abs(second_pilot).^2));

    correlation(idx) = num / denum;
end

% detemining the periodic pilot (can be tuned)
[~, max_corr_idx] = max(correlation);
periodic_pilot = max_corr_idx;

figure;
plot(correlation, 'LineWidth', 1.2, 'HandleVisibility', 'off'); hold on;
xline(periodic_pilot, '--r', 'LineWidth', 2.5);
xlabel('Time (s)', 'FontSize',20);
ylabel('Normalized Correlation', 'FontSize',20);
legend('Highest correlation', 'FontSize',25);
xlim([1 5505]) 
set(gca, 'FontSize', 20, 'LineWidth', 1.5);

% FFT Using Coarse Timing
timing_offset = 5;
fprintf('Timing offset = %d\n', timing_offset);
fft_start_first = periodic_pilot - timing_offset; % assume the coarse timing already in CP
ofdm_symbol_len = Nsc + Ncp; %148
num_symbols = floor((length(r_bb)-fft_start_first+1) / ofdm_symbol_len); % rough number after periodic pilot
fft_output = zeros(Nsc,num_symbols); % Row: subcarrier idx | Col: OFDM symbol idx: 

for i = 1:num_symbols
    fft_start = fft_start_first + ((i-1)*ofdm_symbol_len);
    fft_input = r_bb(fft_start:fft_start+Nsc-1);
    fft_output(:,i) = fft(fft_input, Nsc);
end

%% Channel Estimation from Preamble
% Preamble
rng(100, 'twister');
P = sign(randn(1, Nsc/2));
X_pilot = zeros(Nsc,1);
X_pilot(1:2:end) = 2*P.';

% Channel Estimation from the pilot
Y_pilot = fft_output(:,1);
active_pilot_idx = 1:2:Nsc; % since only odd (matlab) has pilot
symbol_idx = (1:Nsc).';

H_estimated = zeros(Nsc,1);
H_estimated(active_pilot_idx) = Y_pilot(active_pilot_idx) ./ X_pilot(active_pilot_idx);

% Interpolate missing subcarriers
H_estimated_interpolated = interp1(active_pilot_idx.', ...
                H_estimated(active_pilot_idx), ...
                symbol_idx, ...
                'linear', ...
                'extrap');

% convert to time-domain impulse response for periodic pilot evaluation
h_estimated = ifft(H_estimated_interpolated, Nsc);

figure;
stem(abs(h_estimated), 'filled');
grid on;
xlabel('Delay sample');
ylabel('|h[n]|');
title('Estimated Channel Impulse Response');
xline(Ncp, '--r', 'CP length');

%% Header Decoding
% expected length bit message length
if decoding_active == true
    memory_length = 5;
    N_coded_header_bits = 2*(header_length + memory_length);
else
    N_coded_header_bits = header_length;
end

% message-legth symbol length
N_header_mod_symbols = N_coded_header_bits / bits_per_symbol;

Y_header = fft_output(:,2);

% Initial Equalization using channel estimated from preamble
eq_header = Y_header ./ H_estimated_interpolated;

% continuous pilot correction for header
rx_header_pilot = eq_header(pilot_idx);
tx_header_pilot = pilot_value.';

pilot_ratio = rx_header_pilot ./ tx_header_pilot;
phase_error = angle(mean(pilot_ratio));

eq_header = eq_header * exp(-1j*phase_error);

% Extract only header data subcarriers
header_symbols_rx = eq_header(data_idx);

% Remove header OFDM padding
header_symbols_rx = header_symbols_rx(1:N_header_mod_symbols);

% QPSK demodulation
coded_header_bits_rx = qpsk_modulation.demodulate(header_symbols_rx);

% Channel decode header
header_bits_rx = channel_coding.decode(coded_header_bits_rx, header_length, decoding_active);

% Convert 32-bit header to payload length
Nbits_from_header = bi2de(header_bits_rx, 'left-msb');

fprintf('Decoded payload bit length from header: %d bits\n', Nbits_from_header);

% OFDM data symbol length from pilot
if decoding_active == true
    memory_length = 5;
    N_coded_payload_bits = 2*(Nbits_from_header + memory_length);
else
    N_coded_payload_bits = Nbits_from_header;
end

N_payload_mod_symbols = N_coded_payload_bits / bits_per_symbol;
Nofdm_data_symbols_from_header = ceil(N_payload_mod_symbols / N_data);

fprintf('Estimated ofdm payload length from header: %d bits\n', Nofdm_data_symbols_from_header);

%% Equalization and Data Decoding
% OFDM data symbol
Y_data = fft_output(:, 3:end);
Y_data = Y_data(:, 1:Nofdm_data_symbols_from_header);

fprintf('Recovered data OFDM symbols: %d\n', size(Y_data,2));

% Initial Equalization using channel estimated from preamble
eq_data = Y_data ./ repmat(H_estimated_interpolated, 1, size(Y_data,2));

%PLOT - Constellation before estimation
figure;
Yeq_raw = eq_data(:,1);
scatter(real(Yeq_raw(data_idx)), imag(Yeq_raw(data_idx)), 20,'b','filled');
hold on;
scatter(real(Yeq_raw(pilot_idx)), imag(Yeq_raw(pilot_idx)),100,'r','x','LineWidth',2);
title('Constellation BEFORE Pilot Phase Correction','FontSize',14);
xlabel('In-phase','FontSize',12); ylabel('Quadrature','FontSize',12);
legend('Data','Pilots','FontSize',11); grid on; axis equal;

% continuous pilot correction
pilot_ratio_matrix = zeros(4, size(eq_data,2));

for m = 1: size(eq_data,2)
    rx_pilot = eq_data(pilot_idx, m);
    tx_pilot = pilot_value.';

    pilot_ratio_matrix = rx_pilot ./ tx_pilot;
    phase_error = angle(mean(pilot_ratio_matrix)); % average error per symbol
    pilot_ratio_matrix(:, m) = pilot_ratio_matrix; % matrix for saving correction freq and time
    
    % Equalization
    eq_data(:,m) = eq_data(:,m) * exp(-1j*phase_error);
end

%PLOT - Pilot phase eveolution over time
figure;
plot(pilot_ratio_matrix.','LineWidth',1.5);
title('Continuous Pilot Phase Evolution Over Time','FontSize',14);
xlabel('OFDM Symbol Index','FontSize',12);
ylabel('Phase (radians)','FontSize',12);
legend(arrayfun(@(x) sprintf('Pilot sc %d',pilot_idx(x)), ...
    1:length(pilot_idx),'UniformOutput',false),'FontSize',10);
yline(0,'k--','Ideal = 0');
grid on;

%PLOT - Pilot phase vs frequency
figure;
phase_vs_freq = angle(eq_data(pilot_idx,1)./pilot_value.');
plot(pilot_idx, phase_vs_freq,'ro-','LineWidth',2,'MarkerSize',8);
title('Pilot Phase vs Subcarrier Index','FontSize',14);
xlabel('Subcarrier Index','FontSize',12);
ylabel('Phase (radians)','FontSize',12);
yline(0,'k--','Ideal = 0'); grid on;
xlim([0 Nsc+1]);

%% Extract only payload data subcarriers
data_decoded_symbols_matrix = eq_data(data_idx, :);
data_decoded_symbols = reshape(data_decoded_symbols_matrix, 1, []); % Parallel to serial

fprintf('Available QPSK symbols after extraction: %d\n', length(data_decoded_symbols));
fprintf('Expected QPSK symbols: %d\n', N_mod_symbols);

% Remove padded QPSK symbols from transmitter
data_decoded_symbols = data_decoded_symbols(1:N_mod_symbols);

%% QPSK Demodulation
bits_decoded = qpsk_modulation.demodulate(data_decoded_symbols);

%% Viterbi Decoding
data_decoded_bits = channel_coding.decode(bits_decoded, Nbits_from_header, decoding_active);

%% BER Calculation
num_errors = sum(bits ~= data_decoded_bits);
BER = num_errors / Nbits_from_header;

fprintf('Number of bit errors: %d\n', num_errors);
fprintf('BER: %.6f\n', BER);

%% Bit Rate
fprintf('bit rate: %.2f kbps\n', bit_rate/1000);

%% Reconstruct file from decoded bits
% Read 32-bit header to get file length
header_bits = data_decoded_bits(1:32);
file_len_rx = bi2de(double(header_bits), 'left-msb');
fprintf('Decoded file length: %d bytes\n', file_len_rx);

% Extract file data bits (after header)
file_bits = data_decoded_bits(33 : 33 + file_len_rx*8 - 1);

% Convert bits back to bytes
bytes    = reshape(double(file_bits), 8, []).';
file_out = uint8(bi2de(bytes, 'left-msb'));

% Write to output file
fid = fopen('output.txt', 'w');
fwrite(fid, file_out, 'uint8');
fclose(fid);
fprintf('File written to output.txt\n');

% Compare with original
fid  = fopen('input.txt', 'r');
orig = fread(fid, '*char').';
fclose(fid);
m = char(file_out.');
if strcmp(m, orig)
    fprintf('*** SUCCESS - File transferred correctly ***\n');
else
    fprintf('*** File has errors ***\n');
end


%% Constellation Plot
figure;
plot(real(data_decoded_symbols), imag(data_decoded_symbols), '.');
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title('Equalized Data Constellation');
axis equal;

delete('recorded_ofdm.wav');
delete('tx_workspace.mat');

