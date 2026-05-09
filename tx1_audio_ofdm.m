clear; clc; close all;

%% Parameter
params = ofdm_params();

Fsamp = params.Fsamp;
Tsamp = params.Tsamp;
fc = params.fc;
Nsc = params.Nsc;
Ncp = params.Ncp;

T_ofdm = params.T_ofdm;
F_ofdm = params.F_ofdm;
T_cp = params.T_cp;
T_sym = params.T_sym;
L = params.L;
delta_f = params.delta_f;

decoding_active = params.decoding_active;

header_length = params.header_length;

%% Source coding
Nbits = 2e4;
        
if mod(Nbits,2) ~= 0
    error("The number of bits must be even");
end

%bits = randi([0 1], 1, Nbits);
fid = fopen('input.txt','r');
fileData = fread(fid,'*uint8');
fclose(fid);
file_len = length(fileData);

%convert bytes to bits
data_bits = double(reshape(de2bi(fileData, 8, 'left-msb').', 1, []));

% 32-bit header with file length for receiver to know how many bytes to recover
header = double(de2bi(file_len, 32, 'left-msb'));
header = header(:).';

% Complete bit stream: header + file data
bits = [header  data_bits];
Nbits = length(bits);

fprintf('File: %d bytes = %d bits\n', file_len, Nbits);

%% Channel Coding
[coded_bits, N_coded_bits] = channel_coding.encode(bits, decoding_active);

%% QPSK Modulation
[a, N_mod_symbols] = qpsk_modulation.modulate(coded_bits, params);

%% ofdm symbol structure
N_active = params.N_active; % data include continous pilot
N_dc = params.N_dc;
N_outer_zero =  Nsc-N_active-N_dc; % 24
N_zero_each_side = N_outer_zero/2; % 12
pilot_idx = params.pilot_idx;
pilot_value = params.pilot_value;
N_data = N_active-length(pilot_idx); % real data without pilot

left_guard_idx = params.left_guard_idx;
dc_idx = params.dc_idx;
right_guard_idx = params.right_guard_idx;

% active index for data
reserved_idx = [left_guard_idx dc_idx right_guard_idx pilot_idx];
reserved_idx = sort(reserved_idx);

data_idx = setdiff(1:Nsc, reserved_idx);
data_idx = sort(data_idx);
Nofdm_data_symbols = ceil(N_mod_symbols/N_data); % number of OFDM symbols needed to carry the QPSK data symbols
fprintf('Num of OFDM Symbol: %d\n', Nofdm_data_symbols);

% padded the last remain bit with 0 one last symbol
Npad = Nofdm_data_symbols*N_data - length(a);
a_padded = [a zeros(1, Npad)];

% mapping to one OFDM symbols
ofdm_symbols = zeros(Nofdm_data_symbols, Nsc);
for i = 1:Nofdm_data_symbols
    one_ofdm_symbol = zeros(1,Nsc);

    start_idx = (i-1)*N_data + 1;
    end_idx   = i*N_data;

    data_for_this_ofdm = a_padded(start_idx:end_idx);

    one_ofdm_symbol(data_idx) = data_for_this_ofdm;
    one_ofdm_symbol(pilot_idx) = pilot_value;

    ofdm_symbols(i,:) = one_ofdm_symbol;
end

%% Payload IFFT and CP
x = ifft(ofdm_symbols, Nsc, 2);
x_cp = [x(:, end-Ncp+1:end) x];

% Paralel to serial
tx_data_bb = reshape(x_cp.', 1, []);

%% Header Generation
header_bits = de2bi(Nbits, header_length, 'left-msb'); %left most big value

% Channel Coding
[coded_header_bits, N_coded_header_bits] = channel_coding.encode(header_bits, decoding_active);

% QPSK Modulation
[header_symbols, N_header_mod_symbols] = qpsk_modulation.modulate(coded_header_bits, params);

% padded the last remain bit with 0 one last symbol
N_header_pad = N_data - N_header_mod_symbols;
header_symbols_padded = [header_symbols zeros(1, N_header_pad)];

% make one OFDM symbol
header_ofdm_symbol = zeros(1, Nsc);
header_ofdm_symbol(data_idx) = header_symbols_padded;
header_ofdm_symbol(pilot_idx) = pilot_value;

% Header IFFT and CP
x_header = ifft(header_ofdm_symbol, Nsc);
x_header_cp = [x_header(end-Ncp+1:end) x_header];

%% Preamble Generation
X_preamble = zeros(1, Nsc);
rng(100, 'twister');
P = sign(randn(1, Nsc/2));
X_preamble(1:2:end) = 2 * P;

% Preamble IFFT and CP
x_preamble = ifft(X_preamble, Nsc);
x_preamble_cp = [x_preamble(end-Ncp+1:end) x_preamble];

%% combine preamble, header, data, and tail guard
N_tail_guard = params.N_tail_guard;
ofdm_symbol_len = Nsc + Ncp;
tx_bb = [x_preamble_cp, x_header_cp, tx_data_bb, zeros(1, N_tail_guard*ofdm_symbol_len)];

%% Upsampling
tx_bb_up = interp(tx_bb, L);
%[tx_bb_zero, tx_bb_up_vis] = plot_upsampling_visualization(tx_bb, L, 30);

n = 0:length(tx_bb_up)-1;
t = n/Fsamp;

%% Upconversion
I = real(tx_bb_up);
Q = imag(tx_bb_up);

tx_audio = I.*cos(2*pi*fc*t) - Q.*sin(2*pi*fc*t);
% Normalize 
tx_audio = tx_audio / max(abs(tx_audio)) * 0.8;

%% Bit Calculation
N_total_ofdm_symbols = 1 + 1 + Nofdm_data_symbols + N_tail_guard; % preamble + header + data + tail
T_total_theory = N_total_ofdm_symbols * T_sym;
T_total_audio = length(tx_audio)/Fsamp;

bit_rate = Nbits/T_total_audio;

fprintf('Theoretical frame duration: %.3f s\n', T_total_theory);
fprintf('Audio frame duration: %.3f s\n', T_total_audio);
fprintf('bit rate: %.2f kbps\n', bit_rate/1000);

%% save data
save('tx_workspace.mat', ...
     'tx_audio', ...
     'bits', ...
     'N_coded_bits', ...
     'N_mod_symbols', ...
     'bit_rate');

audiowrite('ofdm_tx.wav', tx_audio, Fsamp);