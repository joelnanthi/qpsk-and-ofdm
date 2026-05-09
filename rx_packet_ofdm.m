%function [decoded_bits, crc_ok, packet_id] = rx_packet_ofdm(rx_audio, params)
function [decoded_bits, crc_ok, packet_id, eq_symbols, pilot_phases] = rx_packet_ofdm(rx_audio, params)
% Receive and decode one packet
% rx_audio: recorded row vector
% Returns decoded info bits, CRC pass/fail, and packet ID

Fsamp           = params.Fsamp;
fc              = params.fc;
Nsc             = params.Nsc;
Ncp             = params.Ncp;
L               = params.L;
pilot_idx       = params.pilot_idx;
pilot_value     = params.pilot_value;
data_idx        = params.data_idx;
N_data          = params.N_data;
decoding_active = params.decoding_active;
header_length   = params.header_length;
bits_per_symbol = params.bits_per_symbol;

% Default return values in case of failure
decoded_bits = zeros(1, params.packet_size);
crc_ok       = false;
packet_id    = -1;
eq_symbols   = [];   % ADD THIS LINE
pilot_phases = [];   % ADD THIS LINE

%% Downconvert and LPF
n_rx = 0:length(rx_audio)-1;
t_rx = n_rx / Fsamp;

y_I = rx_audio .* cos(2*pi*fc*t_rx);
y_Q = -rx_audio .* sin(2*pi*fc*t_rx);

[B,A] = butter(8, 0.1);
r_I   = filter(B, A, y_I);
r_Q   = filter(B, A, y_Q);
r_bb  = r_I + 1j*r_Q;
r_bb  = r_bb(1:L:end);

%% CP correlation sync
half_pilot = Nsc/2;
last_idx   = length(r_bb) - 2*half_pilot + 1;
if last_idx < 1
    return;
end

correlation = zeros(1, last_idx);
for idx = 1:last_idx
    a_w = r_bb(idx:idx+half_pilot-1);
    b_w = r_bb(idx+half_pilot:idx+2*half_pilot-1);
    num   = abs(sum(a_w.*conj(b_w)));
    denum = sqrt(sum(abs(a_w).^2)*sum(abs(b_w).^2));
    if denum > 0
        correlation(idx) = num/denum;
    end
end

[~, max_corr_idx] = max(correlation);
timing_offset   = 5;
fft_start_first = max(max_corr_idx - timing_offset, 1);
ofdm_sym_len    = Nsc + Ncp;

num_symbols = floor((length(r_bb) - fft_start_first + 1) / ofdm_sym_len);
if num_symbols < 3
    return;   % need at least preamble + header + 1 data symbol
end

fft_output = zeros(Nsc, num_symbols);
for i = 1:num_symbols
    fs = fft_start_first + (i-1)*ofdm_sym_len;
    if fs+Nsc-1 > length(r_bb), break; end
    fft_output(:,i) = fft(r_bb(fs:fs+Nsc-1), Nsc);
end

%% Channel estimation from preamble (column 1)
rng(100, 'twister');
P_pre = sign(randn(1, Nsc/2));
X_pilot = zeros(Nsc, 1);
X_pilot(1:2:end) = 2*P_pre.';

Y_pilot     = fft_output(:,1);
active_idx  = (1:2:Nsc).';
H_est       = zeros(Nsc, 1);
H_est(active_idx) = Y_pilot(active_idx) ./ X_pilot(active_idx);
H_est = interp1(active_idx, H_est(active_idx), (1:Nsc).', 'linear', 'extrap');

%% Header decoding (column 2)
if decoding_active
    N_coded_hdr = 2*(header_length + 5);
else
    N_coded_hdr = header_length;
end
N_hdr_sym = N_coded_hdr / bits_per_symbol;

eq_hdr = fft_output(:,2) ./ H_est;
ph_hdr = angle(mean(eq_hdr(pilot_idx) ./ pilot_value.'));
eq_hdr = eq_hdr * exp(-1j*ph_hdr);

hdr_syms_rx = eq_hdr(data_idx);
hdr_syms_rx = hdr_syms_rx(1:N_hdr_sym);
coded_hdr_bits = qpsk_modulation.demodulate(hdr_syms_rx.');
hdr_bits = channel_coding.decode(coded_hdr_bits, header_length, decoding_active);

N_coded_payload = bi2de(hdr_bits, 'left-msb');
if N_coded_payload <= 0 || N_coded_payload > 1e5
    return;
end

N_payload_sym    = N_coded_payload / bits_per_symbol;
Nofdm_from_hdr   = ceil(N_payload_sym / N_data);

%% Data decoding (columns 3 to end)
avail = size(fft_output, 2) - 2;
Nofdm_use = min(Nofdm_from_hdr, avail);
if Nofdm_use < 1
    return;
end

Y_data  = fft_output(:, 3:3+Nofdm_use-1);
eq_data = Y_data ./ repmat(H_est, 1, size(Y_data,2));

for m = 1:size(eq_data,2)
    ph = angle(mean(eq_data(pilot_idx,m) ./ pilot_value.'));
    eq_data(:,m) = eq_data(:,m) * exp(-1j*ph);
end

data_syms = eq_data(data_idx, :);
data_syms = reshape(data_syms, 1, []);

N_needed = N_coded_payload / bits_per_symbol;
if length(data_syms) < N_needed
    return;
end
data_syms = data_syms(1:N_needed);

%% QPSK demod + Viterbi
bits_demod  = qpsk_modulation.demodulate(data_syms);

% N_info = packet_id_length + packet_size + crc_length
N_info = params.packet_id_length + params.packet_size + params.crc_length;
all_bits = channel_coding.decode(bits_demod, N_info, decoding_active);

%% Extract packet ID, data, CRC and verify
id_bits   = all_bits(1 : params.packet_id_length);
data_rx   = all_bits(params.packet_id_length+1 : params.packet_id_length+params.packet_size);
crc_rx    = all_bits(params.packet_id_length+params.packet_size+1 : end);

packet_id    = bi2de(double(id_bits), 'left-msb');
expected_crc = compute_crc16([id_bits  data_rx]);
crc_ok       = isequal(double(crc_rx), double(expected_crc));
decoded_bits = data_rx;

%% Return equalised symbols and pilot phases for plotting
% eq_symbols: all equalised data subcarrier symbols across all data OFDM symbols
eq_symbols   = data_syms;   % already computed above

% pilot_phases: phase of each pilot across OFDM symbols (4 x Nofdm_use)
pilot_phases_out = zeros(length(pilot_idx), size(eq_data,2));
for m = 1:size(eq_data,2)
    pilot_phases_out(:,m) = angle(eq_data(pilot_idx,m) ./ pilot_value.');
end
pilot_phases = pilot_phases_out;

end