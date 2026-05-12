function tx_audio = tx_packet_ofdm(packet_bits, params)
% Transmit one packet through OFDM physical layer
% packet_bits: row vector [packet_id(16) | data(500) | crc(16)] = 532 bits
% Returns tx_audio ready to play

Fs          = params.Fs;
fc          = params.fc;
N_sc        = params.N_sc;
N_cp        = params.N_cp;
L           = params.L;
pilot_idx   = params.pilot_idx;
pilot_value = params.pilot_value;
data_idx    = params.data_idx;
N_data      = params.N_data;
decoding_active = params.decoding_active;
header_length   = params.header_length;
ofdm_symbol_len = N_sc + N_cp;

%% Channel coding of packet bits
[coded_bits, N_coded_bits] = channel_coding.encode(packet_bits, decoding_active);

%% QPSK modulation
[a, N_mod_symbols] = qpsk_modulation.modulate(coded_bits, params);

%% OFDM data symbols
Nofdm_data_symbols = ceil(N_mod_symbols / N_data);
Npad    = Nofdm_data_symbols*N_data - length(a);
a_pad   = [a  zeros(1, Npad)];

ofdm_symbols = zeros(Nofdm_data_symbols, N_sc);
for i = 1:Nofdm_data_symbols
    sym = zeros(1, N_sc);
    sym(data_idx)  = a_pad((i-1)*N_data+1 : i*N_data);
    sym(pilot_idx) = pilot_value;
    ofdm_symbols(i,:) = sym;
end

x    = ifft(ofdm_symbols, N_sc, 2);
x_cp = [x(:,end-N_cp+1:end)  x];
tx_data_bb = reshape(x_cp.', 1, []);

%% Header OFDM symbol (contains N_coded_bits so receiver knows length)
header_bits = de2bi(N_coded_bits, header_length, 'left-msb');
[coded_header, ~] = channel_coding.encode(header_bits, decoding_active);
[header_sym, N_hdr_sym] = qpsk_modulation.modulate(coded_header, params);

header_sym_pad = [header_sym  zeros(1, N_data - N_hdr_sym)];
hdr_ofdm = zeros(1, N_sc);
hdr_ofdm(data_idx)  = header_sym_pad;
hdr_ofdm(pilot_idx) = pilot_value;

x_hdr    = ifft(hdr_ofdm, N_sc);
x_hdr_cp = [x_hdr(end-N_cp+1:end)  x_hdr];

%% Preamble
rng(100, 'twister');
P = sign(randn(1, N_sc/2));
X_pre = zeros(1, N_sc);
X_pre(1:2:end) = 2*P;
x_pre    = ifft(X_pre, N_sc);
x_pre_cp = [x_pre(end-N_cp+1:end)  x_pre];

%% Combine: preamble + header + data + tail guard
tx_bb = [x_pre_cp,  x_hdr_cp,  tx_data_bb, ...
         zeros(1, params.N_tail_guard * ofdm_symbol_len)];

%% Upsample and upconvert
tx_bb_up = interp(tx_bb, L);
n        = 0:length(tx_bb_up)-1;
t        = n / Fs;

I = real(tx_bb_up);
Q = imag(tx_bb_up);
tx_audio = I.*cos(2*pi*fc*t) - Q.*sin(2*pi*fc*t);
tx_audio = tx_audio / max(abs(tx_audio)) * 0.8;
end