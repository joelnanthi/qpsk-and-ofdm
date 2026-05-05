clc;
clear;
close all;

Fs = 44100;
T_samp = 1 / Fs;
fc = 4000;
Ts = 2.2676e-3;
Ns = round(Fs*Ts);

%load data
data = load('J:\PWC\Assignment1\Signals_task2/Signal1.mat');

R = data.R;
if isfield(data, 't') %t not recognised for signal1
    t = data.t;
else
    t = (0:length(R)-1)*T_samp; %1/44100
end

R = R(:);
t = t(:);

%figure;
%plot(t,R);

t_pulse = (0:Ns-1)/Fs;
g = sin(pi*t_pulse/Ts);
g = g(:);
g = g / norm(g);
g_fliped = flipud(g);

yI = R.*cos(2*pi*fc*t);
yQ = -R.*sin(2*pi*fc*t);

mI_t = conv(yI, g_fliped);
mQ_t = conv(yQ, g_fliped);
m_t = mI_t + 1j*mQ_t;

%envelope
e_mt = sqrt(mI_t.^ 2 + mQ_t.^2);
t_env = (0:length(e_mt)-1)/Fs; %new length after convolution

%window to slide
window_energy = movmean(e_mt, Ns);

%threshold(trial and error - different for each signal - 0.3 to 04 range)

threshold = 0.4 * max(window_energy);

%find max energy region
idx = find(window_energy > threshold);

first_peak = idx(1);
last_peak = idx(end);

%plot
figure;
plot(t_env, window_energy, 'LineWidth', 1.5); hold on;
yline(threshold, '--r', 'Threshold');
plot(t_env(idx), window_energy(idx), 'r.');
title('Sliding Window Energy (Synchronization)');
xlabel('Time (s)');
ylabel('Energy');
grid on;

%Tsample search
best_offset = 1;
best_score = -inf;
average_score = zeros(1, Ns);

for offset = 1:Ns
    idx = (first_peak+offset-1):Ns:last_peak; %should start from zero offset

    if (isempty(idx))
        continue;
    end

    temp = m_t(idx);
    evaluation_score = mean(abs(temp));
    average_score(offset) = evaluation_score;

    if (evaluation_score > best_score)
        best_score = evaluation_score;
        best_offset = offset;
    end
end
Tsamp_idx = first_peak + best_offset - 1;

%sampling
rx_idx = Tsamp_idx:Ns:last_peak;
num_symbols = length(rx_idx);
rx_symbols = m_t(rx_idx);

%sampled point at transmitted data
figure;
plot(t_env, abs(m_t), 'LineWidth', 1.2); hold on; grid on;
plot(t_env(rx_idx), abs(rx_symbols), 'ro');
xlabel('Time (s)');
ylabel('|m(t)|');
title('Sampled points on matched filter output');

%check
figure;
plot(real(rx_symbols), imag(rx_symbols), '.');
grid on;
axis equal;
xlabel('In-phase');
ylabel('Quadrature');
title('Received sampled symbols');

%Pilot Estimation
pilot = 2 + 2i;
channel_effect = rx_symbols(1) / pilot; %or rx(end), or avg

%Channel Equalization
rx_eq = rx_symbols / channel_effect;

figure;
plot(real(rx_eq), imag(rx_eq), '.');
grid on;
axis equal;
xlabel('In-phase');
ylabel('Quadrature');
title('Equalized symbols');

%pilot Removal
data_syms = rx_eq(2:end-1);

figure;
plot(real(data_syms), imag(data_syms), '.');
grid on;
axis equal;
xlabel('In-phase');
ylabel('Quadrature');
title('Payload symbols after equalization and pilot removal');

%Decoding
Nsymbols = length(data_syms);
detected_symbols = zeros(1,Nsymbols);
detected_bits = zeros(1,2*Nsymbols);

idx = 1;
for i = 1:Nsymbols
    I = real(data_syms(i));
    Q = imag(data_syms(i));
    if (I > 0 && Q > 0)
        detected_symbols(i) = 1 + 1j;
        detected_bits(idx:idx+1) = [0 0];
    elseif (I < 0 && Q > 0)
        detected_symbols(i) = -1 + 1j;
        detected_bits(idx:idx+1) = [1 0];
    elseif (I < 0 && Q < 0)
        detected_symbols(i) = -1 - 1j;
        detected_bits(idx:idx+1) = [1 1];
    elseif (I > 0 && Q < 0)
        detected_symbols(i) = 1 - 1j;
        detected_bits(idx:idx+1) = [0 1];
    end
    idx = idx + 2;
end

%ASCII mapping
wh = 2.^(6:-1:0);
m = char(detected_bits(1:7) * wh');
for l = 2:floor(length(detected_bits)/7)
    m = [m char(detected_bits(7*(l-1)+1 : 7*l) * wh')];
end

%decoded message
disp(m);
