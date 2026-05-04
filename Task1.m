clc;
clear all;
close all;

Fs = 44100;
fc = 4000;
Ts = 2.2676e-3;
Ns = round(Fs*Ts);
mu = 0.00022676;
alpha = 6/18; %Groupnum/18

M = 4;
K = 2;
d2min = 2; %QPSK
Nb = 1e6;
EbN0_dB = 1:2:9;
BER = zeros(size(EbN0_dB));

t = (0:Ns-1)/Fs;
g = sin(pi*t/Ts);
g = g / norm(g); %normalise

g_fliped = flipud(g);

delay = round(mu*Fs);

delay_total = length(g) - 1;
start = delay_total + 1;

for ii = 1:length(EbN0_dB)

    bits = randi([0 1], Nb, 1);

    %qpsk
    bits2 = reshape(bits, 2, []).';
    symbols = zeros(size(bits2,1),1);

    for k = 1:length(symbols)

        b1 = bits2(k,1);
        b2 = bits2(k,2);
        if b1==0 && b2==0
            symbols(k) = 1+1i;
        elseif b1==1 && b2==0
            symbols(k) = -1+1i;
        elseif b1==1 && b2==1
            symbols(k) = -1-1i;
        elseif b1==0 && b2==1
            symbols(k) = 1-1i;
        end
    end

    x = zeros(length(symbols)*Ns,1);
    x(1:Ns:end) = symbols;

    s_bb = conv(x, g);

    t_full = (0:length(s_bb)-1)/Fs;
    I = real(s_bb);
    Q = imag(s_bb);
    s = I .* cos(2*pi*fc*t_full.') - Q .* sin(2*pi*fc*t_full.'); %bb to bp at 4000hz

    %multipath channel
    s_delayed = [zeros(delay,1); s(1:end-delay)];

    r = sqrt(1-alpha^2)*s + alpha*s_delayed;

    E_total = sum(s.^2)/Fs;
    Eb = E_total / Nb;

    %add noise
    EbN0 = 10^(EbN0_dB(ii)/10);
    N0 = Eb / EbN0; %to get N0

    noise = sqrt(N0*Fs/2)*randn(size(r)); %awgn
    r = r + noise;

    %downconvert - receiver
    r_bb = r .* exp(-1i*2*pi*fc*t_full.'); %bp to bb

    y = filter(g_fliped, 1, r_bb);

    num_symbols = length(symbols);

    y_sampled = y(start:Ns:start + (num_symbols-1)*Ns);

    %detection
    detected_bits = zeros(2*num_symbols,1);

    for k = 1:num_symbols
        sym = y_sampled(k);
        I = real(sym);
        Q = imag(sym);
        %decision %check
        if (I > 0 && Q > 0)
            % 1 + j means 00
            detected_bits(2*k-1) = 0;
            detected_bits(2*k) = 0;
        elseif (I < 0 && Q > 0)
            % -1 + j means 10
            detected_bits(2*k-1) = 1;
            detected_bits(2*k) = 0;
        elseif (I < 0 && Q < 0)
            % -1 - j means 11
             detected_bits(2*k-1) = 1;
             detected_bits(2*k) = 1;
        elseif (I >0 && Q < 0)
            % 1 - j means 01
            detected_bits(2*k-1) = 0;
            detected_bits(2*k) = 1;
        end
    end
    %ber
    bits_new = bits(1:length(detected_bits));
    BER(ii) = sum(bits_new ~= detected_bits) / length(bits_new);

    %fprintf('Eb/N0 = %d dB, BER = %e\n', EbN0_dB(ii), BER(ii));

end

%theoretical BER
BER_theory = qfunc(sqrt(d2min*10.^(EbN0_dB/10)));

figure;
semilogy(EbN0_dB, BER, 'o-','LineWidth',2);
hold on;
semilogy(EbN0_dB, BER_theory, '--','LineWidth',2);
grid on;
xlabel('E_b/N_0 (dB)','FontSize',16);
ylabel('BER','FontSize',16);
legend('Simulation','Theory','FontSize',14);
title('QPSK BER over Multipath Channel','FontSize',16);
