clc;
clear;

fs = 44100;
fc = 4000;
Ts = 2.2676e-3;
Ns = round(fs*Ts); 
mu = 0.00022676;
alpha = 7/18; %groupnum/18

M = 4;
K = 2;
Nb = 5e5;
EbN0_dB = 1:1:9;
BER = zeros(size(EbN0_dB));

t = (0:Ns-1)/fs;
g = sin(pi*t/Ts);
%plot(g);
for ii = 1:length(EbN0_dB)

    bits = randi([0 1], Nb, 1);

    %qpsk
    qpsk_bits = reshape(bits,2,[]).';
    symbols = zeros(size(qpsk_bits,1),1);

    for k = 1:length(symbols)
        b1 = qpsk_bits(k,1);
        b2 = qpsk_bits(k,2);

        if(b1 == 0 && b2 == 0)
            symbols(k) = 1+1i;
        elseif(b1 == 1 && b2 == 0)
            symbols(k) = -1+1i;
        elseif(b1 == 1 && b2==1)
            symbols(k) = -1 - 1i;
        elseif(b1 == 0 && b2 == 1)
            symbols(k) = -1-1i;
        end
    end

    x = zeros(length(symbols)*Ns,1);
    x(1:Ns:end) = symbols; 

    s_bb = conv(x,g);  %continuous time baseband

    t_full = (0:length(s_bb)-1/fs);
    I = real(s_bb);
    Q = imag(s_bb);
    s_t = I .* cos(2*pi*fc*t_full.') - Q .* sin(2*pi*fc*t_full.');
    %s_t = real(s_bb.*exp(1i*2*pi*fc*t_full.'));  %shift to 4000hz

    delay = round(mu*fs);
    s_delayed = [zeros(delay,1); s_t(1:end-delay)];

    r = sqrt(1-alpha^2)*s_t + alpha*s_delayed; %multipath channel
    
    E_total = sum(s_t.^2) * (1/fs);
    Eb = E_total/Nb;
    
    EbN0 = 10^(EbN0_dB(ii)/10); %linear conversion
    N0 = Eb/EbN0;  %why
    noise = sqrt(N0*fs/2)*randn(size(r));
   
    r = r + noise;  %add noise to r

    %receiver
    r_bb = r.*exp(-1i*2*pi*fc*t_full.');   %remove carrier

    y = conv(r_bb, flipud(g));

    num_symbols = length(symbols);

    delay_total = length(g) - 1;
    start = delay_total + 1;

    y_sampled = y(start:Ns:start + (num_symbols-1)*Ns);
    
    detected_bits = zeros(2*num_symbols,1);
    for k = 1:num_symbols
        sym = y_sampled(k);

        detected_bits(2*k-1) = real(sym) < 0;
        detected_bits(2*k)   = imag(sym) < 0;
    end

    bits_trimmed = bits(1:length(detected_bits));
    BER(ii) = sum(bits_trimmed ~= detected_bits) / length(bits_trimmed);

    fprintf('Eb/N0 = %d dB, BER = %e\n', EbN0_dB(ii), BER(ii));
end
d2min = 2;
BER_theory = qfunc(sqrt(d2min*10.^(EbN0_dB/10)));

figure;
semilogy(EbN0_dB, BER, 'o-','LineWidth',2); hold on;
semilogy(EbN0_dB, BER_theory, '--','LineWidth',2);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate');
legend('Simulation','Theory');
title('QPSK BER over Multipath Channel');





