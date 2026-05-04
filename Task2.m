clc;
clear;
close all;

%% Parameters (same as Task 1)
Fs = 44100;
fc = 4000;
Ts = 2.2676e-3;
Ns = round(Fs*Ts);

t_pulse = (0:Ns-1)/Fs;
g = sin(pi*t_pulse/Ts);   % half-sine pulse
h = flipud(g);            % matched filter

%% Loop over all signals
for file_idx = 1:13

    filename = ['Signal', num2str(file_idx), '.mat'];
    load(filename);   % assumes variable inside is r

    fprintf('\nProcessing %s...\n', filename);

    %% Ensure column vector
    if size(r,2) > 1
        r = r.';
    end

    %% STEP 1: Downconversion (only if real signal)
    if isreal(r)
        t = (0:length(r)-1)/Fs;
        r_bb = r .* exp(-1i*2*pi*fc*t.');
    else
        r_bb = r;
    end

    %% STEP 2: Matched filter
    y = filter(h, 1, r_bb);

    %% STEP 3: Synchronization (find best offset automatically)
    energies = zeros(Ns,1);

    for offset = 0:Ns-1
        y_test = y(offset+1:Ns:end);
        energies(offset+1) = mean(abs(y_test).^2);
    end

    [~, best_offset] = max(energies);
    best_offset = best_offset - 1;

    fprintf('Best offset = %d\n', best_offset);

    %% STEP 4: Sample at correct timing
    y_sampled = y(best_offset+1:Ns:end);

    %% Optional: Normalize (helps constellation)
    y_sampled = y_sampled / max(abs(y_sampled));

    %% Plot constellation
    figure;
    scatter(real(y_sampled(1:2000)), imag(y_sampled(1:2000)), '.');
    title(['Constellation - Signal ', num2str(file_idx)]);
    xlabel('In-phase');
    ylabel('Quadrature');
    grid on;

    %% STEP 5: Detection (QPSK)
    num_symbols = length(y_sampled);
    detected_bits = zeros(2*num_symbols,1);

    for k = 1:num_symbols

        I = real(y_sampled(k));
        Q = imag(y_sampled(k));

        if (I > 0 && Q > 0)
            detected_bits(2*k-1) = 0;
            detected_bits(2*k)   = 0;

        elseif (I < 0 && Q > 0)
            detected_bits(2*k-1) = 1;
            detected_bits(2*k)   = 0;

        elseif (I < 0 && Q < 0)
            detected_bits(2*k-1) = 1;
            detected_bits(2*k)   = 1;

        elseif (I > 0 && Q < 0)
            detected_bits(2*k-1) = 0;
            detected_bits(2*k)   = 1;
        end
    end

    %% STEP 6: Convert bits → text
    % Ensure multiple of 8
    num_bits = floor(length(detected_bits)/8)*8;
    detected_bits = detected_bits(1:num_bits);

    bits_matrix = reshape(detected_bits, 8, []).';
    message = char(bin2dec(num2str(bits_matrix)));

    fprintf('Decoded message:\n');
    disp(message);

end