%% task_2.m - Stop-and-Wait ARQ over Audio Channel
clear;
clc; 
close all;

params = ofdm_params();
Fs  = params.Fs;

%%Slot timing 
slotDuration = 3.0;   % seconds per slot
                      % must be longer than longest audio packet
                      % data packet ~0.6s + silences ~0.6s + margin
buffer_time  = 0.2;   % seconds before each TX starts

t0   = tic;           % master timer
slot = 0;

fprintf('Slot duration: %.1f s\n', slotDuration);
fprintf('Starting ARQ with tic/toc timing\n\n');
%Load image and convert to bits 
% Use a small grayscale PNG for best demo effect (e.g. 64x64 pixels)
% Create a test image if you do not have one
if ~exist('input_image.png','file')
    img_test = uint8(reshape(0:255, 16, 16));  % small 16x16 gradient
    imwrite(img_test, 'input_image.png');
end

img_orig = imread('input_image.png');
if size(img_orig,3) == 3
    img_orig = rgb2gray(img_orig);   % use grayscale
end
img_size = size(img_orig);
img_vec  = img_orig(:);   % column vector of uint8 pixel values
all_bits = reshape(de2bi(img_vec, 8, 'left-msb').', 1, []);

fprintf('Image: %dx%d = %d pixels = %d bits\n', ...
    img_size(1), img_size(2), numel(img_orig), length(all_bits));

%%Split into packets 
pkt_size    = params.packet_size;     % 500 bits per packet
num_packets = ceil(length(all_bits) / pkt_size);

% Pad to fill last packet
padded_bits = [all_bits  zeros(1, num_packets*pkt_size - length(all_bits))];
packets     = reshape(padded_bits, pkt_size, num_packets).';
fprintf('Packets: %d x %d bits\n', num_packets, pkt_size);

%ARQ state 
pkt_received    = false(1, num_packets);
received_data   = zeros(num_packets, pkt_size);
retransmit_cnt  = zeros(1, num_packets);
total_tx        = 0;

%%Plot accumulators
all_eq_symbols   = [];        % collect all equalised symbols across packets
all_pilot_phases = [];        % collect pilot phases across packets
per_pkt_ber      = zeros(1, num_packets);   % raw BER per packet
retransmit_log   = zeros(1, num_packets);   % final retransmission count

%%Image display 
fig = figure('Name','Task 2 - ARQ Image Reconstruction','NumberTitle','off');
img_rx = zeros(img_size, 'uint8');
subplot(1,2,1); imshow(img_orig); title('Original Image');
ax2 = subplot(1,2,2); imshow(img_rx); title('Received Image (builds up)');
drawnow;

%%Silence padding 
sil_len = round(0.3 * Fs);   % 0.3 s silence before each packet
silence = zeros(1, sil_len);

%ACK and NACK packets 
%% ACK packet: ID=0, alternating data, valid CRC
ack_id     = de2bi(0, 16, 'left-msb');
ack_data   = repmat([1 0], 1, 250);          % 500 bits
ack_crc    = compute_crc16([ack_id  ack_data]);
ACK_packet = double([ack_id  ack_data  ack_crc]);   % 532 bits

%% NACK packet: ID=1, alternating data, valid CRC
nack_id     = de2bi(1, 16, 'left-msb');
nack_data   = repmat([0 1], 1, 250);         % 500 bits
nack_crc    = compute_crc16([nack_id  nack_data]);
NACK_packet = double([nack_id  nack_data  nack_crc]);  % 532 bits

fprintf('ACK packet:  %d bits, CRC computed\n', length(ACK_packet));
fprintf('NACK packet: %d bits, CRC computed\n', length(NACK_packet));

%MAIN ARQ LOOP 
p = 1;
while p <= num_packets

    %%DATA SLOT (even slot)
    %% Wait for slot start
    while toc(t0) < slot * slotDuration
        pause(0.001);
    end
    fprintf('[t=%.2fs] DATA SLOT %d\n', toc(t0), slot);
    pause(buffer_time);   % give mic time to open before speaker starts

    if retransmit_cnt(p) == 0
        fprintf('[TX] Sending packet %d/%d\n', p, num_packets);
    else
        fprintf('[TX] Retransmitting packet %d (attempt %d)\n', ...
            p, retransmit_cnt(p)+1);
    end

    %% Build packet
    id_bits     = de2bi(p, params.packet_id_length, 'left-msb');
    data_bits   = packets(p,:);
    crc_bits    = compute_crc16([id_bits  data_bits]);
    packet_bits = double([id_bits  data_bits  crc_bits]);

    %% Generate audio
    tx_audio   = tx_packet_ofdm(packet_bits, params);
    play_audio = [silence  tx_audio  silence];
    total_tx   = total_tx + 1;

    %% Record and play data packet
    recObj = audiorecorder(Fs, 16, 1);
    record(recObj);
    pause(0.1);
    player = audioplayer(play_audio, Fs);
    playblocking(player);
    pause(0.3);
    stop(recObj);
    rx_audio = getaudiodata(recObj).';
    fprintf('[t=%.2fs] Data slot done\n', toc(t0));

    slot = slot + 1;   % advance to feedback slot

    %% Decode data packet
    [decoded_bits, crc_ok, rx_id, eq_symbols, ~] = ...
        rx_packet_ofdm(rx_audio, params);

    if ~isempty(eq_symbols)
        all_eq_symbols = [all_eq_symbols  eq_symbols];
    end

    %% Decide ACK or NACK
    if crc_ok && rx_id == p
        fb_packet = ACK_packet;
        fprintf('[RX] CRC OK — will send ACK\n');
    else
        fb_packet = NACK_packet;
        if ~crc_ok
            fprintf('[RX] CRC FAIL — will send NACK\n');
        else
            fprintf('[RX] ID mismatch — will send NACK\n');
        end
    end

    %%FEEDBACK SLOT (odd slot) 
    %% Wait for feedback slot start
    while toc(t0) < slot * slotDuration
        pause(0.001);
    end
    fprintf('[t=%.2fs] FEEDBACK SLOT %d\n', toc(t0), slot);
    pause(buffer_time);

    %% Generate and play ACK/NACK audio
    fb_audio  = tx_packet_ofdm(fb_packet, params);
    play_fb   = [silence  fb_audio  silence];

    rec_fb = audiorecorder(Fs, 16, 1);
    record(rec_fb);
    pause(0.1);
    fb_player = audioplayer(play_fb, Fs);
    playblocking(fb_player);
    pause(0.3);
    stop(rec_fb);
    rx_fb_audio = getaudiodata(rec_fb).';
    fprintf('[t=%.2fs] Feedback slot done\n', toc(t0));

    slot = slot + 1;   % advance to next data slot

    %% Decode ACK/NACK
    [~, fb_crc_ok, fb_id] = rx_packet_ofdm(rx_fb_audio, params);

    if fb_crc_ok
        if fb_id == 0
            fprintf('[TX] ACK decoded (ID=0)\n');
            is_ack = true;
        elseif fb_id == 1
            fprintf('[TX] NACK decoded (ID=1)\n');
            is_ack = false;
        else
            fprintf('[TX] Unknown ID=%d — defaulting NACK\n', fb_id);
            is_ack = false;
        end
    else
        fprintf('[TX] Feedback CRC failed — defaulting NACK\n');
        is_ack = false;
    end

    %% ACK/NACK action
    if is_ack
        pkt_received(p)    = true;
        received_data(p,:) = decoded_bits;
        p = p + 1;

        bits_so_far = reshape(received_data(pkt_received,:).', 1, []);
        n_bytes     = floor(length(bits_so_far)/8);
        if n_bytes > 0
            bytes_rx = bi2de(reshape( ...
                bits_so_far(1:n_bytes*8),8,[]).','left-msb');
            valid = min(n_bytes, numel(img_orig));
            img_rx(1:valid) = uint8(bytes_rx(1:valid));
        end
        axes(ax2);
        imshow(img_rx);
        title(sprintf('Received %d/%d packets', ...
            sum(pkt_received), num_packets));
        drawnow;
    else
        retransmit_cnt(p) = retransmit_cnt(p) + 1;
        if retransmit_cnt(p) >= params.max_retransmit
            fprintf('[!!] Max retransmissions — skipping packet %d\n', p);
            p = p + 1;
        end
    end

end

%%FINAL RESULTS
fprintf('\n════════════════════════════════\n');
fprintf('ARQ Complete\n');
fprintf('Packets received:     %d / %d\n', sum(pkt_received), num_packets);
fprintf('Total transmissions:  %d\n', total_tx);
fprintf('Retransmissions:      %d\n', total_tx - num_packets);
fprintf('Packet error rate:    %.4f\n', 1 - sum(pkt_received)/num_packets);
fprintf('Efficiency:           %.2f%%\n', num_packets/total_tx*100);

%% Save final reconstructed image
imwrite(img_rx, 'output_image.png');
fprintf('Saved output_image.png\n');

%% Final side-by-side comparison
figure;
subplot(1,2,1); imshow(img_orig); title('Original');
subplot(1,2,2); imshow(img_rx);   title('Reconstructed via ARQ');

%% Plot 1 — Constellation (most important)
figure;
scatter(real(all_eq_symbols), imag(all_eq_symbols), ...
    3, 'b', 'filled', 'MarkerFaceAlpha', 0.3);
hold on;
ideal = [1+1j  -1+1j  -1-1j  1-1j];
scatter(real(ideal), imag(ideal), 200, 'r', 'x', 'LineWidth', 3);
title('Received QPSK Constellation (All Packets)', 'FontSize', 14);
xlabel('In-phase', 'FontSize', 12);
ylabel('Quadrature', 'FontSize', 12);
legend('Received symbols', 'Ideal QPSK', 'FontSize', 11, 'Location', 'best');
grid on; axis equal;
axis([-3 3 -3 3]);

%% Plot 3 — Retransmissions Per Packet
figure;
bar(retransmit_cnt, 'FaceColor', [0.2 0.6 0.9]);
title('Retransmissions per Packet', 'FontSize', 14);
xlabel('Packet Index', 'FontSize', 12);
ylabel('Number of Retransmissions', 'FontSize', 12);
grid on;
yline(params.max_retransmit, 'r--', ...
    sprintf('Max = %d', params.max_retransmit), 'LineWidth', 1.5);
xlim([0 num_packets+1]);

%% Plot 4 — Cumulative Packets Received Over Time
figure;
cumulative = cumsum(pkt_received);
plot(1:num_packets, cumulative, 'b-', 'LineWidth', 2);
hold on;
plot(1:num_packets, 1:num_packets, 'r--', 'LineWidth', 1.5);
title('Cumulative Packets Received vs Transmitted', 'FontSize', 14);
xlabel('Packet Transmission Attempt', 'FontSize', 12);
ylabel('Packets Successfully Received', 'FontSize', 12);
legend('Actual received', 'Ideal (no errors)', ...
    'FontSize', 11, 'Location', 'northwest');
grid on;

%% Plot 5 — ARQ Efficiency
figure;
labels = {'Useful transmissions', 'Retransmissions'};
values = [num_packets  total_tx-num_packets];
colors = [0.2 0.7 0.3; 0.9 0.3 0.3];
pie_h  = pie(values);
pie_h(1).FaceColor = colors(1,:);
pie_h(3).FaceColor = colors(2,:);
title(sprintf('ARQ Efficiency: %.1f%%  (%d retx out of %d total)', ...
    num_packets/total_tx*100, total_tx-num_packets, total_tx), 'FontSize', 14);
legend(labels, 'FontSize', 11, 'Location', 'best');

%% Plot 6 — Final Image Comparison
figure('Name', 'Final Image Comparison', 'NumberTitle', 'off');
subplot(1,2,1);
imshow(img_orig);
title(sprintf('Original Image (%dx%d)', img_size(1), img_size(2)), ...
    'FontSize', 13);

subplot(1,2,2);
imshow(img_rx);
pkt_success_rate = sum(pkt_received)/num_packets*100;
title(sprintf('Reconstructed Image (%.1f%% packets OK)', pkt_success_rate), ...
    'FontSize', 13);

%% Print summary statistics for report
fprintf('Total packets:          %d\n', num_packets);
fprintf('Packets received OK:    %d (%.1f%%)\n', sum(pkt_received), pkt_success_rate);
fprintf('Total transmissions:    %d\n', total_tx);
fprintf('Retransmissions:        %d\n', total_tx-num_packets);
fprintf('ARQ efficiency:         %.1f%%\n', num_packets/total_tx*100);
fprintf('Max retransmit setting: %d\n', params.max_retransmit);
fprintf('Packet size:            %d bits\n', params.packet_size);
fprintf('Unique QPSK symbols:    %d\n', length(all_eq_symbols));