%ofdm_params.m
function params = ofdm_params()
    
    params.Fs = 44100; %Sampling frequency
    params.Tsamp = 1/params.Fs; %0.0000227

    params.fc = 5000; %carrier frequency
    params.N_sc = 128; 
    params.N_cp = 20; 

    params.T_ofdm = 58e-3;
    params.F_ofdm = round(params.N_sc/params.T_ofdm); %2207hz
    params.T_cp = (params.N_cp/params.N_sc)*params.T_ofdm; %9.06ms
    params.T_sym = params.T_ofdm + params.T_cp; %67ms
    params.L = round(params.Fs/params.F_ofdm); %20 - Upsampling factor
    params.delta_f = 1/params.T_ofdm; %17.24hz

    %Header length
    params.header_length = 32;

    params.decoding_active = true;

    %Modulation
    params.M = 4; %QPSK
    params.bits_per_symbol = log2(params.M);

    %OFDM subcarrier structure
    params.N_active = 96; % data include continous pilot
    params.N_dc = 8; %DC carrier guard
    params.N_outer_zero =  params.N_sc-params.N_active-params.N_dc; %24 - Outer
    params.N_zero_each_side = params.N_outer_zero/2; %12 - Each side

    %Guard subcarrier index positions
    params.left_guard_idx = 1:12;
    params.dc_idx = 61:68;
    params.right_guard_idx = 117:128;

    %Continuous pilots
    params.pilot_idx = [25 45 85 105];
    params.pilot_value = [1 -1 1 -1];

    %Data subcarriers
    reserved_idx = [params.left_guard_idx, ...
                    params.dc_idx, ...
                    params.right_guard_idx, ...
                    params.pilot_idx];

    reserved_idx = sort(reserved_idx);

    params.data_idx = setdiff(1:params.N_sc, reserved_idx);
    params.data_idx = sort(params.data_idx);

    params.N_data = length(params.data_idx);

    %Tail guard
    params.N_tail_guard = 1;

    %Useful lengths
    params.ofdm_symbol_len = params.N_sc + params.N_cp;

    %Task 2 - Packet parameters
    params.packet_size       = 500;   % information bits per packet (≤ 1000)
    params.packet_id_length  = 16;    % bits for packet ID label
    params.crc_length        = 16;    % CRC bits appended per packet
    params.max_retransmit    = 10;    % max retransmissions per packet

    %ACK/NACK bit patterns - alternating pattern
    params.ACK_pattern  = repmat([1 0 1 0 1 0 1 0], 1, 25);  % 200 bits
    params.NACK_pattern = repmat([0 1 0 1 0 1 0 1], 1, 25);  % 200 bits

end