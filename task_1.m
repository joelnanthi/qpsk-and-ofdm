%task_1.m
clear; 
clear; 
clc; 
close all;

%TX audio
disp('STARTING TRANSMITTER');
tx1_audio_ofdm;

load('tx_config.mat', 'tx_audio');

params = ofdm_params();
Fs = params.Fs;

%Play and record through audio channel

% Add silence before and after transmission
silence_time = 0.5; % seconds
silence = zeros(1, round(silence_time*Fs));

tx_audio_play = [silence, tx_audio, silence];

player = audioplayer(tx_audio_play, Fs);
recObj = audiorecorder(Fs, 16, 1);

disp('STARTING RECORDING');
record(recObj);

pause(0.2);

disp('PLAYING OFDM SIGNAL');
playblocking(player);

pause(0.4);

stop(recObj);
disp('STOPPED RECORDING');

rx_audio = getaudiodata(recObj).';

%Save recorded audio for receiver
audiowrite('recorded_ofdm.wav', rx_audio, Fs);

%Running receiver
disp('STARTING RECEIVER');
rx1_audio_ofdm;