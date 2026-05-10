%task_1.m
clear; 
clear; 
clc; 
close all;

%TX audio
disp('Running transmitter...');
tx1_audio_ofdm;

load('tx_config.mat', 'tx_audio');

params = ofdm_params();
Fsamp = params.Fsamp;

%Play and record through audio channel

% Add silence before and after transmission
silence_time = 0.5; % seconds
silence = zeros(1, round(silence_time*Fsamp));

tx_audio_play = [silence, tx_audio, silence];

player = audioplayer(tx_audio_play, Fsamp);
recObj = audiorecorder(Fsamp, 16, 1);

disp('Start recording...');
record(recObj);

pause(0.2);

disp('Playing OFDM signal...');
playblocking(player);

pause(0.4);

stop(recObj);
disp('Recording stopped.');

rx_audio = getaudiodata(recObj).';

%Save recorded audio for receiver
audiowrite('recorded_ofdm.wav', rx_audio, Fsamp);

%Running receiver
disp('Running receiver...');
rx1_audio_ofdm;