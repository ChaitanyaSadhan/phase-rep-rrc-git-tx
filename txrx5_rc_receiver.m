clear all;
clc
trystr = 25;
%git cmds.
%git pull
[status, cmdout] = system('git -C ./ pull'); disp(cmdout);
% pause(1)
%copying the recieved file.
cpcmd = ['cp /home/kcsadhan23/Desktop/Thesis/osfstorage-archive/Experiments/myExperiments/phaserep/phase-rep/audio_tx/call' num2str(trystr) 'tx.mat /home/kcsadhan23/Desktop/Thesis/osfstorage-archive/Experiments/myExperiments/phaserep/audio_tx/'];
[status, cmdout] = system(cpcmd);
pause(1)
%load params.mat
addpath audio_rx/
addpath audio_tx/
load("call"+num2str(trystr)+"rx.mat")
load("call"+num2str(trystr)+"tx.mat")
% Parameters
Fs = 8000;                % Sampling frequency
fc = 1200;                % Carrier frequency
M = 4;                    % QPSK (M = 4)
Ns = 20;                  % Samples per symbol
Rsym = Fs/Ns;             % Symbol rate
phase_offset = pi/15;     % Known phase offset
Nsym = 6;
N = 1000;

data = randi([0 M-1], N, 1);
mod_data = pskmod(data, M, pi/M, 'gray');
tx = upfirdn(mod_data, rcosdesign(0.35, Nsym, Ns, 'sqrt'), Ns, 1);

% Modulation
t = (0:length(tx)-1)'/Fs;
tx_carrier = real(tx .* exp(1j*2*pi*fc*t));

signal_in = tx_carrier./max(tx_carrier)/2/16;
actlen = size(signal_in,1);

%%inserting zeros.

numz = 1000;
zi = 8000;
empty = zeros(numz,1);
tx = [];
for i = 1:ceil(N*Ns/zi)
    start = (i-1)*zi+1;
    
    if i == ceil(N*Ns/zi)
        trim = signal_in(start:end);
        tx=[tx;trim];
    else
        trim = signal_in(start:start+zi-1);
        tx = [tx;trim;empty];
    end
end
duration = (length(tx))/8000

%{
%% pass through AMR channel
audiowrite("temp/txSignal.wav",signal_in,Fs);
[~,~] = system('ffmpeg.exe -y -i temp\txSignal.wav -ar 8000 -ab 2.4k temp\txSignal.amr');
[~,~] = system('ffmpeg.exe -y -i temp\txSignal.amr -ar 8000 -ab 2.4k temp\txSignal_out.wav');
signal_out1 = audioread('temp/txSignal_out.wav');
% Received from AMR channel
%}
dsignal_out = [signal_out1(41:end)*32; zeros(160,1)];
rx_carrier = signal_out(1:size(signal_in,1));

% Demodulation
rx_baseband = rx_carrier .* exp(-1j*2*pi*fc*t);
rx_matched = upfirdn(rx_baseband, rcosdesign(0.35, Nsym, Ns, 'sqrt'), 1, Ns);
rx_matched = rx_matched(Nsym+1:end-Nsym);
rx_corrected = rx_matched * exp(-1j*phase_offset);
rx_data = pskdemod(rx_corrected, M, pi/M, 'gray');

% scatter(real(rx_corrected), imag(rx_corrected), 'filled');
 
% BER Calculation
[num_errors, ber] = biterr(data, rx_data);
disp(['BER: ', num2str(ber)]);
disp(['bitrate: ', num2str(Rsym*log2(M))]);
