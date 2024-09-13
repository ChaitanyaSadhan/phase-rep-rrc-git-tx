% code for testing transmission over a vocoder
clear all;
clc
%load params.mat
trystr = 25;
save_mat_str = "./audio_tx/call"+num2str(trystr)+"tx.mat";

N = 2000;
fcin = [9];
L = 30;
K = length(fcin);
rep = 1; 
Nc = N*rep;
tx_sym = zeros(K,N);
tx_symc = zeros(K,Nc);

%generating qpsk symbols: each value as 1 to 4
for k = 1: K
    tx_sym(k,:) = randi(4,1,N);
    %tx_sym(k,:) = ones(1,N);
    %tx_sym(:) = ceil((1:N)/N*4);%mod(1:N,4)+1;
    tx_symc(k,:) = repelem(tx_sym(k,:), rep);
end

%phase correction.
phase_correction = zeros(K,1);
for i = 1: K
    switch fcin(i)
        case {1,2}
            phase_correction(i) = pi/5.5;
        case {3,4}
            phase_correction(i) = pi/8.5;
        case {5,6,7,8}
            phase_correction(i) = pi/12;
        case {9,10,11,12}
            phase_correction(i) = pi/16;
        otherwise
            phase_correction(i) = 0;
    end
end
fc = (fcin+1)*100/8000; % normalized


% generate tx waveform
signal = zeros(Nc*L,1);
x = zeros(Nc,K);
mapping = [1+1j -1+1j 1-1j -1-1j];
for k = 1: length(fc)
    ctx = [cos(2*pi*fc(k)*[0:L-1]')+sin(2*pi*fc(k)*[0:L-1]') cos(2*pi*fc(k)*[0:L-1]'+pi)+sin(2*pi*fc(k)*[0:L-1]') cos(2*pi*fc(k)*[0:L-1]')+sin(2*pi*fc(k)*[0:L-1]'+pi) cos(2*pi*fc(k)*[0:L-1]'+pi)+sin(2*pi*fc(k)*[0:L-1]'+pi)]';
    signalm = ctx(tx_symc(k,:), :)';
    x(:,k) = mapping(tx_symc(k,:));
    signal = signal + signalm(:);
end
signal_in = signal./max(signal)/2/16;
actlen = size(signal_in,1);

%inserting zeros.
numz = 1000;
zi = 8000;
empty = zeros(numz,1);
tx = [];
for i = 1:ceil(N*L/zi)
    start = (i-1)*zi+1;
    
    if i == ceil(N*L/zi)
        trim = signal_in(start:end);
        tx=[tx;trim];
    else
        trim = signal_in(start:start+zi-1);
        tx = [tx;trim;empty];
    end
end
duration = (length(tx))/8000

%removing zeros.
% rx = [];
% for i = 1:6
%     start = (i-1)*11000+1;
%     rxtrim = tx(start:start+10000-1);
%     rx = [rx; rxtrim];
% end

%saving matfile.
save(save_mat_str,"signal_in", "tx", "tx_sym")
pause(1)
%pushing into git.
[status, cmdout] = system('git add -A');
disp(cmdout);
commitMessage = ['transmitting try ' num2str(trystr)];
[status, cmdout] = system(['git commit -m "' commitMessage '"']);
disp(cmdout);
[status, cmdout] = system('git push');
disp(cmdout);

%txing sound.
soundsc(tx,8000);


%% pass through AMR channel
audiowrite("temp/txSignal.wav",signal_in,8000);
[~,~] = system('ffmpeg.exe -y -i temp\txSignal.wav -ar 8000 -ab 4.75k temp\txSignal.amr');
[~,~] = system('ffmpeg.exe -y -i temp\txSignal.amr -ar 8000 -ab 4.75k temp\txSignal_out.wav');
signal_out1 = audioread('temp/txSignal_out.wav');

% pre-process received waveform
signal_out = [signal_out1(41:end)*32; zeros(160,1)];
signal_rec = signal_out(1:size(signal_in,1));
%signal_rec(1:41) = signal_rec(1:41);
% symbol demodulation

% extract phase
id_sym_rec = zeros(Nc,K);
y = zeros(Nc,K);
for i = 1:Nc
    rxsym = signal_rec(L*(i-1) + (1:L))*2;
    
    for k = 1:K
        I_mixed = rxsym .* cos(2*pi*fc(k)*[0:L-1]'+phase_correction(k));
        Q_mixed = rxsym .* sin(2*pi*fc(k)*[0:L-1]'+phase_correction(k));
        
        lim = round(floor(L*fc(k))/fc(k));

        I = mean(arrayfun(@(i) mean(I_mixed(i:i+lim-1)), 1:(length(I_mixed)-lim+1)));
        Q = mean(arrayfun(@(i) mean(Q_mixed(i:i+lim-1)), 1:(length(Q_mixed)-lim+1)));
        if (k==K)
            plot(I, Q, 'o'); hold on
        end
        y(i,k) = I + 1j*Q;
        id_sym_rec(i,k) = (1*(I<0)) + 2*(Q<0)+1;
    end
end


rx_sym = zeros(K,N);
for k = 1: length(fc)
    rd = reshape(id_sym_rec(:,k), rep, []);
    for i = 1:size(rd, 2)
        rx_sym(k,i) = mode(rd(:, i));
    end
end

% Plot x and y in color
hold on;
k = 1;
y1 = [];
for i = 1:Nc
    switch tx_symc(k,i)
        case 1
            color = 'red';
            y1 = [y(i,k);y1];
        case 2
            color = 'blue';
        case 3
            color = 'green';
        case 4
            color = 'black';
    end
    %plot(real(y(i,k)),imag(y(i,k)),'o',Color=color); 
end

SER = nnz(rx_sym(:)-tx_sym(:))/N;
%SERd = nnz(abs(x-xd)>0.1)/N;
bitrate = K*log2(size(ctx,1))*8000/size(ctx,2)/rep;
if (SER < 0.7)
    fprintf('Bitrate = %1.2f \t SER = %g\n', bitrate, SER)
    return
end

