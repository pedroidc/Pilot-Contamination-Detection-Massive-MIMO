% This simulation intends to simulate a Massive MIMO system with multiple 
% users (Massive MIMOMU). The uplink and downlik stages are considered. In
% other words the Mobile Stations (MS), which are single antenna devices,
% send pilot symbols to the Base Station (BS), which have a large array of
% antennas. The BS, then, estimates the channel to perform beanforming and,
% finally, transmits information to the MS.
% It is also possible to simulate with perfect channel state information at
% the base station, for comparison purposes.
% The beamforming techniques used at the base station for the downlink
% stage are the Zero-Force and the Maximum-Ratio Processing
% 
clear variables; close all;
rng('shuffle')

%% System Parameters
Ntx = 64;   % Number of transmit antennas at the base station (BS)
Nrx = 1;    % Number of receive antennas for each user
Nusers = 4; % Number of users, i.e., mobile stations (MS)
Neve = 1;   % Number of eavesdropper
Np = 100;   % Number of pilot symbols
% fc = 28e9;                          % 28 GHz system
% vLight = physconst('LightSpeed');   % Speed of light
% waveLength = cLight/fc;             % Wave Length

%% Modulation Parameters
% The same modulation is used for all users.
modOrder = 4;       % Modulation order: 4: QPSK, 16: 16QAM, 64: 64QAM, 256: 256QAM
Nbpsym = log2(modOrder);   % Number of bits per symbol
Nc = 64;            % Number of OFDM subcarriers
Ncp = 16;           % Cyclic prefix length in samples

%% Define Modulators and Demodulators
hMod = comm.RectangularQAMModulator( ...
    'ModulationOrder', modOrder, ...
    'BitInput', true, ...
    'NormalizationMethod', 'Average power');
hDemod = comm.RectangularQAMDemodulator( ...
    'ModulationOrder', modOrder, ...
    'BitOutput', true, ...
    'NormalizationMethod', 'Average power');

%% Channel Parameters
% pp = exp(-(0:Nusers-1)/2);  % Channel power profile with exponential decay with factor 2
% pp = pp.'/sum(pp);          % Normalize to make the total power 1;
DF = eye(Ntx);

%% Simulation parameters
totalNumBits = 1e6; % Total number of bits transmited to each user
Nb = 1e3;           % Number of bits transmitted in each frame to each user
Nsymb = Nb/Nbpsym;
Ntrials = totalNumBits/ Nb; % Number of necessary trials
SNR = 0:2:20;
Nsnr = length(SNR);
ber = zeros(Nsnr, 1);
berAutUser = zeros(Nsnr, 1);
berEve = zeros(Nsnr, 1);
MSE = zeros(Nsnr, 1);

% Eveasdropper Power:
Pe = 1;
Pu = 1;

% Flags:
precoderStr = 'zero-forcing'; % Set this to choose between precoders:
                              % 'zero-forcing': Zero-forcing precoder
                              % 'maximum-ratio': Maximum-Ratio Precoder

chanEstStr = 'LS'; % Set this to choose between channele stimation techniques:
                   % 'Perfect-CSI': BS has perfect knowledge of the channel
                   % 'LS': Use Least Squares technique
                   % 'LMMSE': Use Linear Minimum Mean Squared Error
                   %          technique

pilotType = 'Random'; % Set this to choose pilots type
                      % 'Random': for random pilots
                      % 'Orthogonal': for orthogonal pilots
if(strcmp(chanEstStr, 'Perfect-CSI'))
    % If the perfect knowledge of CSI is considered at base station, there
    % us no need to transmit pilot symbols
    pilotType = '';
end

eveCon = 'n'; % Set this to choose between a scenario:
              % 'y': with pilot contamination
              % 'n': without pilot contamination
              
%% Simulation
for iSNR = 1:Nsnr
    
    fprintf('Simulation for SNR %.1f dB running\n', SNR(iSNR))
    
    for iTrial = 1:Ntrials
        
        % fprintf('Trial %.1f out of %.1f\n', iTrial, Ntrials)
        
        % Generating channels for trial ----------------------
        % Generate authentic users channels:
        Haut = sqrt(0.5)*DF*(randn(Ntx, Nusers) + 1j*randn(Ntx, Nusers));
        % Generate eavesdropper channel:
        g = sqrt(0.5)*(randn(Ntx, Neve) + 1j*randn(Ntx, Neve));
        
        % Uplink (MS to BS): ----------------------------------
        % MS pilot generation:
        if(strcmp(pilotType, 'Random'))
            pilotsBits = randsrc(Nusers*Np*Nbpsym, 1, [0 1]); % Random pilots
            xp = step(hMod, pilotsBits);
            xp = reshape(xp, Nusers, Np);
        elseif(strcmp(pilotType, 'Orthogonal'))
            % Not implemented yet
        end
        
        if(strcmp(eveCon, 'y'))
            xpe = sqrt(Pe)*xp(1, :);
            xptx = [xp; xpe];
            H = [Haut, g];
            xptx(1, :) = sqrt(Pu)*xptx(1, :);
%             H(:, 1) = [];
        else
            xptx = xp;
            H = Haut;
        end
        
        % Transmission MS to Bs
        yp = H*xptx;
        yp = awgn(yp, SNR(iSNR));
%         yp = awgn(yp, SNR(iSNR), 'measured');
        
        if(~strcmp(chanEstStr, 'Perfect-CSI'))
            % Channel estimation at the BS:
            if(strcmp(chanEstStr, 'LS'))
%                 Hest = (xp'*xp)\xp'*yp;
                Hest = ((conj(xp)*xp.')\conj(xp)*yp.').';
            elseif(strcmp(chanEstStr, 'LMMSE'))
                % Not implemented yet
            end
        else
            % Perfect CSI at the BS:
            Hest = Haut;
        end
        
        % Downlink (BS to MS): --------------------------------
        % Precoder:
        if(strcmp(precoderStr, 'zero-forcing'))
            WP = sqrt(Ntx - Nusers)*conj(Hest)/(Hest.'*conj(Hest)); % ZF Precoder
        else
            WP = conj(Hest)/sqrt(Ntx*trace(DF)); % Maximum-Ratio Precoder
        end
        
        % BS Data generation:
        bitStream = randsrc(Nusers*Nb, 1, [0 1]);   % Bit stream for every user
        modData = step(hMod, bitStream);            % Modulate
        modData = reshape(modData, Nusers, Nsymb);  % Reorganize symbols, 
                                                    % with each line 
                                                    % containing the symbols 
                                                    % to send to each user.
        x = WP*modData; % Perform precoding
        
        % BS to MS Channels:
        y = H.'*x;
%         y = awgn(y, SNR(iSNR));
        y = awgn(y, SNR(iSNR), 'measured');
        
        % MS Rx:
        demodData = step(hDemod, y(:));
        
        % General Bit error rate:
        autData = demodData(1:Nusers*Nb);
        [~, berTemp] = biterr(bitStream, autData);
        ber(iSNR) = ber(iSNR) + berTemp;
        
        % Bit Error Rate at the Authentic User
        bitsAutUser = bitStream(1:Nb);
        demodAutUser = demodData(1:Nb);
        [~, berTemp] = biterr(bitsAutUser, demodAutUser);
        berAutUser(iSNR) = berAutUser(iSNR) + berTemp;
        
        % Bit Error at the Eavesdropper
        if(strcmp(eveCon, 'y'))
            eveData = y(end, :);
            demodEve = step(hDemod, eveData(:));
            [~, berTemp] = biterr(bitsAutUser, demodEve);
            berEve(iSNR) = berEve(iSNR) + berTemp;
        end
        
        % Channel Mean Squared Error
        MSE(iSNR) = MSE(iSNR) + sum(abs(Haut(:) - Hest(:)).^2)/(Nusers*Ntx);
        
    end
    
end

ber = ber/Ntrials;
berAutUser = berAutUser/Ntrials;
berEve = berEve/Ntrials;
MSE = MSE/Ntrials;

figure
semilogy(SNR, MSE)
hold on;
grid on;
xlabel('SNR')
ylabel('MSE')

figure
semilogy(SNR, ber, 'DisplayName', 'Overal BER')
hold on;
semilogy(SNR, berAutUser, 'DisplayName', 'Authentic User')
semilogy(SNR, berEve, 'DisplayName', 'Eavesdropper User')
grid on;
legend;
xlabel('SNR (dB)')
ylabel('BER')

figure
bar(sqrt(abs(H.'*WP).^2))
ylabel('Gain')
xlabel('Received Node')
legend({'$x_1$', '$x_2$', '$x_3$', '$x_4$'}, ...
    'Interpreter', 'latex')