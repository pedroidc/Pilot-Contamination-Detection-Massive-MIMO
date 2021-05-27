% This simulation intends to simulate the pilot contamination detection
% scheme of [1]. This simulation only analysis if the channel estimation
% is being done correctly when eavesdropper power is set to zero.
%
% [1] HASSAN, M. et al. Pilot contamination attack 
%     detection for multi-cell MU-massive MIMO system. 
%     AEU - International Journal of Electronics and 
%     Communications, v. 113, p. 152945, 2020. 
clear variables; close all;
rng('shuffle')

%% Fixed System Parameters
nAntennas = 200;    % Number of antennas at the base station (BS)
nPilots = 100; 		% Number of pilot symbols
nUsers = 4;         % Number of users
Pe = 0;             % Eavesdropper power
thr = 200;          % Threshold

%% Modulation
modOrder = 4;               % Modulation order: 4: QPSK, 16: 16QAM, 64: 64QAM, 256: 256QAM
Nbpsym = log2(modOrder);    % Number of bits per symbol

hMod = comm.RectangularQAMModulator( ...
    'ModulationOrder', modOrder, ...
    'BitInput', true, ...
    'NormalizationMethod', 'Average power');
hDemod = comm.RectangularQAMDemodulator( ...
    'ModulationOrder', modOrder, ...
    'BitOutput', true, ...
    'NormalizationMethod', 'Average power');

%% Simulation parameters
nBits = nPilots*nUsers*Nbpsym; % Number of bits transmitted in each frame to each user
nTrials = 1e3;  		% Number of necessary trials

%% SNR parameters
SNR = -10:2:20;
nSNR = length(SNR);
N0 = 1./10.^(SNR/10);

%% Initializations
MSE = zeros(nSNR, 1);

%% Flags
pilotType = 'Random'; % Set this to choose pilots type
% 'Random': for random pilots
% 'Orthogonal': for orthogonal pilots

%% Simulation
for iSNR = 1:nSNR
    
    fprintf('Simulation for SNR %.1f dB running\n', SNR(iSNR))
    
    for iTrial = 1:nTrials
        
        % fprintf('Trial %.1f out of %.1f\n', iTrial, Ntrials)
        
        % Generating channels for trial ----------------------
        % Generate authentic users channels:
        Haut = sqrt(0.5)*(randn(nAntennas, nUsers) + 1j*randn(nAntennas, nUsers));
        % Generate eavesdropper channel:
        g = sqrt(0.5)*(randn(nAntennas, 1) + 1j*randn(nAntennas, 1));
        
        % Uplink (MS to BS): ----------------------------------
        % MS pilot generation:
        if(strcmp(pilotType, 'Random'))
            
            pilotsBits = randsrc(nBits, 1, [0 1]); % Random pilots
            xp = step(hMod, pilotsBits);
            xp = reshape(xp, nUsers, nPilots);
            
        elseif(strcmp(pilotType, 'Orthogonal'))
            % Not implemented yet
        end
        
        % Set eavesdropper signal
        x1p = xp(1, :);
        xpe = sqrt(Pe)*x1p;
        
        % Set final transmitted signal and channels
        xptx = [xp; xpe];
        H = [Haut, g];
        
        % Transmission MS to Bs
        yp = H*xptx;
        yp = awgn(yp, SNR(iSNR));
        
        % Least Squares Channel Estimation
        Hest = yp*xp'/(xp*xp');
        h1est = Hest(:, 1);
        
        % MSE evaluation
        h1 = H(:, 1);
        MSE(iSNR) = MSE(iSNR) + sum(abs(h1 - h1est).^2)/nAntennas;
        
    end
    
end

% Averaging:
MSE = MSE/nTrials;

%% Figures
figure
semilogy(SNR, MSE)
grid on;
xlabel('SNR')
ylabel('MSE')

% figure
% plot(SNR, Pd)
% grid on;
% xlabel('SNR (dB)')
% ylabel('P_d')