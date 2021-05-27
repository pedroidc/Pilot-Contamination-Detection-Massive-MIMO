% This simulation intends to simulate the pilot contamination detection
% scheme of [1]. This system consists in a Multi User Massive Mimo system,
% with pilot contamination from neighboor cells and pilot contamination
% from an eavesdropper.
%
% [1] HASSAN, M. et al. Pilot contamination attack
%     detection for multi-cell MU-massive MIMO system.
%     AEU - International Journal of Electronics and
%     Communications, v. 113, p. 152945, 2020.
clear variables; close all; clc;
rng('shuffle')

%% System Parameters
% nAntennasRange = 20:4:100;    % Number of antennas at the base station (BS)
nAntennasRange = [64; 128; 140];
nAntennasTrials = length(nAntennasRange);
nPilotsRange = [12; 100; 128];
nPilotsTrials = length(nPilotsRange);
nUsers = 5;

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
Nb = 1e3;           % Number of bits transmitted in each frame to each user
nTrials =1e3; % Number of necessary trials

%% SNR parameters
SNR = -10:2:20;
nSNR = length(SNR);
N0 = 1./10.^(SNR/10);

%% Eavesdropper Parameters
Pe = 0:0.1:2;
nPe = length(Pe);

%% Initializations
Pd = zeros(nAntennasTrials, nPe, nPilotsTrials, nSNR);

%% Flags
pilotType = 'Random'; % Set this to choose pilots type
% 'Random': for random pilots
% 'Orthogonal': for orthogonal pilots

%% Simulation
for iPe = 1:nPe
    
    statusEve = sprintf('Simulation for eavesdropper power %.1f\n', Pe(iPe));
    fprintf(statusEve)
    
    for iAntennas = 1:nAntennasTrials
        
        % Cast variable for clear code:
        nAntennas = nAntennasRange(iAntennas);
        
        statusAntennas = sprintf('\tSimulation for %.0f Antennas\n', nAntennas);
        fprintf(statusAntennas)
        
        for iPilot = 1:nPilotsTrials
            
            % Cast variable for clear code:
            nPilots = nPilotsRange(iPilot);
            nBits = nUsers*nPilots*Nbpsym;
            
            statusPilots = sprintf('\t\tSimulation for %.0f pilot symbols\n', nPilots);
            fprintf(statusPilots)
            
            for iSNR = 1:nSNR
                
                statusSNR = sprintf('\t\t\t\tSimulation for SNR %.1f dB running\n', SNR(iSNR));
                fprintf(statusSNR)
                
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
                    xpe = sqrt(Pe(iPe))*xp(1, :);
                    
                    % Set final transmitted signal and channels
                    xptx = [xp; xpe];
                    H = [Haut, g];
                    
                    % Transmission MS to Bs
                    yp = H*xptx;
                    yp = awgn(yp, SNR(iSNR));
                    
                    % Least Squares Channel Estimation
                    Hest = yp*xp'/(xp*xp');
                    h1est = Hest(:, 1);
                    
                    % Pilot contamination attack detection
                    sovertau = nAntennas*N0(iSNR)/nPilots;
                    ln = log((2+sovertau)/(1+sovertau));
                    eta = (1 + sovertau)*(2+sovertau)*ln;
                    E = h1est'*h1est/nAntennas;
                    
                    Pd(iAntennas, iPe, iPilot, iSNR) = ...
                        Pd(iAntennas, iPe, iPilot, iSNR) + ...
                        (E >= eta);
                    
                end
                
                fprintf(repmat('\b', 1, numel(statusSNR)));
                
            end
            
            fprintf(repmat('\b', 1, numel(statusPilots)));
            
        end
        
        fprintf(repmat('\b', 1, numel(statusAntennas)));
        
    end
    
    fprintf(repmat('\b', 1, numel(statusEve)));
    
end

% Averaging:
Pd = Pd/nTrials;

%% Save results
fileName = 'results';
save(fileName, 'Pe', ...
    'nAntennasRange', 'nPilotsRange', 'SNR', 'Pd')

%% Figures
% figure
% plot(SNR, ratio)
% grid on;
% xlabel('SNR')
% ylabel('Ratio')
%
% figure
% plot(SNR, Pd)
% grid on;
% xlabel('SNR (dB)')
% ylabel('P_d')
%
% figure
% bar(sqrt(abs(H.'*WP).^2))
% ylabel('Gain')
% xlabel('Received Node')
% legend({'$x_1$', '$x_2$', '$x_3$', '$x_4$'}, ...
%     'Interpreter', 'latex')