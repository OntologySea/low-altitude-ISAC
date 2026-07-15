%% 5G NR PRACH Detection and False Alarm Test
% This example implements the physical random access channel (PRACH) missed
% detection and false alarm conformance tests, as defined in TS 38.141-1.
% You can measure the probability of correct detection of the PRACH
% preamble in the presence of a preamble signal or switch the PRACH
% transmission off to measure the false alarm probability.

% Copyright 2019-2022 The MathWorks, Inc.

%% Introduction
% The PRACH is an uplink transmission used by User Equipment (UE) to
% initiate synchronization with the gNodeB. TS 38.141-1 Section 8.4.1.5
% defines the probability of PRACH detection to be greater than or equal
% to 99% at specific SNR values for a set of PRACH configurations and
% propagation conditions. There are several detection error cases:
%
% * Detecting an incorrect preamble
% * Not detecting a preamble
% * Detecting the correct preamble but with the wrong timing estimation
%
% TS 38.141-1 states that a correct detection is achieved when the
% estimation error of the timing offset of the strongest path is less than
% the time error tolerance given in Table 8.4.1.1-1. For channel
% propagation conditions TDLC300-100 and PRACH preamble format 0, the time
% error tolerance is 2.55 microseconds.
%
% In this example, a PRACH waveform is configured and passed through an
% appropriate channel. At the receiver side, the example performs PRACH
% detection and calculates the PRACH detection probability. The example
% considers the parameters defined in TS 38.141-1 Table 8.4.1.5-1 and Table
% A.6-1. These are: normal mode (i.e., unrestricted set), 2 receive
% antennas, TDLC300-100 channel, normal cyclic prefix, burst format 0, SNR
% -6.0 dB. If you change the PRACH configuration to use one of the other
% PRACH preamble formats listed in Table A.6-1, you need to update the
% values of the time error tolerance and the SNR, according to TS 38.141-1
% Table 8.4.1.1-1 and Tables 8.4.1.5-1 to 8.4.1.5-3, respectively.

clc; clear; %close all;

%% Simulation Configuration
numPRACHSlots = 10;              % Number of PRACH slots to simulate at each SNR
%SNRdB = [-21, -18, -17, -16, -15, -12, -11, -10, -9, -8, -7,-6, -1]; % SNR range in dB
SNRdB = [-18, -15, -10, -6]; % SNR range in dB
foffset = 400.0;                 % Frequency offset in Hz
timeErrorTolerance = 2.55;       % Time error tolerance in microseconds
threshold = [];                  % Detection threshold
prachEnabled = true;             % Enable PRACH transmission. To simulate false alarm test, disable PRACH transmission.

%% ================= 基本参数 =================
c = 3e8;
distances = 100:1000:12000;
numIter = 20;

Nwin = 2;                 % 窗口数（核心参数）

%% Carrier Configuration
% Use the nrCarrierConfig configuration object |carrier| to specify the
% carrier settings. The example considers a carrier characterized by a
% subcarrier spacing of 120 kHz (FR2) and a bandwidth of 66 resource blocks.

carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = 120;
carrier.NSizeGrid = 66;

% Compute the OFDM-related information
ofdmInfo = nrOFDMInfo(carrier);
fs = ofdmInfo.SampleRate;

%% PRACH Configuration
% Table A.6-1 in TS 38.141-1 specifies the PRACH configurations to use for
% the PRACH detection conformance test.
%
% Set the PRACH configuration by using the nrPRACHConfig configuration
% object |prach|, according to Table A.6-1 and Section 8.4.1.4.2 in
% TS 38.141-1.

% Set PRACH configuration
prach = nrPRACHConfig;
prach.FrequencyRange = 'FR2';                    % Frequency range
prach.DuplexMode = 'TDD';                        % Time Division Duplexing
prach.ConfigurationIndex = 135;                   % Configuration index for format B4
prach.SubcarrierSpacing = 120;                   % Subcarrier spacing (must match carrier SCS for FR2 format B4)
prach.SequenceIndex = 0;                         % Logical sequence index
prach.PreambleIndex = 0;                         % Preamble index
prach.RestrictedSet = 'UnrestrictedSet';          % Normal mode
prach.FrequencyStart = 0;                         % Frequency location

% Define the value of ZeroCorrelationZone using the NCS table stored in
% the nrPRACHConfig object
switch prach.Format
    case {'0','1','2'}
        ncsTable = nrPRACHConfig.Tables.NCSFormat012;
        ncsTableCol = (string(ncsTable.Properties.VariableNames) == prach.RestrictedSet);
    case '3'
        ncsTable = nrPRACHConfig.Tables.NCSFormat3;
        ncsTableCol = (string(ncsTable.Properties.VariableNames) == prach.RestrictedSet);
    otherwise
        ncsTable = nrPRACHConfig.Tables.NCSFormatABC;
        ncsTableCol = contains(string(ncsTable.Properties.VariableNames), num2str(prach.LRA));
end
NCS = 0;
zeroCorrelationZone = ncsTable.ZeroCorrelationZone(ncsTable{:,ncsTableCol}==NCS);
prach.ZeroCorrelationZone = zeroCorrelationZone; % Cyclic shift index

%% Propagation Channel Configuration
% Use the nrTDLChannel object to configure the tapped delay line (TDL)
% propagation channel model |channel| as described in TS 38.141-1 Table
% 8.4.1.1-1.

channel = nrTDLChannel;
channel.DelayProfile = "TDLC300";           % Delay profile
channel.MaximumDopplerShift = 100.0;        % Maximum Doppler shift in Hz
channel.SampleRate = ofdmInfo.SampleRate;   % Input signal sample rate in Hz
channel.MIMOCorrelation = "Low";            % MIMO correlation
channel.TransmissionDirection = "Uplink";   % Uplink transmission
channel.NumReceiveAntennas = 2;             % Number of receive antennas
channel.NormalizePathGains = true;          % Normalize delay profile power
channel.Seed = 42;                          % Channel seed. Change this for different channel realizations
channel.NormalizeChannelOutputs = true;     % Normalize for receive antennas

% Get the channel characteristic information
channelInfo = info(channel);

%% Loop for SNR Values
% Use a loop to run the simulation for the set of SNR points given by the
% vector |SNRdB|.

% Initialize variables storing detection probability at each SNR
pDetection = zeros(size(SNRdB));

% Store the configuration parameters needed to generate the PRACH waveform
waveconfig.NumSubframes = prach.SubframesPerPRACHSlot;
waveconfig.Windowing = [];
waveconfig.Carriers = carrier;
waveconfig.PRACH.Config = prach;
waveconfig.PRACH.Enable = prachEnabled;

% The temporary variables 'prach_init', 'waveconfig_init', 'ofdmInfo_init',
% and 'channelInfo_init' are used to create the temporary variables
% 'prach', 'waveconfig', 'ofdmInfo', and 'channelInfo' within the SNR loop
% to create independent instances in case of parallel simulation
prach_init = prach;
waveconfig_init = waveconfig;
ofdmInfo_init = ofdmInfo;
channelInfo_init = channelInfo;

for snrIdx = 1:numel(SNRdB) % comment out for parallel computing
% parfor snrIdx = 1:numel(SNRdB) % uncomment for parallel computing
% To reduce the total simulation time, you can execute this loop in
% parallel by using the Parallel Computing Toolbox. Comment out the 'for'
% statement and uncomment the 'parfor' statement. If the Parallel Computing
% Toolbox(TM) is not installed, 'parfor' defaults to normal 'for' statement
    
    % Display progress in the command window
    timeNow = char(datetime('now','Format','HH:mm:ss'));
    fprintf([timeNow ': Simulating SNR = %+5.1f dB...'], SNRdB(snrIdx));

    % Set the random number generator settings to default values
    rng('default');
    
    % Initialize variables for this SNR point, required for initialization
    % of variables when using the Parallel Computing Toolbox
    prach = prach_init;
    waveconfig = waveconfig_init;
    ofdmInfo = ofdmInfo_init;
    channelInfo = channelInfo_init;
    
    % Reset the channel so that each SNR point will experience the same
    % channel realization
    reset(channel);
    
    % Normalize noise power to account for the sampling rate, which is a
    % function of the IFFT size used in OFDM modulation. The SNR is defined
    % per carrier resource element for each receive antenna.
    SNR = 10^(SNRdB(snrIdx)/10);
    N0 = 1/sqrt(2.0*channel.NumReceiveAntennas*double(ofdmInfo.Nfft)*SNR);
    
    % Detected preamble count
    detectedCount = 0;
    
    % Window step size (0.5 µs, tunable)
    win_shift = round(fs * 0.5e-6);

    % Loop for each PRACH slot
    numActivePRACHSlots = 0;
    for nSlot = 0:numPRACHSlots-1
        
        % Generate PRACH waveform for the current slot
        prach.NPRACHSlot = nSlot;
        waveconfig.PRACH.Config.NPRACHSlot = nSlot;
        [waveform,~,winfo] = hNRPRACHWaveformGenerator(waveconfig);
        
        % Set PRACH timing offset in microseconds as per TS 38.141-1 Figure
        % 8.4.1.4.2-2 and Figure 8.4.1.4.2-3
        if prach.LRA==139 % short preamble, values as in Figure 8.4.1.4.2-2
            baseOffset = ((winfo.WaveformResources.PRACH.Resources.PRACHSymbolsInfo.NumCyclicShifts/2)/prach.LRA)/prach.SubcarrierSpacing*1e3; % (microseconds)
            timingOffset = baseOffset + mod(nSlot,10)/10+1; % (microseconds)
        else % Long preamble, values as in Figure 8.4.1.4.2-3
            baseOffset = 0; % (microseconds)
            timingOffset = baseOffset + mod(nSlot,9)/10; % (microseconds)
        end
        sampleDelay = fix(timingOffset / 1e6 * ofdmInfo.SampleRate);
        
        % Generate transmit waveform
        txwave = [zeros(sampleDelay,1); waveform];

        % Pass data through channel model. Append zeros at the end of the
        % transmitted waveform to flush channel content. These zeros take
        % into account any delay introduced in the channel. This is a mix
        % of multipath delay and implementation delay. This value may
        % change depending on the sampling rate, delay profile and delay
        % spread
        rxwave = channel([txwave; zeros(channelInfo.MaximumChannelDelay, size(txwave,2))]);

        % Add noise
        noise = N0*complex(randn(size(rxwave)), randn(size(rxwave)));
        rxwave = rxwave + noise;

        % Skip this slot if the PRACH is inactive.
        % Skip the detection of this slot after advancing the channel to
        % make sure that the channel is always synchronized with the
        % current slot.
        % If the PRACH is inactive in this slot, the receiver should not
        % expect any PRACH transmission and thus should not even try to
        % detect a PRACH. Skipping the detection of an inactive slot is
        % particularly important when performing a conformance test. If the
        % PRACH is inactive, the reference waveform computed internally in
        % the |nrPRACHDetect| function is empty. This leads to an empty
        % correlation and thus to an empty detected preamble. This empty
        % preamble leads to an incorrect value of the detection
        % probability.
        if isempty(winfo.WaveformResources.PRACH.Resources.PRACHSymbols)
            continue;
        end
        numActivePRACHSlots = numActivePRACHSlots + 1;
        

        %% ================= Proposed（窗口移位） =================
        idetected = false;

        for n = 1:Nwin

            shift = (n-1)*win_shift;

            if shift + length(txwave) > length(rxwave)
                break;
            end

            segment = rxwave(shift+1 : shift+length(txwave));

            % Remove the implementation delay of the channel filter
            segment = segment((channelInfo.ChannelFilterDelay + 1):end, :);
        
            % Apply frequency offset
            t = ((0:size(segment, 1)-1)/channel.SampleRate).';
            segment = segment .* repmat(exp(1i*2*pi*foffset*t), 1, size(segment, 2));

            [ind2, ~] = nrPRACHDetect(carrier, prach, segment, 'DetectionThreshold', threshold);

            if ~isempty(ind2)
                idetected = true;
                break;
            end
        end

        % Remove the implementation delay of the channel filter
        rxwave = rxwave((channelInfo.ChannelFilterDelay + 1):end, :);
        
        % Apply frequency offset
        t = ((0:size(rxwave, 1)-1)/channel.SampleRate).';
        rxwave = rxwave .* repmat(exp(1i*2*pi*foffset*t), 1, size(rxwave, 2));

        % PRACH detection for all cell preamble indices
        %[detected, offsets] = nrPRACHDetect(carrier, prach, rxwave, 'DetectionThreshold', threshold);
        [detected, offsets] = nrPRACHDetect(carrier, prach, rxwave);
        
        % Test for preamble detection
        % if (length(detected)== 1)
        if (isscalar(detected))
            if ~prachEnabled
                % For the false alarm test, any preamble detected is wrong
                detectedCount = detectedCount + 1;
            else
                % Test for correct preamble detection
                if (detected==prach.PreambleIndex)

                    % Calculate timing estimation error
                    trueOffset = timingOffset/1e6; % (s)
                    measuredOffset = offsets(1)/channel.SampleRate;
                    timingerror = abs(measuredOffset-trueOffset);

                    % Test for acceptable timing error
                    if (timingerror<=timeErrorTolerance/1e6)
                        detectedCount = detectedCount + 1; % Detected preamble
                    end
                end
            end
        end
        
    end % of nSlot loop
    
    % Compute final detection probability for this SNR
    pDetection(snrIdx) = detectedCount/numActivePRACHSlots;

    % Display the detection probability for this SNR
    fprintf('Detection probability: %d%%\n', pDetection(snrIdx)*100);
    
end % of SNR loop

%% Results
% At the end of the SNR loop, the example plots the calculated detection
% probabilities for each SNR value against the target probability.

% Plot detection probability
figure('Name','Detection Probability');
plot(SNRdB,pDetection,'b-o','LineWidth',2,'MarkerSize',7);
title(['Detection Probability for ', num2str(numPRACHSlots) ' PRACH Slot(s)'] );
xlabel('SNR (dB)'); ylabel('Detection Probability');
grid on; hold on;
% Plot target probability
if prachEnabled
    % For a missed detection test, detection probability should be >= 99%
    pTarget = 99;
else
    % For a false alarm test, detection probability should be < 0.1%
    pTarget = 0.1; %#ok<UNRCH>
end
%plot(-6.0,pTarget/100,'rx','LineWidth',2,'MarkerSize',7);
plot(-10.0,pTarget/100,'rx','LineWidth',2,'MarkerSize',7);
legend('Simulation Result', ['Target ' num2str(pTarget) '% Probability'],'Location','best');
minP = 0;
if(~isnan(min(pDetection)))
    minP = min([pDetection(:); pTarget]);
end
axis([SNRdB(1)-0.1 SNRdB(end)+0.1 minP-0.05 1.05])

%% References
%
% # 3GPP TS 38.141-1. "NR; Base Station (BS) conformance testing. Part 1:
% Conducted conformance testing." _3rd Generation Partnership Project;
% Technical Specification Group Radio Access Network_.
% # 3GPP TS 38.104. "NR; Base Station (BS) radio transmission and
% reception." _3rd Generation Partnership Project; Technical Specification
% Group Radio Access Network_.
