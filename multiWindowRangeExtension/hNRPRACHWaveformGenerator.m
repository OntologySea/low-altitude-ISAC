%hNRPRACHWaveformGenerator Generate a 5G NR PRACH waveform  
%   [WAVEFORM,GRIDSET,INFO] = hNRPRACHWaveformGenerator(WAVECONFIG) 
%   generates a 5G NR physical random access channel (PRACH) waveform
%   WAVEFORM given input WAVECONFIG parameters. The function also returns
%   two structure arrays, GRIDSET and INFO.
%   GRIDSET is a structure array containing the following fields for each
%   carrier in WAVECONFIG.Carriers: 
%   ResourceGrid            - PRACH resource grid 
%   Info                    - Structure with information corresponding to
%                             the PRACH OFDM modulation. If the PRACH is
%                             configured for FR2 or the PRACH slot for the
%                             current configuration spans more than one
%                             subframe, some of the OFDM-related
%                             information may be different between PRACH
%                             slots. In this case, the info structure is an
%                             array of the same length as the number of
%                             PRACH slots in the waveform.
%   INFO is a structure containing the following field: 
%   WaveformResources - Structure containing the following field:
%       PRACH - Structure containing the following field:
%           Resources - Structure array containing the following fields:
%               NPRACHSlot       - PRACH slot numbers of each allocated
%                                  PRACH preamble
%               PRACHSymbols     - PRACH symbols corresponding to each
%                                  allocated PRACH slot
%               PRACHSymbolsInfo - Additional information associated with
%                                  PRACH symbols
%               PRACHIndices     - PRACH indices corresponding to each
%                                  allocated PRACH slot
%               PRACHIndicesInfo - Additional information associated with
%                                  PRACH indices
% 
%   WAVECONFIG is a structure containing the following fields:
%   NumSubframes            - Number of 1ms subframes in generated waveform
%   DisplayGrids            - Display the resource grid after signal
%                             generation (0,1) (optional - default: 0)
%   Windowing               - Number of time-domain samples over which
%                             windowing and overlapping of OFDM symbols is
%                             applied. If absent or set to [], the default
%                             value is used, as detailed in <a
%                             href="matlab: doc('nrPRACHOFDMModulate')"
%                             >nrPRACHOFDMModulate</a>. (optional)
%   Carriers - Carrier-specific configuration object, as described in
%              <a href="matlab:help('nrCarrierConfig')">nrCarrierConfig</a> with these properties:
%       SubcarrierSpacing   - Subcarrier spacing in kHz
%       CyclicPrefix        - Cyclic prefix
%       NSizeGrid           - Number of resource blocks
%   PRACH    - Structure containing the following fields:
%       Enable              - Enable or disable this PRACH configuration
%                             (0,1) (optional - default: 1)
%       AllocatedPreambles  - Index (0-based) of the allocated PRACH
%                             preambles to transmit. This field considers
%                             only the active PRACH preambles. Set this
%                             value to 'all' to include all the active
%                             PRACH preambles in the waveform.
%       Power               - PRACH power scaling in dB. This parameter
%                             represents beta_PRACH (in dB) in TS 38.211
%                             Section 6.3.3.2.
%       Config              - PRACH-specific configuration object, as described in
%                             <a href="matlab:help('nrPRACHConfig')">nrPRACHConfig</a> with these properties:
%           FrequencyRange      - Frequency range (used in combination with
%                                 DuplexMode to select a configuration table
%                                 from TS 38.211 Table 6.3.3.2-2 to 6.3.3.2-4)
%           DuplexMode          - Duplex mode (used in combination with
%                                 FrequencyRange to select a configuration table
%                                 from TS 38.211 Table 6.3.3.2-2 to 6.3.3.2-4)
%           ConfigurationIndex  - Configuration index, as defined in TS 38.211
%                                 Tables 6.3.3.2-2 to 6.3.3.2-4
%           SubcarrierSpacing   - PRACH subcarrier spacing in kHz
%           SequenceIndex       - Logical root sequence index
%           PreambleIndex       - Scalar preamble index within cell
%           RestrictedSet       - Type of restricted set
%           ZeroCorrelationZone - Cyclic shift configuration index
%           RBOffset            - Starting resource block (RB) index of the
%                                 initial uplink bandwidth part (BWP) relative
%                                 to carrier resource grid
%           FrequencyStart      - Frequency offset of lowest PRACH transmission
%                                 occasion in frequency domain with respect to
%                                 PRB 0 of the initial uplink BWP
%           FrequencyIndex      - Index of the PRACH transmission occasions in
%                                 frequency domain
%           TimeIndex           - Index of the PRACH transmission occasions in
%                                 time domain
%           ActivePRACHSlot     - Active PRACH slot number within a subframe or
%                                 a 60 kHz slot
%           NPRACHSlot          - PRACH slot number
%
%   Example 1:
%   % Generate a 10ms PRACH waveform for the default values for
%   % nrPRACHConfig and nrCarrierConfig. Display the PRACH-related OFDM
%   % information.
%
%   waveconfig.NumSubframes = 1;
%   waveconfig.Carriers = nrCarrierConfig;
%   waveconfig.PRACH.Config = nrPRACHConfig;
%   [waveform,gridset] = hNRPRACHWaveformGenerator(waveconfig);
%   disp(gridset.Info)
%
%   Example 2:
%   % Generate a 10ms PRACH waveform for a single carrier for the default
%   % values for nrPRACHConfig and nrCarrierConfig, considering only the
%   % first and third active PRACH preambles. Display the PRACH-related
%   % OFDM information and check which preambles are transmitted.
%
%   waveconfig.NumSubframes = 10;
%   waveconfig.Carriers = nrCarrierConfig;
%   waveconfig.PRACH.Config = nrPRACHConfig;
%   waveconfig.PRACH.AllocatedPreambles = [0 2];
%   [waveform,gridset,info] = hNRPRACHWaveformGenerator(waveconfig);
%   disp([info.WaveformResources.PRACH.Resources.NPRACHSlot])
%
%   See also nrPRACHOFDMModulate, nrPRACHOFDMInfo, nrPRACHConfig,
%   nrPRACHGrid, nrPRACH, nrPRACHIndices, nrCarrierConfig.

%   Copyright 2019-2022 The MathWorks, Inc.

function [waveform,gridset,winfo] = hNRPRACHWaveformGenerator(waveconfig)
    
    % Unbundle carrier and PRACH specific parameter structure for easier access
    carrier = waveconfig.Carriers;
    prach = waveconfig.PRACH.Config;
    numSubframes = waveconfig.NumSubframes;
    
    % Defaulting for the grid plotting
    if ~isfield(waveconfig,'DisplayGrids')
        waveconfig.DisplayGrids = 0;
    end
    
    % Windowing
    if isfield(waveconfig,'Windowing')
        windowing = waveconfig.Windowing;
    else
        windowing = [];
    end
    
    % Defaulting for PRACH enable field
    if isfield(waveconfig.PRACH,'Enable')
        prachEnabled = waveconfig.PRACH.Enable;
    else
        prachEnabled = 1;
    end
    
    % Allocated PRACH preambles to transmit
    if isfield(waveconfig.PRACH,'AllocatedPreambles')
        allocatedPreambles = waveconfig.PRACH.AllocatedPreambles;
    else
        allocatedPreambles = 'all';
    end
    
    % PRACH power scaling
    if isfield(waveconfig.PRACH,'Power')
        prachPower = waveconfig.PRACH.Power;
    else
        prachPower = 0;
    end

    % Get the OFDM information
    prachOFDMInfo = nrPRACHOFDMInfo(carrier, prach, 'Windowing', windowing);
    
    % Total number of PRACH slots in the generated waveform.
    % Note that a PRACH preamble can last more than one subframe. In this
    % case, there may be a non-integer number of PRACH slots within the
    % chosen number of subframes. For this reason, the maximum number of
    % full PRACH preambles within the time resources given by
    % |numSubframes| is considered. The generated waveform is then padded
    % with zeros to cover the remaining subframes up to |numSubframes|.
    numPRACHSlots = floor(numSubframes / prach.SubframesPerPRACHSlot);
    if mod(numSubframes, prach.SubframesPerPRACHSlot)
        warnTxt = "The chosen number of subframes (" + string(numSubframes) + ...
            ") is not an integer multiple of the number of subframes spanned by a nominal PRACH slot for the chosen PRACH configuration (" + ...
            string(prach.SubframesPerPRACHSlot) + "). The generated waveform contains PRACH for the first " + string(numPRACHSlots) + " PRACH slots only.";
        warning('nr5g:hNRPRACHWaveformGenerator:NumSubframesLessThanPRACHSlot',warnTxt);
    end
    
    % Get the array of allocated active PRACH preambles to transmit
    if (strcmpi(allocatedPreambles,'all'))
        allocatedPreambles = 0:numPRACHSlots;
    end
    
    % Initialize outputs and internal resources
    waveform = [];
    gridset = struct('ResourceGrid',[],'Info',prachOFDMInfo);
    resourceGrid = [];
    resources = repmat(struct('NPRACHSlot',[],'PRACHSymbols',[],'PRACHSymbolsInfo',[],...
        'PRACHIndices',[],'PRACHIndicesInfo',[]),numPRACHSlots,1);

    % Initialize the active PRACH preamble counter (0-based)
    preambleCount = 0;

    % Set the number of antennas
    P = 1;
    
    % Loop over all PRACH slots
    firstSlot = prach.NPRACHSlot;
    for nSlot = firstSlot + (0:numPRACHSlots-1)
        
        % Update PRACH slot number
        prach.NPRACHSlot = nSlot;
        
        % Generate an empty PRACH resource grid
        prachGrid = nrPRACHGrid(carrier, prach, P);
        
        % Create the PRACH symbols
        [prachSymbols, prachInfoSym] = nrPRACH(carrier, prach);
        
        % If this preamble is to be transmitted, generate the waveform
        if ~isempty(prachSymbols) && any(allocatedPreambles==preambleCount)
            % Create the PRACH indices and retrieve PRACH information
            [prachIndices, prachInfoInd] = nrPRACHIndices(carrier, prach);
            
            % Map the PRACH symbols into the grid
            prachGrid(prachIndices) = prachSymbols * db2mag(prachPower) * prachEnabled;
            
            % Capture resource info for this PRACH instance
            resources(nSlot-firstSlot+1).PRACHSymbols = prachSymbols * db2mag(prachPower);
            resources(nSlot-firstSlot+1).PRACHIndices = prachIndices;
            resources(nSlot-firstSlot+1).PRACHIndicesInfo = prachInfoInd;
        end
        resources(nSlot-firstSlot+1).NPRACHSlot = prach.NPRACHSlot;
        resources(nSlot-firstSlot+1).PRACHSymbolsInfo = prachInfoSym;
        
        % Generate the PRACH waveform for this slot and append it to the
        % existing waveform
        if prachEnabled
            [wave, prachOFDMInfo] = nrPRACHOFDMModulate(carrier, prach, prachGrid, 'Windowing', windowing);
        else
            % If the PRACH transmission is disabled, manually create a
            % waveform of zeros and only compute the PRACH OFDM info.
            prachOFDMInfo = nrPRACHOFDMInfo(carrier, prach, 'Windowing', windowing);
            N = sum(prachOFDMInfo.SymbolLengths) + prachOFDMInfo.OffsetLength;
            wave = complex(zeros(N, P));
        end
        waveform = [waveform; wave]; %#ok<AGROW>
        
        % Capture the OFDM modulation info
        if prach.SubframesPerPRACHSlot > 1 || strcmpi(prach.FrequencyRange,'FR2')
            % For long preambles that span more than 1 subframe, some of
            % the fields in the OFDM information (e.g.,
            % CyclicPrefixLengths, GuardLengths, and SymbolLengths) can
            % vary for each PRACH preamble
            % For FR2, the length of the cyclic prefix can be different
            % between PRACH slots, depending on how many times the PRACH
            % occasion crosses time instants 0 and 0.5 ms.
            gridset.Info(nSlot-firstSlot+1) = prachOFDMInfo;
        else
            % For all the preambles that span at most 1 subframe, the OFDM
            % information is the same in each PRACH slot
            gridset.Info = prachOFDMInfo;
        end
        
        % Update the resource grid
        resourceGrid = [resourceGrid prachGrid]; %#ok<AGROW>
        
        % Update PRACH preamble count if PRACH preamble is active in
        % this slot
        preambleCount = preambleCount + ~isempty(prachSymbols);
    end
    
    % Capture all resources info for this PRACH configuration
    waveinfo.PRACH.Resources = resources;
    
    % Zero-padding of the generated waveform to guarantee the correct
    % length given by |numSubframes|
    lengthPadding = (prachOFDMInfo.SampleRate/1e3)*(numSubframes) - length(waveform);
    waveform = [waveform; complex(zeros(lengthPadding, P))];

    % Capture the resource grid
    gridset.ResourceGrid = resourceGrid;

    % Plot the resource grid
    if waveconfig.DisplayGrids
        plotResourceGrid(gridset, prachPower);
    end
    
    winfo.WaveformResources = waveinfo;
end

% Plot the resource grid
function plotResourceGrid(gridset, prachPower)

    % Display the resource grid plot in a new figure
    cmap = parula(64);
    figure('Name','Resource Grid');
    axRG = axes;
    image(axRG,40*abs(gridset.ResourceGrid)/db2mag(prachPower));
    axis(axRG,'xy');
    colormap(axRG,cmap);
    title(axRG,sprintf('PRACH Resource Grid (Size [%s])',strjoin(string(size(gridset.ResourceGrid)),' ')));
    xlabel(axRG,'Symbols'); ylabel(axRG,'Subcarriers');
end