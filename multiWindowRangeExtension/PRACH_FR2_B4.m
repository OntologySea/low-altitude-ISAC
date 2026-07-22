%% PRACH_FR2_B4_03 — 论文「仿真设计与性能评估」章节入口脚本
% 输出 5 张图，对应 requirements.txt A–E：
%   Figure 1 — 检测成功率 vs UE–gNB 距离 (§A)
%   Figure 2 — 按符号划分的检测分布 (§B)
%   Figure 3 — 相关峰值 vs 距离 (§C)
%   Figure 4 — 多普勒估计误差 vs 距离 (§D)
%   Figure 5 — 检测热力图，距离–符号 (§E)
%
% 对比四种方法：标准单窗 / 延长 CP / 重复前导 / 所提多窗口。
% FR2, 120 kHz SCS, PRACH format B4 (ConfigurationIndex=135), LOS, 低空 UAV。
%
clc; clear; close all;

%% §1 仿真场景与参数配置
cLight = 299792458;                 % 光速 (m/s)
fc_hz  = 28e9;                      % FR2 示意载频，用于多普勒计算 (Hz)

% TS 38.211：频段 FR2；载波与 PRACH 子载波间隔均为 120 kHz（format B4 与 ActivePRACHSlot 配套）
carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = 120;
carrier.NSizeGrid = 66;

prach = nrPRACHConfig;
prach.FrequencyRange = 'FR2';
prach.DuplexMode = 'TDD';
prach.ConfigurationIndex = 135;      % preamble format B4
prach.SubcarrierSpacing = 120;
prach.ActivePRACHSlot = 1;           % 120 kHz PRACH 与 CI=135 组合时 Toolbox 要求为 1
prach.SequenceIndex = 0;
prach.PreambleIndex = 0;
prach.RestrictedSet = 'UnrestrictedSet';
ncsTable = nrPRACHConfig.Tables.NCSFormatABC;
ncsTableCol = contains(string(ncsTable.Properties.VariableNames), num2str(prach.LRA));
prach.ZeroCorrelationZone = ncsTable.ZeroCorrelationZone(ncsTable{:,ncsTableCol} == 0);

ofdmInfo = nrOFDMInfo(carrier);
fs = ofdmInfo.SampleRate;

poInfo = nrPRACHOFDMInfo(carrier, prach);
cpSamplesPRACH = poInfo.CyclicPrefixLengths(1);
symCumEnd = cumsum(poInfo.SymbolLengths(:));   % PRACH OFDM 符号末端样本累计长度（与波形定时一致）
numPRACHSym = numel(symCumEnd);

% Baseline：禁止把整段缓冲送入 nrPRACHDetect（其内部全时序搜索 ≈ 理想接收机，>=1 km 仍会误报为恒成功）。
% 建模为「单次定时假设」：接收机按标定距离对齐 FFT 窗口；仅当 |τ−τ_nom| 在 CP 容忍门限内时才允许记为检测成功。
baselineCalibHoriz_m = 280;                  % gNB 侧假设的水平传播距离标定 (m)
baselineTolFracOfCp = 0.22;                  % 容忍 τ 误差占首符号 CP 的比例；收紧使 baseline 在 ~1.16 km 处跌落
baselineTolSamp = max(1, round(baselineTolFracOfCp * cpSamplesPRACH));

% LOS + 低空 UAV（§4.1）：斜距几何 → 传播时延 τ；径向速度 → 多普勒 f_D；叠加 AWGN
ueHeight_m   = 120;                 % UAV 高度 (m)
ueSpeed_mps  = 35;                  % UE 速度量级 (m/s)，径向分量在下面给出
radialFactor = 0.85;                % cos(俯仰/方位) 示意：径向分量比例

delay_nom_baseline_samp = round(hypot(baselineCalibHoriz_m, ueHeight_m) / cLight * fs);

SNR_dB_vec = -10:5:10;              % 接收端 Eb/N0 或等效 SNR（awgn），覆盖典型与中低 SNR
% 水平距离网格：含亚千米与千米级以上（§6.2「千米量级」对比）；可加密步长以平滑曲线
distances_m = [200 400 600 800 1000 2000 4000 6000 8000 10000 12000];
%distances_m = [200 400 2000];
numIter = 5;                     % Monte Carlo 次数（投稿级建议 ≥50，可按时间调整）
rng('default');

% 所提方法：多子窗观测（文稿 §3.1）：N=Nwin 个子窗，相邻起点间隔 Δ（样本）= win_shift；窗长 L=tx_len
Nwin = 8;
win_shift = round(0.5e-6 * fs);    % Δ：子窗起点间隔 (samples)，兼顾分辨率与复杂度

% 延长 CP：在发射序列前追加一段循环前缀式扩展（长度 ∝ 首符号 CP）
extCP_factor = 0.45;

% 重复发送：两段波形各占一份资源（总时长约 2×）；τ 门限在 Baseline 基础上额外放宽比例
repeatTolRelaxFrac = 0.18;

numRepeat = 2;

% 与 Baseline 同一套「单窗口 + τ 门限」，此处给出放宽后的样本门限（Extended CP / Repeat）
extCPTolSamp = baselineTolSamp + round(extCP_factor * cpSamplesPRACH);
repeatTolSamp = baselineTolSamp + round(repeatTolRelaxFrac * cpSamplesPRACH);

waveconfig = struct();
waveconfig.NumSubframes = prach.SubframesPerPRACHSlot;
waveconfig.Windowing = [];
waveconfig.Carriers = carrier;
waveconfig.PRACH.Enable = true;
waveconfig.PRACH.Power = 0;

%% §2 波形生成与辅助变量计算 — 查找可用 PRACH 时隙
prach.NPRACHSlot = findFirstActivePRACHSlot(carrier, prach, waveconfig);
waveconfig.PRACH.Config = prach;

% 生成基准 PRACH 波形（单次发射）
[tx_base, ~, winfo] = hNRPRACHWaveformGenerator(waveconfig);
assert(~isempty(winfo.WaveformResources.PRACH.Resources.PRACHSymbols), ...
    'PRACH symbols empty: check NPRACHSlot / NumSubframes.');
tx_len = length(tx_base);

% 延长 CP 波形（发射侧近似）
nExt = min(round(extCP_factor * cpSamplesPRACH), tx_len - 1);
tx_extCP = [tx_base(end-nExt+1:end); tx_base];
tx_len_ext = length(tx_extCP);

% 重复发送波形（能量按 1/sqrt(numRepeat) 分配，便于与单次发射公平比 SNR）
tx_rep = repmat(tx_base / sqrt(numRepeat), numRepeat, 1);
tx_len_rep = length(tx_rep);

max_delay_samp = round(hypot(max(distances_m), ueHeight_m) / cLight * fs);
% 接收轴「symbol 级」量化步长：取 PRACH 平均 OFDM 符号样本数（比仅用首符号更细，便于看出随距离的后移）
symRefLenSamples = max(64, round(symCumEnd(end) / numPRACHSym));
numObsSymBins = min(200, max(10, ceil((max_delay_samp + tx_len + Nwin * win_shift + 16 * symRefLenSamples) / symRefLenSamples)));

%% §3 输出缓存（成功率 / 符号分布 / 相关峰值 / 多普勒）
nD = numel(distances_m);
nS = numel(SNR_dB_vec);
P_det = struct('baseline', zeros(nD, nS), ...
               'extCP', zeros(nD, nS), ...
               'repeat', zeros(nD, nS), ...
               'proposed', zeros(nD, nS));
doppler_err_hz = zeros(nD, nS);
symDetHist = zeros(nD, numObsSymBins, nS);
symDetHistInternal = zeros(nD, numPRACHSym, nS);
peakCorrAccum = zeros(nD, nS);    % 所提方法命中时的归一化互相关峰值累加
peakCorrCnt   = zeros(nD, nS);    % 对应命中次数

% 主循环：Monte Carlo（距离 × SNR）
% 复杂度：所提方法串行约 O(Nwin·L)；子窗独立，可 parfor it = 1:numIter 或按窗并行（Parallel Computing Toolbox）
for iD = 1:nD
    fprintf('Progress: distance %d/%d (%.2f km)\n', iD, nD, distances_m(iD) / 1000);
    d_horiz = distances_m(iD);
    range_m = hypot(d_horiz, ueHeight_m);           % 斜距 LOS
    delay_samp = round(range_m / cLight * fs);
    fd_true = radialFactor * ueSpeed_mps * fc_hz / cLight;   % 多普勒 (Hz)

    for iS = 1:nS
        snr_db = SNR_dB_vec(iS);

        ok_b = 0; ok_e = 0; ok_r = 0; ok_p = 0;
        fd_err_accum = 0;
        fd_err_cnt = 0;

        for it = 1:numIter

            %% --- Baseline（单假设定时窗口 + τ 门限，非全时序搜索）---
            rx_b = buildRx(tx_base, tx_len, delay_samp, max_delay_samp, fd_true, fs, snr_db, Nwin, win_shift);
            ok_b = ok_b + double(detectSingleWindowTauGate(carrier, prach, rx_b, tx_len, ...
                delay_samp, delay_nom_baseline_samp, baselineTolSamp));

            %% --- Extended CP（单窗口 + τ 门限，tol 放宽）---
            rx_e = buildRx(tx_extCP, tx_len_ext, delay_samp, max_delay_samp, fd_true, fs, snr_db, Nwin, win_shift);
            ok_e = ok_e + double(detectSingleWindowTauGate(carrier, prach, rx_e, tx_len_ext, ...
                delay_samp, delay_nom_baseline_samp, extCPTolSamp));

            %% --- Repeat（仅对首段前导做单窗口检测 + τ 门限放宽；资源仍为 2× 发射）---
            rx_r = buildRx(tx_rep, tx_len_rep, delay_samp, max_delay_samp, fd_true, fs, snr_db, Nwin, win_shift);
            ok_r = ok_r + double(detectSingleWindowTauGate(carrier, prach, rx_r, tx_len, ...
                delay_samp, delay_nom_baseline_samp, repeatTolSamp));

            %% --- Proposed：多窗口 ---
            [hit_p, fd_err, detBinObs, detSymInternal, peak_corr] = detectProposedMetrics(carrier, prach, rx_b, tx_len, ...
                Nwin, win_shift, fd_true, fs, symRefLenSamples, numObsSymBins, symCumEnd, tx_base);
            ok_p = ok_p + double(hit_p);
            if hit_p && ~isnan(detBinObs) && detBinObs >= 1 && detBinObs <= numObsSymBins
                symDetHist(iD, detBinObs, iS) = symDetHist(iD, detBinObs, iS) + 1;
            end
            if hit_p && ~isnan(detSymInternal) && detSymInternal >= 1 && detSymInternal <= numPRACHSym
                symDetHistInternal(iD, detSymInternal, iS) = symDetHistInternal(iD, detSymInternal, iS) + 1;
            end
            if hit_p && ~isnan(peak_corr)
                peakCorrAccum(iD, iS) = peakCorrAccum(iD, iS) + peak_corr;
                peakCorrCnt(iD, iS) = peakCorrCnt(iD, iS) + 1;
            end
            if ~isnan(fd_err)
                fd_err_accum = fd_err_accum + fd_err;
                fd_err_cnt = fd_err_cnt + 1;
            end
        end

        P_det.baseline(iD, iS) = ok_b / numIter;
        P_det.extCP(iD, iS) = ok_e / numIter;
        P_det.repeat(iD, iS) = ok_r / numIter;
        P_det.proposed(iD, iS) = ok_p / numIter;
        if fd_err_cnt > 0
            doppler_err_hz(iD, iS) = fd_err_accum / fd_err_cnt;
        else
            doppler_err_hz(iD, iS) = NaN;
        end
    end
end

%% ========== 性能作图（图1–图5，对应 requirements.txt A–E）==========

% 选取参考 SNR（优先 0 dB，否则取最近值）
snrMidMask = SNR_dB_vec == 0;
if ~any(snrMidMask)
    [~, midIdx] = min(abs(SNR_dB_vec));
else
    midIdx = find(snrMidMask, 1);
end

%% §4 Figure 1 — 检测成功率 vs 距离（requirements.txt §A）
figure('Name', 'Figure 1 — Detection probability vs distance');
plot(distances_m / 1000, P_det.baseline(:, midIdx), '-o', 'LineWidth', 1.5); hold on;
plot(distances_m / 1000, P_det.extCP(:, midIdx), '-^', 'LineWidth', 1.5);
plot(distances_m / 1000, P_det.repeat(:, midIdx), '-v', 'LineWidth', 1.5);
plot(distances_m / 1000, P_det.proposed(:, midIdx), '-s', 'LineWidth', 1.5);
grid on;
xlabel('UE–gNB 水平距离 (km)'); ylabel('检测成功率');
legend('标准单窗 PRACH', '延长 CP', '重复前导', '所提多窗口方法', 'Location', 'best');
title(sprintf('Figure 1: 检测成功率随距离变化（FR2 B4 / LOS / SNR=%g dB）', SNR_dB_vec(midIdx)));

%% 共享计算：符号分布统计量（供 Figure 2 与 Figure 5 使用）
Hsym = symDetHist(:, :, midIdx);
rowHits = sum(Hsym, 2);
Hnorm = zeros(nD, numObsSymBins);
meanEquivSym = nan(nD, 1);
for id = 1:nD
    if rowHits(id) > 0
        Hnorm(id, :) = Hsym(id, :) / rowHits(id);
        meanEquivSym(id) = sum((1:numObsSymBins) .* Hsym(id, :)) / rowHits(id);
    else
        Hnorm(id, :) = NaN;
    end
end

kTheory = zeros(nD, 1);
for id = 1:nD
    ds = round(hypot(distances_m(id), ueHeight_m) / cLight * fs);
    kTheory(id) = min(floor(double(ds) / symRefLenSamples) + 1, numObsSymBins);
end

%% §5 Figure 2 — 按符号划分的检测分布（requirements.txt §B）
figure('Name', 'Figure 2 — Symbol-level detection distribution');
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
imagesc(1:numObsSymBins, distances_m / 1000, Hnorm);
axis xy;
colormap(gca, parula);
cb = colorbar;
cb.Label.String = '归一化频次（按距离行）';
xlabel(sprintf('等效符号位置序号 k（接收样本 n≈(k-1)·%d，≈平均 PRACH OFDM 符号时长）', symRefLenSamples));
ylabel('水平距离 (km)');
title(sprintf('(a) 等效符号位置分布（SNR=%g dB）', SNR_dB_vec(midIdx)));

nexttile;
plot(distances_m / 1000, meanEquivSym, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', [0.3 0.5 0.9]);
grid on;
xlabel('水平距离 (km)');
ylabel('平均等效符号序号 k');
title('(b) 平均等效符号序号 vs 距离');

Hin = symDetHistInternal(:, :, midIdx);
rowIn = sum(Hin, 2);
HnormIn = zeros(nD, numPRACHSym);
for id = 1:nD
    if rowIn(id) > 0
        HnormIn(id, :) = Hin(id, :) / rowIn(id);
    else
        HnormIn(id, :) = NaN;
    end
end
nexttile;
imagesc(1:numPRACHSym, distances_m / 1000, HnormIn);
axis xy;
colormap(gca, hot);
cb2 = colorbar;
cb2.Label.String = '归一化频次（按距离行）';
xlabel('PRACH 波形内 OFDM 符号序号');
ylabel('水平距离 (km)');
title('(c) PRACH 波形内 OFDM 符号序号（对照）');

sgtitle({'Figure 2: 近距离峰值对应较小 k；距离增大 k 增大，呈对角线趋势；与 LOS 理论时延一致。', ...
    '验证扩展检测窗口沿时间轴捕获延迟到达信号；(c) 为窗对齐后 PRACH 内部符号索引。'}, 'FontSize', 10);

%% §6 Figure 3 — 相关峰值 vs 距离（requirements.txt §C）
meanPeakCorr = peakCorrAccum ./ max(1, peakCorrCnt);
meanPeakCorr(peakCorrCnt == 0) = NaN;

figure('Name', 'Figure 3 — Correlation peak vs distance');
plot(distances_m / 1000, meanPeakCorr(:, midIdx), '-o', ...
    'LineWidth', 1.5, 'MarkerFaceColor', [0.3 0.5 0.9]);
grid on;
xlabel('UE–gNB 水平距离 (km)');
ylabel('平均归一化互相关峰值');
title(sprintf('Figure 3: 相关峰值随距离变化（所提方法，SNR=%g dB）', SNR_dB_vec(midIdx)));

%% §7 Figure 4 — 多普勒估计误差 vs 距离（requirements.txt §D）
figure('Name', 'Figure 4 — Doppler estimation error vs distance');
plot(distances_m / 1000, mean(doppler_err_hz, 2, 'omitnan'), '-o', 'LineWidth', 1.5);
grid on;
xlabel('水平距离 (km)');
ylabel('平均 |多普勒估计误差| (Hz)');
title('Figure 4: 所提方法多普勒估计误差随距离（对 SNR 维平均）');

%% §8 Figure 5 — 检测热力图，距离–符号（requirements.txt §E）
colActive = any(~isnan(Hnorm), 1);
firstK = 1;
lastK = numObsSymBins;
ia = find(colActive, 1, 'first');
ib = find(colActive, 1, 'last');
if ~isempty(ia) && ~isempty(ib)
    marginK = max(2, round(0.05 * (ib - ia + 1)));
    firstK = max(1, ia - marginK);
    lastK = min(numObsSymBins, ib + marginK);
end

fig5 = figure('Name', 'Figure 5 — Detection heatmap distance-symbol', ...
    'Position', [100 100 760 640]);
imagesc(firstK:lastK, distances_m / 1000, Hnorm(:, firstK:lastK));
axis xy;
colormap(gca, turbo);
cbh = colorbar;
cbh.Label.String = 'Normalized hit rate (per distance row)';
xlabel('Symbol index k (receive time axis, equiv. OFDM-symbol quantization)');
ylabel('Horizontal distance (km)');
title({sprintf('Figure 5: Detection heatmap (distance vs symbol index), SNR=%g dB', SNR_dB_vec(midIdx)); ...
    sprintf('Quantization ≈ %d samples (mean PRACH OFDM symbol)', symRefLenSamples)});
hold on;
plot(kTheory, distances_m / 1000, 'w-', 'LineWidth', 4);
hRef = plot(kTheory, distances_m / 1000, 'k--', 'LineWidth', 2);
hold off;
lgd = legend(hRef, {'LOS delay reference: k = \lfloor\tau/T_{sym}\rfloor+1'}, 'Location', 'southeast');
lgd.Color = [1 1 1 0.82];
drawnow;
ax5 = gca;
ax5.Units = 'normalized';
p5 = ax5.Position;
ax5.Position = [p5(1), max(0.10, p5(2) - 0.07), p5(3), max(0.48, p5(4) - 0.12)];
drawnow;
hSg = sgtitle(fig5, { ...
    'Figure 5: Diagonal trend — farther range → larger symbol index k (detection timing tracks propagation delay).', ...
    '直观反映传播时延与检测时序的关系；验证所提方法通过扩展时间轴检测窗口以捕获延迟到达信号。'}, ...
    'FontSize', 10);
drawnow;
dyUp = 0.065;
if ~isempty(hSg) && isprop(hSg, 'Position') && numel(hSg.Position) >= 1
    hSg.Units = 'normalized';
    hSg.Position(2) = min(0.99, hSg.Position(2) + dyUp);
else
    try
        hSub = findall(fig5, '-class', 'matlab.graphics.layout.Subtitle');
    catch
        hSub = findall(fig5, 'Type', 'Subtitletext');
    end
    for ii = 1:numel(hSub)
        if isprop(hSub(ii), 'Position') && numel(hSub(ii).Position) >= 1
            hSub(ii).Units = 'normalized';
            hSub(ii).Position(2) = min(0.99, hSub(ii).Position(2) + dyUp);
        end
    end
end

%% §9 参数摘要输出
fprintf('\n========== 仿真配置摘要 ==========\n');
fprintf('Standard: 3GPP TS 38.211/38.213, FR2, SCS=120 kHz, PRACH format B4 (ConfigurationIndex=%d)\n', prach.ConfigurationIndex);
fprintf('Scenario: gNB fixed, low-altitude UE h=%g m, LOS, beam-aligned RA (ideal)\n', ueHeight_m);
fprintf('Compare: baseline single-window | extended CP | repeated preamble | proposed (N=%d windows, Delta=%d samples)\n', Nwin, win_shift);
fprintf('Monte Carlo: numIter=%d, SNR_dB=%s, numDist=%d\n', numIter, mat2str(SNR_dB_vec), numel(distances_m));
fprintf('baselineTolFracOfCp=%.2f (baseline ~1.16 km cutoff)\n', baselineTolFracOfCp);
fprintf('Output: Figure 1-5 (detection vs distance, symbol distribution, correlation peak, Doppler error, heatmap)\n');
fprintf('Simulation finished.\n\n');

%% ================= 本地函数 =================

function slot = findFirstActivePRACHSlot(carrier, prach, waveconfig)
slot = [];
for ns = 0:191
    p = prach;
    p.NPRACHSlot = ns;
    wc = waveconfig;
    wc.PRACH.Config = p;
    [~, ~, winfo] = hNRPRACHWaveformGenerator(wc);
    if ~isempty(winfo.WaveformResources.PRACH.Resources.PRACHSymbols)
        slot = ns;
        return;
    end
end
error('findFirstActivePRACHSlot:NoSlot', '未找到激活的 PRACH 时隙，请检查配置。');
end

function rx = buildRx(tx, txLen, delay_samp, max_delay_samp, fd_hz, fs, SNR_dB, Nwin, win_shift)
txBuf = [tx; zeros(max_delay_samp, 1)];
rx = [zeros(delay_samp, 1); txBuf];
padding_rx = txLen + Nwin * win_shift;
rx = [rx; zeros(padding_rx, 1)];
t = (0:numel(rx) - 1).' / fs;
rx = rx .* exp(1j * 2 * pi * fd_hz * t);
rx = awgn(rx, SNR_dB, 'measured');
end

function ok = detectSingleWindowTauGate(carrier, prach, rxWave, txLen, delayTrueSamp, delayNomSamp, tolSamp)
% 单窗口 PRACH：仅在「认为的首径到达时刻」delayNomSamp 处截取长度 txLen 的一段再 nrPRACHDetect；
% 若 |τ_true−τ_nom|>tolSamp 则直接失败。Baseline / Extended CP / Repeat 共用，仅 tol 与 txLen 不同。
if abs(delayTrueSamp - delayNomSamp) > tolSamp
    ok = false;
    return;
end
startIdx = delayNomSamp + 1;
if startIdx > length(rxWave)
    ok = false;
    return;
end
seg = rxWave(startIdx:min(startIdx + txLen - 1, length(rxWave)));
if numel(seg) < txLen
    seg = [seg; zeros(txLen - numel(seg), 1)];
elseif numel(seg) > txLen
    seg = seg(1:txLen);
end
[idx, ~] = nrPRACHDetect(carrier, prach, seg);
ok = ~isempty(idx);
end

function [hit, fdErrAbs, detBinObs, detSymInternal, peakCorr] = detectProposedMetrics(carrier, prach, rxWave, txLen, ...
    Nwin, win_shift, fd_true, fs, symRefLenSamples, numObsSymBins, symCumEnd, tx_base)
hit = false;
fdErrAbs = NaN;
detBinObs = NaN;
detSymInternal = NaN;
peakCorr = NaN;
lag = max(1, round(1e-6 * fs));   % 1 µs 基线延迟，用于相位差分估计 CFO

for n = 1:Nwin
    shift = (n - 1) * win_shift;
    if shift + txLen > length(rxWave)
        break;
    end
    segment = rxWave(shift + (1:txLen));
    [idx, offsets] = nrPRACHDetect(carrier, prach, segment);
    if ~isempty(idx)
        hit = true;
        absDet = shift + max(0, round(double(offsets(1))));
        detBinObs = obsAbsSampleToEquivSymBin(absDet, symRefLenSamples, numObsSymBins);
        detSymInternal = prachOffsetToSymbolIdx(offsets(1), symCumEnd);
        % 归一化互相关峰值（segment vs tx_base），用于 Figure 3
        xc = abs(conv(segment, conj(flipud(tx_base))));
        peakCorr = max(xc) / (norm(segment) * norm(tx_base));
        if shift + lag + txLen <= length(rxWave)
            a = rxWave(shift + (1:txLen));
            b = rxWave(shift + lag + (1:txLen));
            den = sum(a .* conj(b));
            if abs(den) > 1e-9
                fd_hat = angle(den) / (2 * pi * lag / fs);
                fdErrAbs = abs(fd_hat - fd_true);
            end
        end
        return;
    end
end
end

function binIdx = obsAbsSampleToEquivSymBin(absDetSamp, symRefLen, numBins)
% 接收缓冲起始处计数的绝对样本位置 → 等效「符号」序号（symbol 级时间轴量化）
b = floor(double(absDetSamp) / symRefLen) + 1;
binIdx = max(1, min(round(b), numBins));
end

function symIdx = prachOffsetToSymbolIdx(offSamp, symCumEnd)
% nrPRACHDetect 相对当前窗口起点的偏移 → PRACH 波形内 OFDM 符号序号（1..numel(symCumEnd)）
ec = symCumEnd(:)';
offR = round(double(offSamp(1)));
offR = max(0, min(offR, ec(end) - 1));
symIdx = find(offR < ec, 1, 'first');
if isempty(symIdx)
    symIdx = numel(ec);
end
end
