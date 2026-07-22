% =========================================================================
% IEEE ISAC PRACH: Multi-Window Receiver Simulation
% =========================================================================
clear; clc; close all;

%% --- 1. System & PRACH Parameters (系统与前导参数) ---
c = 3e8;                    % Speed of light (m/s)
fc = 28e9;                  % Carrier frequency FR2 (Hz)
lambda = c / fc;            % Wavelength (m)
SCS = 120e3;                % Subcarrier spacing (Hz)
T_sym = 1 / SCS;            % OFDM symbol duration (approx.)
L_zc = 139;                 % ZC sequence length (Short format)
F_samp = L_zc * SCS;        % Baseband sampling rate (Hz)
T_samp = 1 / F_samp;        % Sampling interval (s)

N_sym = 12;                 % Number of repeated symbols (e.g., Format B4)
N_windows = 8;              % Number of observation windows (Multi-window!)
N_obs = N_windows * L_zc;   % Total observation samples per symbol

% Generate ZC Root Sequence (Local Replica)
u = 1; % Root index
n = 0:L_zc-1;
zc_seq = exp(-1j * pi * u * n .* (n + 1) / L_zc).'; 

%% --- 2. Target Scenario (目标感知场景) ---
% 设定一个远超传统单窗口限制的远距离目标
true_dist = 4500;           % Target distance (4.5 km)
true_vel = 40;              % Target radial velocity (40 m/s)

% 物理量转换为基带参数 (两径雷达模型: 往返延迟)
delay_sec = 2 * true_dist / c; 
delay_samp = round(delay_sec / T_samp); 
fd = 2 * true_vel / lambda; % Doppler frequency shift
phase_shift_per_sym = 2 * pi * fd * T_sym;

%% ========================================================================
% Figure 1: 2D Range-Doppler Heatmap (单次快照，直观证明多窗口能力)
% ========================================================================
SNR_dB_demo = -5; % 演示热力图所用的 SNR

% 1. 构造接收信号张量 (Fast-time x Slow-time)
rx_tensor = zeros(N_obs, N_sym);
for m = 1:N_sym
    % 延迟信号 (使用零填充模拟)
    delayed_sig = zeros(N_obs, 1);
    if delay_samp + L_zc <= N_obs
        delayed_sig(delay_samp + 1 : delay_samp + L_zc) = zc_seq;
    end
    
    % 加入多普勒相位旋转 (Slow-time)
    doppler_phase = exp(1j * phase_shift_per_sym * (m - 1));
    
    % 加入高斯白噪声
    noise = (randn(N_obs, 1) + 1j * randn(N_obs, 1)) / sqrt(2);
    noise_power = 10^(-SNR_dB_demo/10);
    
    % 组合
    rx_tensor(:, m) = delayed_sig * doppler_phase + sqrt(noise_power) * noise;
end

% 2. 接收端处理：频域/时域相关 (Pulse Compression)
% 将接收信号与本地 ZC 序列进行滑动相关
corr_tensor = zeros(N_obs, N_sym);
for m = 1:N_sym
    % 为了简单高效，这里使用 filter 模拟相关滤波，然后截取正确的部分
    corr_out = filter(flipud(conj(zc_seq)), 1, rx_tensor(:, m));
    % 对齐延迟峰值
    corr_tensor(:, m) = circshift(corr_out, -L_zc + 1);
end

% 3. 多普勒处理 (Slow-time FFT)
N_fft_doppler = 256; % 增加 FFT 点数使频谱更平滑
RD_map = fftshift(fft(corr_tensor, N_fft_doppler, 2), 2);
RD_power = 10 * log10(abs(RD_map).^2);

% 计算坐标轴刻度
dist_axis = (0:N_obs-1) * T_samp * c / 2 / 1000; % 距离 (km)
vel_axis = linspace(-0.5, 0.5, N_fft_doppler) * (1/T_sym) * lambda / 2; % 速度 (m/s)

% --- Plot Fig 1 ---
figure('Name', 'Fig 1: Range-Doppler Map', 'Position', [100, 100, 600, 500]);
imagesc(vel_axis, dist_axis, RD_power);
colormap('jet'); colorbar;
set(gca, 'YDir', 'normal');
xlabel('Velocity (m/s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Distance (km)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Proposed Multi-Window ISAC PRACH (SNR = %d dB)', SNR_dB_demo), 'FontSize', 13);
% 标出真实目标位置
hold on;
plot(true_vel, true_dist/1000, 'ro', 'MarkerSize', 12, 'LineWidth', 2);
text(true_vel+5, true_dist/1000+0.5, 'Detected Target', 'Color', 'r', 'FontSize', 11, 'FontWeight', 'bold');

% 标出传统单窗口的截断红线
yline(L_zc * T_samp * c / 2 / 1000, 'w--', 'LineWidth', 2);
text(min(vel_axis)+5, L_zc * T_samp * c / 2 / 1000 + 0.3, 'Standard 1x Window Limit', 'Color', 'w', 'FontSize', 11, 'FontWeight', 'bold');

%% ========================================================================
% Figure 2: RMSE vs SNR (Monte Carlo 定量分析)
% ========================================================================
SNR_vec = -20:2:0;          % 扫描的信噪比范围
N_trials = 100;             % 蒙特卡洛仿真次数 (正式跑图建议设为 1000)

rmse_dist = zeros(size(SNR_vec));
rmse_vel = zeros(size(SNR_vec));

fprintf('Starting Monte Carlo Simulation for RMSE...\n');
for idx = 1:length(SNR_vec)
    snr = SNR_vec(idx);
    noise_power = 10^(-snr/10);
    
    err_d_sq = 0;
    err_v_sq = 0;
    
    for trial = 1:N_trials
        % 重复生成与 Fig 1 相似的接收张量
        rx = zeros(N_obs, N_sym);
        for m = 1:N_sym
            ds = zeros(N_obs, 1);
            if delay_samp + L_zc <= N_obs
                ds(delay_samp + 1 : delay_samp + L_zc) = zc_seq;
            end
            dp = exp(1j * phase_shift_per_sym * (m - 1));
            n_w = (randn(N_obs, 1) + 1j * randn(N_obs, 1)) / sqrt(2);
            rx(:, m) = ds * dp + sqrt(noise_power) * n_w;
        end
        
        % 相关处理
        c_tensor = zeros(N_obs, N_sym);
        for m = 1:N_sym
            tmp = filter(flipud(conj(zc_seq)), 1, rx(:, m));
            c_tensor(:, m) = circshift(tmp, -L_zc + 1);
        end
        
        % 参数联合估计 (2D Peak Search)
        RD = fft(c_tensor, N_fft_doppler, 2); % 不 shift 方便取下标
        [~, max_idx] = max(abs(RD(:)));
        [row, col] = ind2sub(size(RD), max_idx);
        
        % 提取距离与速度
        est_delay_samp = row - 1;
        est_dist = est_delay_samp * T_samp * c / 2;
        
        % 恢复多普勒频率 (考虑 FFT bin 的折叠)
        k = col - 1; 
        if k >= N_fft_doppler/2
            k = k - N_fft_doppler;
        end
        est_fd = k * (1 / (N_fft_doppler * T_sym));
        est_vel = est_fd * lambda / 2;
        
        % 累加误差
        err_d_sq = err_d_sq + (est_dist - true_dist)^2;
        err_v_sq = err_v_sq + (est_vel - true_vel)^2;
    end
    
    % 计算 RMSE
    rmse_dist(idx) = sqrt(err_d_sq / N_trials);
    rmse_vel(idx) = sqrt(err_v_sq / N_trials);
    
    fprintf('SNR = %3d dB | Dist RMSE = %6.2f m | Vel RMSE = %6.2f m/s\n', snr, rmse_dist(idx), rmse_vel(idx));
end

% --- Plot Fig 2 ---
figure('Name', 'Fig 2: RMSE vs SNR', 'Position', [750, 100, 800, 400]);

% 子图 1: 距离估计误差
subplot(1,2,1);
semilogy(SNR_vec, rmse_dist, '-o', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'b');
grid on;
xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Distance RMSE (m)', 'FontSize', 12, 'FontWeight', 'bold');
title('Ranging Performance', 'FontSize', 13);
set(gca, 'FontSize', 11);

% 子图 2: 速度估计误差
subplot(1,2,2);
semilogy(SNR_vec, rmse_vel, '-s', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'Color', 'r');
grid on;
xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Velocity RMSE (m/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Velocity Estimation Performance', 'FontSize', 13);
set(gca, 'FontSize', 11);