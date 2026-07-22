import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import lfilter

# =========================================================================
# IEEE ISAC PRACH: Multi-Window Receiver Simulation (Python Version)
# =========================================================================

# --- 1. System & PRACH Parameters (系统与前导参数) ---
c = 3e8                     # Speed of light (m/s)
fc = 28e9                   # Carrier frequency FR2 (Hz)
lambda_c = c / fc           # Wavelength (m)
SCS = 120e3                 # Subcarrier spacing (Hz)
T_sym = 1 / SCS             # OFDM symbol duration (approx.)
L_zc = 139                  # ZC sequence length (Short format)
F_samp = L_zc * SCS         # Baseband sampling rate (Hz)
T_samp = 1 / F_samp         # Sampling interval (s)

N_sym = 12                  # Number of repeated symbols (e.g., Format B4)
N_windows = 8               # Number of observation windows (Multi-window!)
N_obs = N_windows * L_zc    # Total observation samples per symbol

# Generate ZC Root Sequence (Local Replica)
u = 1 # Root index
n = np.arange(L_zc)
# 注意：Python中虚数单位是 1j
zc_seq = np.exp(-1j * np.pi * u * n * (n + 1) / L_zc)

# --- 2. Target Scenario (目标感知场景) ---
true_dist = 4500            # Target distance (4.5 km)
true_vel = 40.0             # Target radial velocity (40 m/s)

# 物理量转换为基带参数 (两径雷达模型: 往返延迟)
delay_sec = 2 * true_dist / c
delay_samp = int(np.round(delay_sec / T_samp))
fd = 2 * true_vel / lambda_c # Doppler frequency shift
phase_shift_per_sym = 2 * np.pi * fd * T_sym

# ========================================================================
# Figure 1: 2D Range-Doppler Heatmap (单次快照，直观证明多窗口能力)
# ========================================================================
SNR_dB_demo = -5 # 演示热力图所用的 SNR
noise_power_demo = 10**(-SNR_dB_demo / 10)

# 1. 构造接收信号张量 (Fast-time x Slow-time)
rx_tensor = np.zeros((N_obs, N_sym), dtype=complex)
for m in range(N_sym):
    delayed_sig = np.zeros(N_obs, dtype=complex)
    if delay_samp + L_zc <= N_obs:
        delayed_sig[delay_samp : delay_samp + L_zc] = zc_seq
        
    doppler_phase = np.exp(1j * phase_shift_per_sym * m)
    noise = (np.random.randn(N_obs) + 1j * np.random.randn(N_obs)) / np.sqrt(2)
    rx_tensor[:, m] = delayed_sig * doppler_phase + np.sqrt(noise_power_demo) * noise

# 2. 接收端处理：频域/时域相关 (Pulse Compression)
corr_tensor = np.zeros((N_obs, N_sym), dtype=complex)
# 匹配滤波器的脉冲响应：翻转并取共轭
h = np.flipud(np.conjugate(zc_seq))

for m in range(N_sym):
    # 使用 lfilter 进行滑动相关
    corr_out = lfilter(h, [1.0], rx_tensor[:, m])
    # 将峰值对齐到真实延迟位置 (回退 L_zc - 1 个样本)
    corr_tensor[:, m] = np.roll(corr_out, -(L_zc - 1))

# 3. 多普勒处理 (Slow-time FFT)
N_fft_doppler = 256
# 沿慢时间维度 (axis=1) 做 FFT 并将零频移到中心
RD_map = np.fft.fftshift(np.fft.fft(corr_tensor, n=N_fft_doppler, axis=1), axes=1)
RD_power = 10 * np.log10(np.abs(RD_map)**2 + 1e-12) # 加极小值防止 log(0)

# 计算坐标轴刻度
dist_axis = np.arange(N_obs) * T_samp * c / 2 / 1000 # 距离 (km)
vel_axis = np.linspace(-0.5, 0.5, N_fft_doppler) * (1/T_sym) * lambda_c / 2 # 速度 (m/s)

# --- Plot Fig 1  (Grayscale/Black & White Version) ---
fig1 = plt.figure(figsize=(7, 5.5))
# 使用 pcolormesh 绘制热力图 (注意轴的对应关系)
X, Y = np.meshgrid(vel_axis, dist_axis)
plt.pcolormesh(X, Y, RD_power, cmap='gray_r', shading='auto')
plt.colorbar(label='Power (dB)')
plt.xlabel('Velocity (m/s)', fontsize=12, fontweight='bold')
plt.ylabel('Distance (km)', fontsize=12, fontweight='bold')
plt.title(f'Proposed Multi-Window ISAC PRACH (SNR = {SNR_dB_demo} dB)', fontsize=13)

# 标出真实目标位置：将原本的红圈和红字改为纯黑 ('black')
plt.plot(true_vel, true_dist/1000, marker='o', color='black', markersize=10, markerfacecolor='none', markeredgewidth=2)
plt.text(true_vel + 5, true_dist/1000 + 0.5, 'Detected Target', color='black', fontsize=11, fontweight='bold')

# 标出传统单窗口的截断红线：将原本的白色虚线和白字改为纯黑 ('black')
standard_limit_km = L_zc * T_samp * c / 2 / 1000
plt.axhline(y=standard_limit_km, color='black', linestyle='--', linewidth=2)
plt.text(np.min(vel_axis) + 5, standard_limit_km + 0.3, 'Standard 1x Window Limit', color='black', fontsize=11, fontweight='bold')
plt.tight_layout()

# 紧跟在 Fig 1 画完之后，针对 fig1 句柄进行保存
fig1.savefig('fig1_heatmap.pdf', format='pdf', bbox_inches='tight')

# ========================================================================
# Figure 2: RMSE vs SNR (Monte Carlo 定量分析)
# ========================================================================
SNR_vec = np.arange(-20, 2, 2)
N_trials = 1000 # 蒙特卡洛仿真次数 (正式跑图建议设为 500-1000)

rmse_dist = np.zeros(len(SNR_vec))
rmse_vel = np.zeros(len(SNR_vec))

print('Starting Monte Carlo Simulation for RMSE...')
for idx, snr in enumerate(SNR_vec):
    noise_power = 10**(-snr / 10)
    err_d_sq = 0.0
    err_v_sq = 0.0
    
    for trial in range(N_trials):
        # 1. 生成接收信号
        rx = np.zeros((N_obs, N_sym), dtype=complex)
        for m in range(N_sym):
            ds = np.zeros(N_obs, dtype=complex)
            if delay_samp + L_zc <= N_obs:
                ds[delay_samp : delay_samp + L_zc] = zc_seq
            dp = np.exp(1j * phase_shift_per_sym * m)
            n_w = (np.random.randn(N_obs) + 1j * np.random.randn(N_obs)) / np.sqrt(2)
            rx[:, m] = ds * dp + np.sqrt(noise_power) * n_w
            
        # 2. 相关处理
        c_tensor = np.zeros((N_obs, N_sym), dtype=complex)
        for m in range(N_sym):
            tmp = lfilter(h, [1.0], rx[:, m])
            c_tensor[:, m] = np.roll(tmp, -(L_zc - 1))
            
        # 3. 参数联合估计 (2D Peak Search)
        # 注意：寻找峰值时不需要 fftshift，直接按索引找，找完再解算频率
        RD = np.fft.fft(c_tensor, n=N_fft_doppler, axis=1)
        row, col = np.unravel_index(np.argmax(np.abs(RD)), RD.shape)
        
        # 提取距离
        est_delay_samp = row
        est_dist = est_delay_samp * T_samp * c / 2
        
        # 恢复多普勒频率 (处理 FFT 折叠问题)
        k = col
        if k >= N_fft_doppler / 2:
            k = k - N_fft_doppler
        est_fd = k * (1 / (N_fft_doppler * T_sym))
        est_vel = est_fd * lambda_c / 2
        
        # 累加误差平方
        err_d_sq += (est_dist - true_dist)**2
        err_v_sq += (est_vel - true_vel)**2
        
    # 计算 RMSE
    rmse_dist[idx] = np.sqrt(err_d_sq / N_trials)
    rmse_vel[idx] = np.sqrt(err_v_sq / N_trials)
    print(f"SNR = {snr:3d} dB | Dist RMSE = {rmse_dist[idx]:6.2f} m | Vel RMSE = {rmse_vel[idx]:6.2f} m/s")

# --- Plot Fig 2 ---
fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5)) # 这里把原来的 fig 改成了 fig2

# 子图 1: 距离估计误差
ax1.semilogy(SNR_vec, rmse_dist, '-ob', linewidth=2, markersize=7)
ax1.grid(True, which="both", ls="--")
ax1.set_xlabel('SNR (dB)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Distance RMSE (m)', fontsize=12, fontweight='bold')
ax1.set_title('(a) Ranging Performance', fontsize=13)

# 子图 2: 速度估计误差
ax2.semilogy(SNR_vec, rmse_vel, '-sr', linewidth=2, markersize=7)
ax2.grid(True, which="both", ls="--")
ax2.set_xlabel('SNR (dB)', fontsize=12, fontweight='bold')
ax2.set_ylabel('Velocity RMSE (m/s)', fontsize=12, fontweight='bold')
ax2.set_title('(b) Velocity Estimation Performance', fontsize=13)

plt.tight_layout()

# 针对 fig2 句柄进行保存
fig2.savefig('fig2_rmse.pdf', format='pdf', bbox_inches='tight')

# ========================================================================
# Figure 4: RMSE of Range vs Number of Windows (N_w)
# ========================================================================
Nw_vec = np.arange(1, 9) # 1 to 8
SNR_fixed = -5 # Fixed SNR for ablation study
noise_power_fixed = 10**(-SNR_fixed / 10)
N_trials_nw = 500

rmse_dist_nw = np.zeros(len(Nw_vec))
print('\nStarting Monte Carlo Simulation for RMSE vs Nw...')

for idx, nw in enumerate(Nw_vec):
    n_obs_nw = nw * L_zc
    err_d_sq_nw = 0.0
    
    for trial in range(N_trials_nw):
        rx_nw = np.zeros((n_obs_nw, N_sym), dtype=complex)
        for m in range(N_sym):
            ds_nw = np.zeros(n_obs_nw, dtype=complex)
            if delay_samp < n_obs_nw:
                end_idx = min(delay_samp + L_zc, n_obs_nw)
                ds_nw[delay_samp : end_idx] = zc_seq[0 : end_idx - delay_samp]
                
            dp = np.exp(1j * phase_shift_per_sym * m)
            n_w = (np.random.randn(n_obs_nw) + 1j * np.random.randn(n_obs_nw)) / np.sqrt(2)
            rx_nw[:, m] = ds_nw * dp + np.sqrt(noise_power_fixed) * n_w
            
        c_tensor_nw = np.zeros((n_obs_nw, N_sym), dtype=complex)
        for m in range(N_sym):
            tmp = lfilter(h, [1.0], rx_nw[:, m])
            c_tensor_nw[:, m] = np.roll(tmp, -(L_zc - 1))
            
        RD_nw = np.fft.fft(c_tensor_nw, n=N_fft_doppler, axis=1)
        row, col = np.unravel_index(np.argmax(np.abs(RD_nw)), RD_nw.shape)
        
        est_dist_nw = row * T_samp * c / 2
        err_d_sq_nw += (est_dist_nw - true_dist)**2
        
    rmse_dist_nw[idx] = np.sqrt(err_d_sq_nw / N_trials_nw)
    print(f"Nw = {nw} | Dist RMSE = {rmse_dist_nw[idx]:6.2f} m")

# --- Plot Fig 4 ---
fig4 = plt.figure(figsize=(6, 4.5))
plt.plot(Nw_vec, rmse_dist_nw, '-og', linewidth=2, markersize=8)
plt.grid(True, which="both", ls="--")
plt.xlabel('Number of Windows ($N_w$)', fontsize=12, fontweight='bold')
plt.ylabel('Distance RMSE (m)', fontsize=12, fontweight='bold')
plt.title(f'Ranging Performance vs. $N_w$ (SNR = {SNR_fixed} dB)', fontsize=13)

min_nw = int(np.ceil((delay_samp + L_zc) / L_zc))
plt.axvline(x=min_nw, color='red', linestyle='--', linewidth=2)
plt.text(min_nw + 0.2, np.max(rmse_dist_nw)*0.5, f'Min $N_w$ required', color='red', fontsize=11, fontweight='bold')
plt.tight_layout()

fig4.savefig('fig4_rmse_nw.pdf', format='pdf', bbox_inches='tight')

# 最后再让两张图在屏幕上显示出来 (可选)
# plt.show()
