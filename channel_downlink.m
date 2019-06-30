function [H_f] = channel_downlink(N_bs, N_ms, fc, Nc, Np, sigma_2_alpha, sigma, tau_max, fs, K)
%Inputs:
%   N_bs：基站天线数
%   N_ms：单天线用户
%   fc： 子载波频率
%   Nc： 簇
%   Np： 每一簇中路径数
%   sigma_2_alpha: average power of path, generally, sigma_2_alpha = 1;
%   sigma： 角度扩展
%   tau_max： 最大延时
%   fs： 采样率
%   K： 子载波数量
%Outputs:
%   H_f：信道矩阵 空频域 N_ms*N_bs*K

%% 设置参数
lambda = 3e8/fc;
d_ant = lambda/2;
Lp = Nc*Np;
angle = sigma*pi/180;
ang_min = 2*angle;
ang_max = pi-angle*2;

%% 生成导引矢量
n_MS = (0:(N_ms-1)).';	% transmitter for A_MS
n_BS = (0:(N_bs-1)).';	% receiver for A_BS
phi_BS = ang_min + (ang_max-ang_min)*rand(Nc,1); % 生成簇主径上的角度
phi = phi_BS*ones(1,Np) + ones(Nc,1)*((rand(1,Np)-0.5)*2*angle); % 生成簇上路径的角度
phi = reshape(phi.',[1, Lp]);

theta_MS = ang_min + (ang_max-ang_min)*rand(Nc,1);	
theta = theta_MS*ones(1,Np) +ones(Nc,1)* ((rand(1,Np)-0.5)*2*angle); % 生成簇上路径的角度
theta = reshape(theta.', [1, Lp]);
miu_BS = 2*pi*d_ant/lambda*cos(phi);
miu_MS = 2*pi*d_ant/lambda*cos(theta);
A_BS = exp(1i*n_BS*miu_BS)/sqrt(N_bs);
A_MS = exp(1i*n_MS*miu_MS)/sqrt(N_ms);

%% 生成时延 τ 路径增益 α
tau = tau_max*rand(1, Lp);
tau = sort(tau);
miu_tau = 2*pi*tau*fs/K;
alpha_temp = sqrt(sigma_2_alpha/2)*(randn(1,Lp) + 1i*randn(1, Lp));
alpha = sort(alpha_temp, 'descend');

%% 生成空频域信道
H_f = zeros(K, N_ms, N_bs);
for k = 1:K
    D_diag = sqrt(N_ms*N_bs/Lp)*diag(alpha.*exp(1i*(k-1)*miu_tau));
    H_f(k,:,:) = A_MS * D_diag * A_BS';
end
end