%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%载波间隔13kHz,1024路，共13.31M左右带宽，Ng=Nfft/4,用QPSK为21.3Mbps
%7.69*10^-5s=76.9us/OFDM符号   13k个OFDM符号/s 
%76.9/1027=0.075us每采样点即Ts=0.07512us
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555%%%%%%%%%%%%%
clear all
close all
clc
NgType=1; % NgType=1/2 for cyclic prefix/zero padding|对于CP或ZP，NgType=1或2
nt='CP';  
Ch=1;  % Ch=0/1 for AWGN/multipath channel|对于AWGN/多径瑞利信道，channelCh=0/1 
chType='CH'; Target_neb=500; 
%% 
PowerdB=[-4 -7.5 -9.5 -11 -15 -26 -30 -30 ]; % Channel tap power profile 'dB'|信道抽头功率特性 'dB'
Delay=[0 1 3 8 15 28 55 108];          % Channel delay 'sample'|信道时延（采样点）
%Delay=[0 1 2 3 4 5 6 7 ];
Power=10.^(PowerdB/10);     % Channel tap power profile 'linear scale'|信道抽头功率特性（线性尺度）
Ntap=length(PowerdB);       % Chanel tap number|信道抽头数
Lch=Delay(end)+1;           % Channel length|信道长度
%Lch=150;
%% 
Nbps=2; M=2^Nbps;  % Modulation order=2/4/6 for QPSK/16QAM/64QAM|调制阶数=2/4/6：QPSK/16QAM/64QAM
Nfft=1024;           % FFT size|FFT大小
Ng=Nfft/8;         % Ng=0: Guard interval length|GI的长度，没有保护间隔时，Ng=0
%Ng=3;  
Nsym=Nfft+Ng;      % Symbol duration|符号周期
norms=[1 sqrt(2) 0 sqrt(10) 0 sqrt(42)];
%% 
%%%%%%%导频设计
Nps=8;%导频间隔为4
Np=Nfft/Nps;%导频数目
Nd=Nfft-Np;%数据子信道数
%% 
SNR=[0:5:40];    % EbN0     % Number of iterations for each EbN0|对于每一个EbN0的迭代次数
Nframe=200;         % Number of symbols per frame|每一帧的符号数
file_name=['OFDM_BER_' chType '_' nt '_' 'GL' num2str(Ng) '.dat'];
fid=fopen(file_name, 'w+');
%% 
for i=1:length(SNR)
       randn('seed',1); rand('seed',1); 
       Ber=0;%误码率
       Neb=0; Ntb=0; % Initialize the number of error/total bits|初始化错误比特数/总比特数
       sigPow=0;         % Signal power initialization|初始信号功率
       MSE=zeros(1,9);
      for nsym=1:Nframe  %一帧有Nframe个符号
        X= randi([0,M-1],1,Nd); % bit: integer vector
        Xp = 2*(randn(1,Np)>0)-1; %产生导频
        X_mod= qammod(X,M)/norms(Nbps);
        x_GI=zeros(1,Nsym); %欲留带导频的发送序列
      %%   形成一个OFDM符号 
         ip=0;  pilot_loc=[];
         X_symbol=zeros(1,Nfft);%X_symbol为一个OFDM符号
         for n=1:Nfft %在特定位置插入导频
             if mod(n,Nps)==1   
                 X_symbol(n)=Xp(floor(n/Nps)+1) ;pilot_loc=[pilot_loc,n];ip=ip+1;
             else X_symbol(n)=X_mod(n-ip);
             end
         end
         x= ifft(X_symbol,Nfft); %点数为Nfft
         x_GI = [x(Nfft-Ng+1:Nfft) x];                  % Add CP|加循环前缀
         %% 信道仿真
            channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2); %8个抽头复信道
            h=zeros(1,Lch); h(Delay+1)=channel; % cir: channel impulse response|信道脉冲响应
           % h=h(1:Lch);%取前7径
            H_orignal= fft(h,Nfft); channel_length = length(h); % True channel and its time-domain length|实际信道和它的长度
           
            H_orignal=H_orignal/norm(H_orignal,'fro');
           h_temp=ifft(H_orignal,Nfft);
            h=h_temp(1:Lch);
            
            H_power_dB = 10*log10(abs(H_orignal.*conj(H_orignal)));     % True channel 
            snr = SNR(i)+10*log10(Nbps);%时域信噪比
         %  y = conv(x_GI,h_temp(1:Lch));   %这里用时域卷积 
           y = conv(x_GI,h);
      %end
          sigPow = mean(y.*conj(y));
          y_GI=awgn(y,SNR(i),'measured');
          noisePow=sigPow*10^(-SNR(i)/10);
          y_data_pilot=y_GI(Ng+1:Nsym);%去除保护间隔
          Y_data_pilot=fft(y_data_pilot);         %先fft,再频域均衡
     %% 不同估计方法比较
              for q=1:3
                  if q==1
                      [H_est,H_est2] = LS_CE( Y_data_pilot,Xp,pilot_loc,Nfft,Nps,'linear'); method='LS-linear'; % LS estimation with linear interpolation
                  elseif q==2
                      [H_est,H_est2] = LS_CE( Y_data_pilot,Xp,pilot_loc,Nfft,Nps,'spline'); method='LS-spline'; % LS estimation with spline interpolation
                  else
                      H_est3 = MMSE_CE( Y_data_pilot,Xp,pilot_loc,Nfft,Nps,h,10*log10(1/noisePow)); method='MMSE'; % MMSE estimation
                  end
                     htt=ifft(H_est2,Np);
                     hh=[htt zeros(1,Nfft-Np)];
                     H_LSdft=fft(hh,Nfft);
                     
                 H_est_power_dB = 10*log10(abs(H_est.*conj(H_est)));
                 h_est = ifft(H_LSdft); 
                 if(q==3)
                       h_est3=ifft(H_est3);
                      h_DFT = h_est3(1:Lch); 
                 else
                     h_DFT=h_est(1:Lch);
                 end
                 H_DFT = fft(h_DFT,Nfft); % DFT-based channel estimation
                 H_DFT_power_dB = 10*log10(abs(H_DFT.*conj(H_DFT)));
                 %H_DFT=H_LSdft;
%                  if (nsym==1&i==1)%只画第idB信噪比时
%                    figure(i+10), subplot(319+2*q), plot(H_power_dB,'b','linewidth',1); grid on; hold on;
%                    plot(H_est_power_dB,'r','Markersize',4,'linewidth',1);
%                    title(method); xlabel('Subcarrier Index'); ylabel('Power[dB]');
%                    legend('True Channel',method);  set(gca,'fontsize',9)
%                    title(['EbN0=',num2str(SNR(i)),'dB']);
%                    subplot(320+2*q), plot(H_power_dB,'b','linewidth',1); grid on; hold on;
%                    plot(H_DFT_power_dB,'r','Markersize',4,'linewidth',1);
%                    title([method ' with DFT']); xlabel('Subcarrier Index'); ylabel('Power[dB]');
%                    legend('True Channel',[method ' with DFT']);  set(gca,'fontsize',9)
%                     title(['SNR=',num2str(SNR(i)),'dB']);
%                  end
                 %MSE(q) = MSE(q) + (H_orignal-H_est)*(H_orignal-H_est)';
                 %MSE(q+3) = MSE(q+3) + (H_orignal-H_DFT)*(H_orignal-H_DFT)';
                    MSE(q) = MSE(q) + norm(H_orignal-H_est,'fro');
                    MSE(q+3) = MSE(q+3) + norm(H_orignal-H_LSdft,'fro');
                    MSE(q+6) = MSE(q+6) + norm(H_orignal-H_DFT,'fro');
              end
                H_shift=H_DFT; % Channel frequency response|信道频率响应
                Y_shift=Y_data_pilot;
                Xmod_r= Y_shift./H_shift;  % Equalizer - channel compensation|均衡器
             
      data=zeros(1,Nd);%存解调数据 
      ip = 0;
        for k=1:Nfft
         if mod(k,Nps)==1, ip=ip+1;  else  data(k-ip)=Xmod_r(k);  end
        end
       X_r=qamdemod(data*norms(Nbps),M);
      Neb=Neb+sum(sum(de2bi(X_r,Nbps)~=de2bi(X,Nbps)));
      Ntb=Ntb+Nfft*Nbps;  %[Ber,Neb,Ntb]=ber(bit_Rx,bit,Nbps); 
   end   %对应Nframe次个符号一帧  
   MSEs(i,:) = MSE/(Nframe);  
   
       if i==0
         sigPow= sigPow/Nsym/Nframe/N_iter;
         fprintf('Signal power= %11.3e\n', sigPow);
         fprintf(fid,'%%Signal power= %11.3e\n%%EbN0[dB]       BER\n', sigPow);
        else
         Ber = Neb/Ntb;     
         fprintf('SNR=%3d[dB], BER=%4d/%8d =%11.3e\n',SNR(i), Neb,Ntb,Ber)
         fprintf(fid, '%d\t%11.3e\n', SNR(i), Ber);
         if Ber<1e-6,  break;  end
       end
  end      
 %%%%%信道估计归一化均方误差
    for count=1:length(SNR)
        NMSE(count,:)=20*log10(   MSEs(count,:)./(  norm(H_orignal,'fro') ) );
    end  %计算归一化均方误差
       
%     figure(3)
%     plot(SNR,NMSE(:,1),'r--d');hold on
%     plot(SNR,NMSE(:,4),'b--s');hold on
%     legend('LS linear interpolate','LS linear interpolate whith DFT');
%     xlabel('SNR');ylabel('NMSE(dB)');title('LS linear interpolate for channel estimation');
%     figure(4)
%      plot(SNR,NMSE(:,2),'r--d');hold on
%      plot(SNR,NMSE(:,3),'b--h');hold on
%     plot(SNR,NMSE(:,5),'r--s');hold on
%     plot(SNR,NMSE(:,6),'b--+');
     figure(3)
    plot(SNR,NMSE(:,1),'r--d');hold on
    plot(SNR,NMSE(:,2),'b--d');hold on
    plot(SNR,NMSE(:,9),'k--d');hold on
   % plot(SNR,NMSE(:,3),'k--s');hold on
    plot(SNR,NMSE(:,4), 'g--d');hold on
   % plot(SNR,NMSE(:,7),'r--s');hold on
    legend('LS linear interpolate','LS spline interpolate ','MMSE with dft ','LS DFT interpolate');
    xlabel('SNR');ylabel('NMSE(dB)');title('LS linear interpolate for channel estimation');
 
if (fid~=0),  fclose(fid);   end
disp('Simulation is finished');
figure
plot_ber(file_name,Nbps);