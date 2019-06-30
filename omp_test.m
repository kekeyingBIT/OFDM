clc
clear all
close all
Nfft=1024;
Nps=[4 8 16 32 64];
Np_sel=Nfft./Nps;
%信道
%K=8稀疏
Ntap=8; 
PowerdB=[-4 -7.5 -9.5 -11 -15 -26 -30 -30 ];
Delay=[0 1 3 8 15 28 55 108];  
PowerdB=PowerdB(1:Ntap);
Delay=Delay(1:Ntap);

Power=10.^(PowerdB/10); 
Lch=Delay(end)+1; 

channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2); 
h=zeros(1,Nfft); h(Delay+1)=channel;
iter=100;
SNR=[0:5:40];
NMSE=zeros(length(Nps),length(SNR));
for rate=1:length(Nps)
    Np=Np_sel(rate);
    for term=1:length(SNR)
            MSE_omp=0;MSE_dft=0;
            for i=1:iter
                Xp=2*randi([0,1],1,Np)-1;
                %OMP
                Psi=dftmtx(Nfft);
                H=Psi*h';
                Phi=zeros(Np,Nfft);%
                sequence=randperm(Nfft);
                pilot_index=sort(sequence(1:Np));

                for ii=1:Np
                    Phi(ii,pilot_index(ii))=Xp(ii);
                end
                %Phi=randn(Np,Nfft);
                A=Phi*Psi;
                y=Phi*H;   %y为传输的导频信号
                y_noise=awgn(y,SNR(term),'measured');

                theta=CS_OMP(y_noise,A,Ntap);
                %OMP恢复
                MSE_omp=MSE_omp+norm(theta-h',2);

                 disp(['rate=',num2str(Nps(rate)),' ','SNR=',num2str(SNR(term)),'  ',num2str(i),'/',num2str(iter)])
            end
            MSE_omp=MSE_omp/(iter);

          NMSE(rate,term)=20*log10(MSE_omp./norm(h',2) ); 

    end
end
save NMSE NMSE
S = ['r--d';'y--d';'g--d';'b--d';'m--d'];
for j=1:length(Nps);
    plot(SNR,NMSE(j,:),S(j,:));hold on
end   
plot(SNR,-SNR,'k-s');hold on
legend('ratio=1/4','ratio=1/8','ratio=1/16','ratio=1/32','ratio=1/64','reference');
xlabel('SNR(dB)');ylabel('NMSE(dB)');
title('Nfft=1024');