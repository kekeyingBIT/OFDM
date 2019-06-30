Ntap=4; 
PowerdB=[-4 -7.5 -9.5 -11 -15 -26 -30 -30 ];
Delay=[0 1 3 8 15 28 55 108];  
PowerdB=PowerdB(1:Ntap);
Delay=Delay(1:Ntap);
Power=10.^(PowerdB/10); 
Lch=Delay(end)+1; 
SNR=[0:5:30];
Nfft=1024;
channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2); 
h=zeros(1,Lch); h(Delay+1)=channel;
%H=fft(h,Nfft);
%H_orignal=H*norm(H,2);
%h=ifft(H_orignal,Nfft);
rolloff = 0.25;
span = 6;
sps =4;
M=4;
b=rcosdesign(rolloff,span,sps,'normal');
%h_up=upsample(h,4);
d= randi([0,M-1],1,100);
X_mod= qammod(d,M)/sqrt(2);

d_p=upfirdn(X_mod,b,sps);
d_up=conv(d_p,h);
%d_up=d_up*(mean(abs(d_p)))/(mean(abs(d_up)));
%y=awgn(d_up,100,'measured');
d_r=upfirdn(d_up,b,1,sps);

figure(1)
stem(abs(X_mod));hold on;

% X=d_r(span+1:span+length(X_mod));
% X=X/mean(abs(X));
stem(abs(d_p));hold on;
stem(abs(d_up)); hold on;
stem(abs(d_r)); hold on;
figure(2)
plot(real(d_p));
