Ntap=6; 
PowerdB=[-4 -7.5 -9.5 -11 -15 -26 -30 -30 ];
Delay=[0 1 3 8 15 28 55 108];  
PowerdB=PowerdB(1:Ntap);
Delay=Delay(1:Ntap);
Power=10.^(PowerdB/10); 
Lch=Delay(end)+1; 
SNR=[0:5:30];
Nfft=1024;
nmse=zeros(1,length(SNR));

    for iter=1:1
            channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2); 
            h=zeros(1,Lch); h(Delay+1)=channel;

                rolloff = 0.25;
                span = 6;
                sps =4;
                M=4;
                b=rcosdesign(rolloff,span,sps,'sqrt');

%                 d= randi([0,M-1],1,100);
%                 X_mod= qammod(d,M)/sqrt(2);

                d_x=upsample(h,4);
                d_p=[0  0 d_x(1:length(d_x)-2)];
                figure(10)
                stem(abs(d_x));hold on;
                stem(abs(d_p));
                title('改变抽头位置到非整数采样点');
                d_up=conv(d_p,b);   %脉冲成型后
                %加噪声
                d_up=awgn(d_up,100,'measured');
                d_re=conv(d_up,b);
                d_d=downsample(d_re,4);

            delay=sps;
            d_down=zeros(delay,length(d_re));
            d_recover=zeros(1,length(d_re));
            for ii=1:delay
                len=length(downsample(d_re(ii:length(d_re)),4));
                d_down(ii,1:len)=downsample(d_re(ii:length(d_re)),4);
            end
            count=1;
            for jj=1:delay:length(d_re)
                count=count+1;
                d_recover(jj:jj+delay-1)=d_down(:,count);
            end
        if(iter==1)
        figure(1)
        stem(abs(h));hold on;
        temp=d_d(span+1:span+length(h));
        stem(abs(temp));
        %nmse=norm(temp-abs(h),2)./norm(abs(h),2)

        title('原信道与下采样信道')
          figure(2)
         stem(abs(d_re)); hold on; %升余弦后接收
         title('包括成型、匹配的真实信道');
         figure(3)
         stem(abs(d_recover));
         title('估计的下采样信道重构真实信道');
         figure(4)
         difference=abs(d_re)-abs(d_recover);
         stem(difference);
         title('完全正确估计下的重构误差');
        end
     nmse(term)=nmse(term)+(norm(difference,2)./norm(abs(d_re),2));
    end
    nmse(term)=20*log10(nmse(term));

