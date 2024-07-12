%整体S变换
function [ST,ymax,maxi,ymean,Tim,TimQS,Cyc,Fre]=CommainS(ftdata,ttimej,fttimej,index2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fm=480;%计算的高频范围
chc=2880;%计算窗长（分钟），不要更改
buc=2520;%滑动步长（分钟），不要更改
tdc=(chc-buc)/2;%两侧剔除结果（分钟），不要更改
djs=360;%重采样窗长（分钟）
%参数说明，计算窗长2天的话，fm=480对应周期为2880/480=6分钟，这样既保证频率分量足够也能显著降低计算量，否则默认计算到2分种，计算存储量会变为3倍。
%每次循环计算重叠部分为chc-buc=360分钟，这是因为S变换汉宁窗衰减为边界5%的数据点即2880*5%=144分钟，这部分结果不可靠
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(ftdata);%数据长度(分钟)
m=floor((n-chc)/buc)+1;%迭代次数
cs=(chc-2*tdc)/djs;
ST=zeros(fm+1,m*cs);
disp(['需要迭代',num2str(m),'次']);
for i=1:1:m
    h=ftdata(1+(i-1)*buc:chc+(i-1)*buc);
    [st_h,~,Fre] = stm1(h,0,fm);
    tmp=st_h(:,tdc+1:chc-tdc);
    tmpa=reshape(tmp,fm+1,djs,cs);
    ST(:,1+(i-1)*cs:i*cs)=reshape(max(tmpa,[],2),fm+1,cs);
%     for jj=1:1:cs
%         ST(:,jj+(i-1)*cs)=max(tmp(:,1+(jj-1)*djs:djs*jj),[],2);%以重采样时窗内最大值作为重采样后各频率成分对应值
%     end
    %disp(['第',num2str(i),'次完成']);
end
[ymax,maxi]=max(abs(ST));%每时刻的最大能量及对应的频率坐标
ymean=mean(abs(ST));%每时刻的平均能量
Tim=unique(floor(fttimej/100));
Tim=Tim(tdc/60+1:djs/60:end);%
Tim=Tim(1:length(ymax));
TimQS=unique(floor(ttimej(index2)/100));%存在缺数的时间
Cyc=1./Fre;%信号周期，单位同采样率单位
disp('S变换完成');
%计算所读数据及其差分的超限率
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end