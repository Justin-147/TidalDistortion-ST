%计算超限率
function [zQDcxl,zSLcxl,timecx,timeqs,bs,umean,ustd]=ComCXL(ftdata,ttimej,fttimej,index2)
tmp=sort(ftdata,'descend');
tmp=tmp(length(tmp)/20+1:end);%去掉5%的坏点
umean=mean(tmp);%均值
ustd=std(tmp);%标准差
%bs=1:0.5:2.5;%阈值
bs=2;%阈值
tmp=reshape(ftdata,1440,length(ftdata)/1440);
zQDcxl=[];
zSLcxl=[];
timeqs=unique(floor(ttimej(index2)/10000));%存在缺数的天数
for kk=1:1:length(bs)
    yz=bs(kk)*ustd;%不同的阈值
    tmp1=abs(tmp-umean)-yz;
    tmp1(tmp1<=0)=0;
    QDcxl=sum(tmp1,1);%强度超限率
    SLcxl=sum(tmp1>0,1);%数量超限率
    zQDcxl=[zQDcxl;QDcxl];
    zSLcxl=[zSLcxl;SLcxl];
end
timecx=unique(floor(fttimej/10000));
end