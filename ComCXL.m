%���㳬����
function [zQDcxl,zSLcxl,timecx,timeqs,bs,umean,ustd]=ComCXL(ftdata,ttimej,fttimej,index2)
tmp=sort(ftdata,'descend');
tmp=tmp(length(tmp)/20+1:end);%ȥ��5%�Ļ���
umean=mean(tmp);%��ֵ
ustd=std(tmp);%��׼��
%bs=1:0.5:2.5;%��ֵ
bs=2;%��ֵ
tmp=reshape(ftdata,1440,length(ftdata)/1440);
zQDcxl=[];
zSLcxl=[];
timeqs=unique(floor(ttimej(index2)/10000));%����ȱ��������
for kk=1:1:length(bs)
    yz=bs(kk)*ustd;%��ͬ����ֵ
    tmp1=abs(tmp-umean)-yz;
    tmp1(tmp1<=0)=0;
    QDcxl=sum(tmp1,1);%ǿ�ȳ�����
    SLcxl=sum(tmp1>0,1);%����������
    zQDcxl=[zQDcxl;QDcxl];
    zSLcxl=[zSLcxl;SLcxl];
end
timecx=unique(floor(fttimej/10000));
end