%����S�任
function [ST,ymax,maxi,ymean,Tim,TimQS,Cyc,Fre]=CommainS(ftdata,ttimej,fttimej,index2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fm=480;%����ĸ�Ƶ��Χ
chc=2880;%���㴰�������ӣ�����Ҫ����
buc=2520;%�������������ӣ�����Ҫ����
tdc=(chc-buc)/2;%�����޳���������ӣ�����Ҫ����
djs=360;%�ز������������ӣ�
%����˵�������㴰��2��Ļ���fm=480��Ӧ����Ϊ2880/480=6���ӣ������ȱ�֤Ƶ�ʷ����㹻Ҳ���������ͼ�����������Ĭ�ϼ��㵽2���֣�����洢�����Ϊ3����
%ÿ��ѭ�������ص�����Ϊchc-buc=360���ӣ�������ΪS�任������˥��Ϊ�߽�5%�����ݵ㼴2880*5%=144���ӣ��ⲿ�ֽ�����ɿ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(ftdata);%���ݳ���(����)
m=floor((n-chc)/buc)+1;%��������
cs=(chc-2*tdc)/djs;
ST=zeros(fm+1,m*cs);
disp(['��Ҫ����',num2str(m),'��']);
for i=1:1:m
    h=ftdata(1+(i-1)*buc:chc+(i-1)*buc);
    [st_h,~,Fre] = stm1(h,0,fm);
    tmp=st_h(:,tdc+1:chc-tdc);
    tmpa=reshape(tmp,fm+1,djs,cs);
    ST(:,1+(i-1)*cs:i*cs)=reshape(max(tmpa,[],2),fm+1,cs);
%     for jj=1:1:cs
%         ST(:,jj+(i-1)*cs)=max(tmp(:,1+(jj-1)*djs:djs*jj),[],2);%���ز���ʱ�������ֵ��Ϊ�ز������Ƶ�ʳɷֶ�Ӧֵ
%     end
    %disp(['��',num2str(i),'�����']);
end
[ymax,maxi]=max(abs(ST));%ÿʱ�̵������������Ӧ��Ƶ������
ymean=mean(abs(ST));%ÿʱ�̵�ƽ������
Tim=unique(floor(fttimej/100));
Tim=Tim(tdc/60+1:djs/60:end);%
Tim=Tim(1:length(ymax));
TimQS=unique(floor(ttimej(index2)/100));%����ȱ����ʱ��
Cyc=1./Fre;%�ź����ڣ���λͬ�����ʵ�λ
disp('S�任���');
%�����������ݼ����ֵĳ�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end