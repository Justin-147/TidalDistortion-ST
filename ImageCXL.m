function hf=ImageCXL(ttimej,stdata,fttimej,ftdata,timecx,timeqs,zQDcxl,zSLcxl,QS,nEvent,jwdEvent,zjEvent,jwd,bs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������
FSL=12;
FWL='Normal';
fps=[294 125 764 722];%figure��Χ
kbl=0.80;%ÿͼ�����
gbl=0.20;%ÿͼ�߱���
zwz=0.115;%ͼ�����
xwz=0.075;%ͼ������
tjg=(1-xwz)/4;

bsb1=4;%ͼ2y�᷽���޶�
bsb2=2;%ͼ3y�᷽���޶�
bsb3=4;%ͼ4y�᷽���޶�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hf=figure();
set(hf,'position',fps);
set(hf,'PaperPositionMode','auto');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��һ��
tt=datenum(floor(ttimej/10^8),mod(floor(ttimej/10^6),100),mod(floor(ttimej/10^4),100),mod(floor(ttimej/100),100),mod(ttimej,100),zeros(length(ttimej),1));
stdata(stdata==QS)=NaN;
subplot(4,1,1)
plot(tt,stdata,'r','linewidth',0.75);
xlim([tt(1) tt(end)]);
ym=[min(stdata) max(stdata)];
if ym(1)~=ym(2)
    ylim(ym);
end
tlabel();
set(gca,'tickdir','out');
set(gca,'xticklabel',[]);
set(gca,'FontSize',FSL,'FontWeight',FWL);
set(gca,'Position',[zwz xwz+3*tjg kbl gbl]);
ylabel('ԭʼ����/�۲ⵥλ','FontSize',FSL,'FontWeight',FWL);
%��һ��
%�����
dylim=get(gca,'ylim');
if ~isnan(jwd(1))
    if ~isempty(nEvent)
        [~,dz,in]=intersect(ttimej,nEvent);
        if ~isempty(in)
            dis=num2str(floor(vdist(repmat(jwd(2),length(in),1),repmat(jwd(1),length(in),1),jwdEvent(in,2),jwdEvent(in,1))/1000));
            dis=[dis,repmat('km',size(dis,1),1)];
            MM=num2str(zjEvent(in));
            ddy=(dylim(2)-dylim(1))/20;
            hold on
            yy1=repmat(dylim(2),length(dz),1);
            plot(tt(dz),yy1-0.5*ddy,'rv','markerfacecolor','none','markeredgecolor','b','Markersize',7);
            yy2=yy1;
            yy2(2:2:end)=yy2(2:2:end)-ddy;
            yy3=yy1-2*ddy;
            yy3(2:2:end)=yy3(2:2:end)-ddy;
            text(tt(dz),yy2,MM,'FontSize',FSL-3,'FontWeight',FWL);
            text(tt(dz),yy3,dis,'FontSize',FSL-3,'FontWeight',FWL,'rotation',90,'verticalalignment','middle','horizontalalignment','right');
            hold off
        end
    end
end
box off;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�ڶ���
tt2=datenum(floor(fttimej/10^8),mod(floor(fttimej/10^6),100),mod(floor(fttimej/10^4),100),mod(floor(fttimej/100),100),mod(fttimej,100),zeros(length(fttimej),1));
fmean=mean(ftdata);
fstd=std(ftdata);
ftdata(ftdata==QS)=NaN;
subplot(4,1,2)
plot(tt2,ftdata,'k','linewidth',0.75);
xlim([tt(1) tt(end)]);
ym=[fmean-bsb1*fstd fmean+bsb1*fstd];
if ym(1)~=ym(2)
    ylim(ym);
end
tlabel();
set(gca,'tickdir','out');
set(gca,'xticklabel',[]);
set(gca,'FontSize',FSL,'FontWeight',FWL);
set(gca,'Position',[zwz xwz+2*tjg kbl gbl]);
ylabel('��Ƶ����/�۲ⵥλ','FontSize',FSL,'FontWeight',FWL);
box off;
%�ڶ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startcx=1;%ֻ�����ֳ�����
LET=num2str(bs');
LET=[LET,repmat('����������ֵ',length(bs),1)];
LET=LET(startcx:end,:);
%������
tt3=datenum(floor(timecx/10^4),mod(floor(timecx/100),100),mod(timecx,100));
%cc='gbkr';
cc='k';
subplot(4,1,3)
hold on;
[~,iin,~]=intersect(timecx,timeqs);
tmpstem=zQDcxl;
tmpstem(:,iin)=0;
for kk=startcx:1:size(tmpstem,1)
    stem(tt3,tmpstem(kk,:)','linewidth',2.6-kk*0.4,'color',cc(kk),'marker','none');
end
xlim([tt(1) tt(end)]);
ttmpstem=tmpstem(startcx:end,:);
ttmpstem=ttmpstem(:);
%ylim([0 max(ttmpstem)]);
ym=[0 mean(ttmpstem)+bsb2*std(ttmpstem)];
if ym(1)~=ym(2)
    ylim(ym);
end
dylim3=get(gca,'ylim');
if ~isnan(jwd(1))
    if ~isempty(nEvent)
        [~,dz3,~]=intersect(timecx,floor(nEvent/10000));
        if ~isempty(dz3)
            yy3=repmat(dylim3(1),length(dz3),1);
            plot(tt3(dz3),yy3,'rv','markerfacecolor','none','markeredgecolor','b','Markersize',7);
        end
    end
end
hold off;
tlabel();
set(gca,'tickdir','out');
set(gca,'xticklabel',[]);
set(gca,'FontSize',FSL,'FontWeight',FWL);
set(gca,'Position',[zwz xwz+tjg kbl gbl]);
haa=legend(LET,'location','best');
set(haa,'FontSize',FSL-2);
ylabel('ǿ�ȳ�����/�۲ⵥλ','FontSize',FSL,'FontWeight',FWL);
box off;
%������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���Ķ�
subplot(4,1,4)
hold on;
tmpstem=zSLcxl;
tmpstem(:,iin)=0;
for kk=startcx:1:size(tmpstem,1)
    stem(tt3,tmpstem(kk,:)','linewidth',2.6-kk*0.4,'color',cc(kk),'marker','none');
end
xlim([tt(1) tt(end)]);
ttmpstem=tmpstem(startcx:end,:);
ttmpstem=ttmpstem(:);
%ylim([0 max(ttmpstem)]);
ym=[0 mean(ttmpstem)+bsb3*std(ttmpstem)];
if ym(1)~=ym(2)
    ylim(ym);
end
dylim3=get(gca,'ylim');
if ~isnan(jwd(1))
    if ~isempty(nEvent)
        if ~isempty(dz3)
            yy3=repmat(dylim3(1),length(dz3),1);
            plot(tt3(dz3),yy3,'rv','markerfacecolor','none','markeredgecolor','b','Markersize',7);
        end
    end
end
hold off;
tlabel();
set(gca,'tickdir','out');
set(gca,'FontSize',FSL,'FontWeight',FWL);
set(gca,'Position',[zwz xwz kbl gbl]);
%legend(LET,'location','best');
ylabel('����������/��','FontSize',FSL,'FontWeight',FWL);
xlabel('����','FontSize',FSL,'FontWeight',FWL);
box off;
%���Ķ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end