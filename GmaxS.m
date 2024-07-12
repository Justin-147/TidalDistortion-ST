function hf=GmaxS(ttimej,stdata,fttimej,ftdata,Tim,TimQS,nEvent,jwdEvent,zjEvent,jwd,QS,ST,ymax,ymean,Cyc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
bsb3=2;%ͼ4y�᷽���޶�
cm=[];
if isempty(cm)
    tmpst=abs(ST(:));
    mst=max(tmpst);
    cm(1)=mean(tmpst)-bsb2*std(tmpst);
    cm(2)=mean(tmpst)+bsb2*std(tmpst);
    if cm(1)<0
        cm(1)=0;
    end
    if cm(2)>mst
        cm(2)=mst;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hf=figure();
set(hf,'position',fps);
set(hf,'PaperPositionMode','auto');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%������
tt3=datenum([floor(Tim/10^6),mod(floor(Tim/10^4),100),mod(floor(Tim/10^2),100),mod(Tim,100),zeros(length(Tim),2)]);
[~,iin,~]=intersect(Tim,TimQS);
ST(:,iin)=0;
subplot(4,1,3)
hold on;
imagesc(tt3,1:size(ST,1),abs(ST));
set(gca,'ydir','normal');
box off;
colormap('jet');
if ~isempty(cm)
    caxis(cm);
end
yt=49:120:409;
ytl=round(Cyc(yt));
set(gca,'ytick',yt);
set(gca,'yticklabel',ytl);
xlim([tt(1) tt(end)]);
ylim([1 size(ST,1)]);
hold off;
tlabel();
set(gca,'tickdir','out');
set(gca,'xticklabel',[]);
set(gca,'FontSize',FSL,'FontWeight',FWL);
set(gca,'Position',[zwz xwz+tjg kbl gbl]);
ylabel('�ź�����/����','FontSize',FSL,'FontWeight',FWL);
hbar=colorbar;
set(hbar,'position',[0.9203 0.3256 0.0115 0.1580]);
htext=annotation('textbox','string',['����',char(13),char(10),'/�۲ⵥλ']);
set(htext,'Position',[0.9293 0.5051 0.0535 0.0450],'linestyle','none','FontSize',FSL-2,'FontWeight',FWL,'HorizontalAlignment','center');
%������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���Ķ�
subplot(4,1,4)
cc='kr';
hold on;
tmpstem=[ymax;ymean];
tmpstem(:,iin)=0;
for kk=1:1:size(tmpstem,1)
    stem(tt3,tmpstem(kk,:)','linewidth',2.6-kk*0.4,'color',cc(kk),'marker','none');
end
xlim([tt(1) tt(end)]);
ttmpstem=tmpstem(:);
%ylim([0 max(ttmpstem)]);
ym=[0 mean(ttmpstem)+bsb3*std(ttmpstem)];
if ym(1)~=ym(2)
    ylim(ym);
end
hold off;
tlabel();
set(gca,'tickdir','out');
set(gca,'FontSize',FSL,'FontWeight',FWL);
set(gca,'Position',[zwz xwz kbl gbl]);
LET=['������';'ƽ������'];
legend(LET,'location','best');
ylabel('S�任��ֵ','FontSize',FSL,'FontWeight',FWL);
xlabel('����','FontSize',FSL,'FontWeight',FWL);
box off;
%���Ķ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
