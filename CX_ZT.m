%Flag=1ֻ�����ʣ�Flag=2ֻʱƵ��Flag=3ͬʱ�����ʺ�ʱƵ
%FileName,PathName
%NFZ�ļ�����
%QSȱ�����
%FlagStep������
%daterange���ݵ�ʱ�䷶Χ
%%
function daterange=CX_ZT(Flag,FileName,PathName,NFZ,QS,FlagStep)
load('Station.mat');%���ڲ�ѯ̨վ�����Ϣ
nEvent=[];
jwdEvent=[];
zjEvent=[];

dname1='�����ʼ�ʱƵ����-����';
if exist(dname1,'dir')~=7
    mkdir(dname1);
end
dname2='�����ʼ�ʱƵ����-ͼ��';
if exist(dname2,'dir')~=7
    mkdir(dname2);
end

%% �����˲������˲���ǰ��ͨ��FirLowFilter.m������Ʋ���
load('firlowfilterW.mat','b');
delelen=2;%ɾ���˲������������ǰ���delelen���ݣ���λ�죬һ��Ϊ�˲������ȵ�һ��
%% ����ļ����д���
if NFZ==1%һ���ļ������ļ���ȡ��Щ���
    FileName={FileName};
end
ALLTime=0;
for iiNFZ=1:1:NFZ
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    STime=0;
    %��������
    tic
    FF=FileName{iiNFZ};
    dbfile=[PathName,FF];
    disp(['�ļ�����:',num2str(NFZ),' Ŀǰ���ڴ����',num2str(iiNFZ),'�ļ���Ϊ��',FileName{iiNFZ}]);
    drsj=load(dbfile); %��������
    [~,N]=size(drsj);
    %��������������ݣ��������ļ�
    if N~=2
        continue;
    end
    %���������������ļ�
    ttimej=drsj(:,1);
    tdata=drsj(:,2);
    if length(FF)>=12
        tkkx=strmatch(strcat(FF(1:5),FF(7)),[TZDM,CDBH]);%�ж�վ����Ϣ�Ƿ����
        if ~isempty(tkkx)
            jwd=JWD(tkkx(1),:);
            wname=deblank(TZM(tkkx(1),:));
        else
            jwd=[NaN NaN];
            wname='';
        end
    else
        jwd=[NaN NaN];
        wname='';
    end
    daterange(iiNFZ,1)=ttimej(1);
    daterange(iiNFZ,2)=ttimej(end);
    timeDR=toc;
    disp(['�������ݺ�ʱ:',num2str(timeDR),'��']);
    STime=STime+timeDR;
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
    %˵����tdataԭʼ���ݣ�stdataȥ̨�����ݣ�qtdataȥ3�������ݣ�itdata��ֵ�����ݣ�ftdata��ֵ���˲����ݣ�ttimejʱ��,fttimejȥ�����ɿ����ݵ�ʱ��
    tic
    [tdata,ttimej]=FillGap(tdata,ttimej,QS); %�����,�ļ�����Ȼ����ȱ�����
    kswz=find(mod(ttimej,10000)==0,1);%��0��0�����ݿ�ʼ
    jswz=find(mod(ttimej,10000)==2359,1,'last');
    tdata=tdata(kswz:jswz);
    ttimej=ttimej(kswz:jswz);
    if FlagStep==1
        stdata=EraStep(tdata,ttimej,QS,4,0);%�չ��㣬����̨�׺���
    else
        stdata=tdata;%������
    end
    tmpx=1:1:length(stdata);
    tmpy=stdata;
    indd=(stdata==QS);
    tmpx(indd)=[];
    tmpy(indd)=[];
    pp=polyfit(tmpx',tmpy,3);
    tmpyy=polyval(pp,tmpx);%������
    %3�׶���ʽ���ȥ���ƣ����ֵͽ����������Ƶ��й©���������ܻ�Ӱ��Ŀ��Ƶ�������Ҫ���˲�ǰ��ȥ��
    qtdata=stdata;
    qtdata(tmpx)=qtdata(tmpx)-tmpyy';
    %���Բ�ֵ����ȥͻ��
    [itdata,index2]=RepInvalid(qtdata,QS,2);
    if length(index2)==length(tdata)%ȫ��Ϊȱ��
        continue;
    end
    %�˲�,ȥ�����ɿ�����
    ftdata=itdata-filtfilt(b,1,itdata);
    fttimej=ttimej(delelen*1440+1:length(ttimej)-delelen*1440);
    ftdata=ftdata(delelen*1440+1:length(ttimej)-delelen*1440);
    timeLB=toc;
    disp(['�����˲���ʱ:',num2str(timeLB),'��']);
    STime=STime+timeLB;
    %% ���㳬����
    if Flag==1||Flag==3
        tic
        %���㳬����
        [zQDcxl,zSLcxl,timecx,timeqs,bs,~,~]=ComCXL(ftdata,ttimej,fttimej,index2);
        %title1=strsplit(FileName{iiNFZ},'.');
        %title1=[title1{1},'_',wname];
        FF=FileName{iiNFZ};
        f_nn=find(FF=='.',1,'last')-1;
        title1=[FF(1:f_nn),'_',wname];
        %% ��ͼ
        hf=ImageCXL(ttimej,stdata,fttimej,ftdata,timecx,timeqs,zQDcxl,zSLcxl,QS,nEvent,jwdEvent,zjEvent,jwd,bs);
        saveas(hf,[pwd,'\',dname2,'\',title1,'_cxl'],'tif');
        close(hf);
        timeCX=toc;
        disp(['�����ʺ�ʱ:',num2str(timeCX),'��']);
        %��������
        fm=['%',num2str(length(num2str(timecx(1)))),'i %.5f\n'];
        outname=strcat(pwd,'\',dname1,'\',title1,'_zQDcxl','.txt');
        zQDcxl(isnan(zQDcxl))=QS;
        fidof=fopen(outname,'wt');
        fprintf(fidof,fm,[timecx';zQDcxl]);
        fclose(fidof);
        outname=strcat(pwd,'\',dname1,'\',title1,'_zSLcxl','.txt');
        zSLcxl(isnan(zSLcxl))=QS;
        fidof=fopen(outname,'wt');
        fprintf(fidof,fm,[timecx';zSLcxl]);
        fclose(fidof);
        STime=STime+timeCX;
    end    
    %% ʱƵ����
    if Flag==2||Flag==3
        tic
        %S�任
        [ST,ymax,maxi,ymean,Tim,TimQS,Cyc,Fre]=CommainS(ftdata,ttimej,fttimej,index2);
        %title1=strsplit(FileName{iiNFZ},'.');
        %title1=[title1{1},'_',wname];
        FF=FileName{iiNFZ};
        f_nn=find(FF=='.',1,'last')-1;
        title1=[FF(1:f_nn),'_',wname];
        %% ��ͼ
        hf=GmaxS(ttimej,stdata,fttimej,ftdata,Tim,TimQS,nEvent,jwdEvent,zjEvent,jwd,QS,ST,ymax,ymean,Cyc);
        saveas(hf,[pwd,'\',dname2,'\',title1,'_ST'],'tif');
        close(hf);
        timeST=toc;
        disp(['S�任��ʱ:',num2str(timeST),'��']);
        %��������
        fm=['%',num2str(length(num2str(Tim(1)))),'i %.5f\n'];
        outname=strcat(pwd,'\',dname1,'\',title1,'_STmax','.txt');
        ymax(isnan(ymax))=QS;
        fidof=fopen(outname,'wt');
        fprintf(fidof,fm,[Tim';ymax]);
        fclose(fidof);
        outname=strcat(pwd,'\',dname1,'\',title1,'_STmean','.txt');
        ymean(isnan(ymean))=QS;
        fidof=fopen(outname,'wt');
        fprintf(fidof,fm,[Tim';ymean]);
        fclose(fidof);
        STime=STime+timeST;
    end
    disp(['���ļ���ʱ:',num2str(STime),'��']);
    ALLTime=ALLTime+STime;
end
disp(['�ܺ�ʱ:',num2str(ALLTime/60),'����']);
end