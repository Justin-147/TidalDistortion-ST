%Flag=1只超限率；Flag=2只时频；Flag=3同时超限率和时频
%FileName,PathName
%NFZ文件个数
%QS缺数标记
%FlagStep归零标记
%daterange数据的时间范围
%%
function daterange=CX_ZT(Flag,FileName,PathName,NFZ,QS,FlagStep)
load('Station.mat');%用于查询台站相关信息
nEvent=[];
jwdEvent=[];
zjEvent=[];

dname1='超限率及时频分析-数据';
if exist(dname1,'dir')~=7
    mkdir(dname1);
end
dname2='超限率及时频分析-图件';
if exist(dname2,'dir')~=7
    mkdir(dname2);
end

%% 读入滤波器，滤波器前期通过FirLowFilter.m进行设计测试
load('firlowfilterW.mat','b');
delelen=2;%删除滤波器处理后数据前后各delelen数据，单位天，一般为滤波器长度的一半
%% 逐个文件进行处理
if NFZ==1%一个文件与多个文件读取有些差别
    FileName={FileName};
end
ALLTime=0;
for iiNFZ=1:1:NFZ
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    STime=0;
    %导入数据
    tic
    FF=FileName{iiNFZ};
    dbfile=[PathName,FF];
    disp(['文件总数:',num2str(NFZ),' 目前正在处理第',num2str(iiNFZ),'文件名为：',FileName{iiNFZ}]);
    drsj=load(dbfile); %导入数据
    [~,N]=size(drsj);
    %如果不是两列数据，则跳过文件
    if N~=2
        continue;
    end
    %读入了整个数据文件
    ttimej=drsj(:,1);
    tdata=drsj(:,2);
    if length(FF)>=12
        tkkx=strmatch(strcat(FF(1:5),FF(7)),[TZDM,CDBH]);%判断站点信息是否存在
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
    disp(['导入数据耗时:',num2str(timeDR),'秒']);
    STime=STime+timeDR;
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
    %说明：tdata原始数据，stdata去台阶数据，qtdata去3阶项数据，itdata插值后数据，ftdata插值后滤波数据，ttimej时间,fttimej去掉不可靠数据的时间
    tic
    [tdata,ttimej]=FillGap(tdata,ttimej,QS); %填补断数,文件中仍然存在缺数标记
    kswz=find(mod(ttimej,10000)==0,1);%从0点0分数据开始
    jswz=find(mod(ttimej,10000)==2359,1,'last');
    tdata=tdata(kswz:jswz);
    ttimej=ttimej(kswz:jswz);
    if FlagStep==1
        stdata=EraStep(tdata,ttimej,QS,4,0);%日归零，日内台阶忽略
    else
        stdata=tdata;%不归零
    end
    tmpx=1:1:length(stdata);
    tmpy=stdata;
    indd=(stdata==QS);
    tmpx(indd)=[];
    tmpy(indd)=[];
    pp=polyfit(tmpx',tmpy,3);
    tmpyy=polyval(pp,tmpx);%三阶项
    %3阶多项式拟合去趋势，这种低阶项会在整个频带泄漏能量，可能会影响目标频带，因此要在滤波前先去除
    qtdata=stdata;
    qtdata(tmpx)=qtdata(tmpx)-tmpyy';
    %线性插值，不去突跳
    [itdata,index2]=RepInvalid(qtdata,QS,2);
    if length(index2)==length(tdata)%全部为缺数
        continue;
    end
    %滤波,去掉不可靠数据
    ftdata=itdata-filtfilt(b,1,itdata);
    fttimej=ttimej(delelen*1440+1:length(ttimej)-delelen*1440);
    ftdata=ftdata(delelen*1440+1:length(ttimej)-delelen*1440);
    timeLB=toc;
    disp(['数据滤波耗时:',num2str(timeLB),'秒']);
    STime=STime+timeLB;
    %% 计算超限率
    if Flag==1||Flag==3
        tic
        %计算超限率
        [zQDcxl,zSLcxl,timecx,timeqs,bs,~,~]=ComCXL(ftdata,ttimej,fttimej,index2);
        %title1=strsplit(FileName{iiNFZ},'.');
        %title1=[title1{1},'_',wname];
        FF=FileName{iiNFZ};
        f_nn=find(FF=='.',1,'last')-1;
        title1=[FF(1:f_nn),'_',wname];
        %% 画图
        hf=ImageCXL(ttimej,stdata,fttimej,ftdata,timecx,timeqs,zQDcxl,zSLcxl,QS,nEvent,jwdEvent,zjEvent,jwd,bs);
        saveas(hf,[pwd,'\',dname2,'\',title1,'_cxl'],'tif');
        close(hf);
        timeCX=toc;
        disp(['超限率耗时:',num2str(timeCX),'秒']);
        %保存数据
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
    %% 时频分析
    if Flag==2||Flag==3
        tic
        %S变换
        [ST,ymax,maxi,ymean,Tim,TimQS,Cyc,Fre]=CommainS(ftdata,ttimej,fttimej,index2);
        %title1=strsplit(FileName{iiNFZ},'.');
        %title1=[title1{1},'_',wname];
        FF=FileName{iiNFZ};
        f_nn=find(FF=='.',1,'last')-1;
        title1=[FF(1:f_nn),'_',wname];
        %% 画图
        hf=GmaxS(ttimej,stdata,fttimej,ftdata,Tim,TimQS,nEvent,jwdEvent,zjEvent,jwd,QS,ST,ymax,ymean,Cyc);
        saveas(hf,[pwd,'\',dname2,'\',title1,'_ST'],'tif');
        close(hf);
        timeST=toc;
        disp(['S变换耗时:',num2str(timeST),'秒']);
        %保存数据
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
    disp(['此文件耗时:',num2str(STime),'秒']);
    ALLTime=ALLTime+STime;
end
disp(['总耗时:',num2str(ALLTime/60),'分钟']);
end