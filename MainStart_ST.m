%数据时频分析处理，2015-11-29，刘琦
%说明：
%为了使单个功能更加独立化模块化，单个功能都是从数据读取开始，到输出图件或数据结束，因此由于重复读取数据会导致一定的时间损耗。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%基本参数
clear;clc;close all;
QS=999999;%缺数标记
iflag=3;
%超限率及时频分析:1只超限率；2只时频；3同时超限率和时频
FlagStep=1;%归零
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%导入数据的路径
[FileName,PathName] = uigetfile({'*.txt';'*.dat';'*.*'},'Select the files','MultiSelect', 'on');
if iscell(FileName)
    NFZ=length(FileName);
elseif FileName==0  %如果没有打开文件，则跳出程序
    return;
else
    NFZ=1;
end     %NFZ为待处理文件个数

if iflag~=0
    clc;
    disp('开始超限率、时频分析：')
    daterange=CX_ZT(iflag,FileName,PathName,NFZ,QS,FlagStep);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
