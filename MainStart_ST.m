%����ʱƵ��������2015-11-29������
%˵����
%Ϊ��ʹ�������ܸ��Ӷ�����ģ�黯���������ܶ��Ǵ����ݶ�ȡ��ʼ�������ͼ�������ݽ�������������ظ���ȡ���ݻᵼ��һ����ʱ����ġ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��������
clear;clc;close all;
QS=999999;%ȱ�����
iflag=3;
%�����ʼ�ʱƵ����:1ֻ�����ʣ�2ֻʱƵ��3ͬʱ�����ʺ�ʱƵ
FlagStep=1;%����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�������ݵ�·��
[FileName,PathName] = uigetfile({'*.txt';'*.dat';'*.*'},'Select the files','MultiSelect', 'on');
if iscell(FileName)
    NFZ=length(FileName);
elseif FileName==0  %���û�д��ļ�������������
    return;
else
    NFZ=1;
end     %NFZΪ�������ļ�����

if iflag~=0
    clc;
    disp('��ʼ�����ʡ�ʱƵ������')
    daterange=CX_ZT(iflag,FileName,PathName,NFZ,QS,FlagStep);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
