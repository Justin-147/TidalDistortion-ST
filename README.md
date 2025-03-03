# 基于S变换研究潮汐畸变的程序
# Program for studying tidal distortion based on S-transform (TidalDistortion-ST)
## Language

- [English](#english)
- [中文](#中文)

---

## English
### Start MATLAB and run MainStart_ST.m. A selection dialog will pop up for choosing input data. After selection, click OK to proceed. The program supports batch processing of multiple files and will automatically create two folders in the program directory: "Over-limit Rate and Time-Frequency Analysis - Figures" (for storing result plots) and "Over-limit Rate and Time-Frequency Analysis - Data" (for storing result data).
### Input data must be two-column minute-value formatted data.
### Workflow: 
1. Perform basic data preprocessing (including gap filling, daily zeroing, third-order polynomial detrending, and linear interpolation).
2. Apply digital filtering to extract high-pass information.
3. Calculate intensity over-limit rate (zQDcxl) and quantity over-limit rate (zSLcxl) by removing a specified proportion of outliers and performing statistical analysis under given thresholds.
4. Generate and save plots, and export result data. 
5. Perform S-transform on digitally filtered data using overlapping sliding windows with specified window lengths and step sizes.
6. Extract maximum amplitude values and average amplitude values across frequencies at given time points from S-transform results, outputting them as time series.
### Author: Liu Qi

---

## 中文
### 启动matlab，运行MainStart_ST.m，会弹出选择框选择输入数据，选择后点击确定即可。支持多文件批量选择处理，程序会自动在程序所在目录下新建“超限率及时频分析-图件”文件夹（用于存放结果图件）和“超限率及时频分析-数据”文件夹（用于存放结果数据）。
### 输入数据需要是两列格式的分钟值数据。
### 程序流程：
1. 进行简单的数据预处理（包括补断数、日归零、三阶多项式去趋势、线性插值）；
2. 进行数字滤波提取高通信息；
3. 计算强度超限率zQDcxl和数量超限率zSLcxl（即去掉一定比例的坏点后在给定阈值下进行强度和数量的统计）；
4. 成图保存图件和输出结果数据
5. 基于数字滤波后的数据可以进行S变换，选择的方式是给定窗长步长进行重叠滑动；
6. 基于S变换结果在给定时刻下提取各频率最大幅值和平均幅值，作为时间序列输出。
### 作者：刘琦

---
