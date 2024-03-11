

%     .................................................
%             ____  _       _   ____  _____   _        
%            |  _ \| |     |_| |  _ \|  ___| |_|       
%            | |_) | |___   _  | |_) | |___   _        
%            |  _ /|  _  | | | |  _ /|___  | | |       
%            | |   | | | | | | | |    ___| | | |       
%            |_|   |_| |_| |_| |_|   |_____| |_|       
%     .................................................
%     PhiPsi:     a general-purpose computational      
%                 mechanics program written in Fortran.
%     Website:    http://phipsi.top                    
%     Author:     Fang Shi  
%     Contact me: shifang@ustc.edu.cn     

%-------------------------------------------------------------------
%--------------------- PhiPsi_Post_Plot ----------------------------
%-------------------------------------------------------------------

%---------------- Start and define global variables ----------------
clear all; close all; clc; format compact;  format long;
global Key_Dynamic Version Num_Gauss_Points 
global Filename Work_Dirctory Full_Pathname num_Crack
global Num_Processor Key_Parallel Max_Memory POST_Substep
global tip_Order split_Order vertex_Order junction_Order    
global Key_PLOT Key_POST_HF Num_Crack_HF_Curves
global Plot_Aperture_Curves Plot_Pressure_Curves Num_Step_to_Plot
global Key_TipEnrich Elements_to_be_Ploted Key_Plot_Elements
global Key_Plot_EleGauss EleGauss_to_be_Ploted Crack_to_Plot

% Number of Gauss points of enriched element (default 64) for integral solution 2.
Num_Gauss_Points = 64;       
Filename=[];
%-------------------------- Settings -------------------------------
% Add path of source files.
addpath('src_fcw')
addpath('src_geom3d')
addpath('src_meshes3d')
addpath('src_phipsi_post_animate')
addpath('src_phipsi_post_cal')
addpath('src_phipsi_post_main')
addpath('src_phipsi_post_plot')
addpath('src_phipsi_post_read')
addpath('src_phipsi_post_tool')

% Set default figure colour to white.
set(0,'defaultfigurecolor','w')

% Set default figure visible off.
set(0,'DefaultFigureVisible','off')

% Output information of matlab command window to log file.
diary('Command Window.log');        
diary on;
% rotate3d on;
% Display welcome information.

set(gcf,'renderer','opengl')

Welcome  
       
tic;
Tclock=clock;
Tclock(1);

disp([' >> Start time is ',num2str(Tclock(2)),'/',num2str(Tclock(3)),'/',num2str(Tclock(1))...
     ,' ',num2str(Tclock(4)),':',num2str(Tclock(5)),':',num2str(round(Tclock(6))),'.'])
disp(' ') 

% Make the "patch" method supported by "getframe", added in version 4.8.10
% See more: http://www.mathworks.com/support/bugreports/384622
% opengl('software') 
opengl hardware     

%----------------------- Pre-Processing ----------------------------
disp(' >> Reading input file....') 

% -------------------------------------
%   Set color and font
% -------------------------------------                           
PhiPsi_Color_and_Font_Settings    

% -------------------------------------
%   Start Post-processor.      
% -------------------------------------   
Key_PLOT   = zeros(7,15);                                   % Initialize the Key_PLOT

%###################################################################################
%                                                                                  #
%                                                                                  #
%                              User defined part.                                  #
%                                                                                  #
%                                                                                  #
%###################################################################################
%----------------------   phipsi_verification  --------------------------
Filename='hydraulic_fracturing_3d';Work_Dirctory='X:\PhiPsi_Project\phipsi_verification\hydraulic_fracturing_3d';

Num_Step_to_Plot      =  -999       ;       %-999  %后处理结果计算步号
Crack_to_Plot         = 0         ;       %仅绘制给定的裂缝，=0则绘制全部裂缝. 支持3D. 2022-09-28.
% Defor_Factor          = 5.0e0      ;       %变形放大系数
Defor_Factor          = 1.0e0      ;       %变形放大系数
% Defor_Factor          = 0.0000001      ;       %变形放大系数
% Defor_Factor          = 1.0e-180   ;       %变形放大系数

% 第1行,有限元网格: Plot Mesh(=1),
%                   2:Node(=1,only nodes;=2,plot nodes + node number)
%                   3:Element(=1,plot solid element;=2,plot only element number;=3,plot element number + element center),Gauss points(4),
%                   5: 裂缝及裂缝坐标点等
%                      =1,绘制裂缝面;                =2,绘制裂缝面及坐标点;           =3,绘制裂缝面及坐标点编号;
%                      =4,绘制裂缝面离散线;          =5,绘制裂缝面及离散单元编号),
%                   6: 计算点(流体节点)及其编号
%                      =1,流体单元和节点;     =2,流体单元和节点+节点编号;    =3,流体单元和节点+单元编号;  =4,流体单元及编号+计算点及编号
%                   7: 流体单元计算点以及Gauss点相关信息
%                      = 1, 增强节点载荷(矢量);        = 2, 计算点净水压(矢量);     = 3, 计算点流量(矢量);       = 4, 计算点开度(矢量);
%                      =11, 面力对应的节点载荷(矢量);  =12, 流体单元净水压(云图);   =13, 流体单元流量(云图);     =14, 流体单元开度(云图),
%                      =21, 流体单元Gauss点(单点积分)外法线向量;    =22, 流体单元Gauss点(单点积分)局部坐标系;    =25, 流体单元Gauss点(单点积分)接触力
%                   8:增强节点(=1,增强节点;               =2,增强节点+增强单元(纯增强单元)网格;     =3,增强节点+增强节点对应的FEM节点号),
%                   9:网格线(=1,全部网格线),10:Blank,
%                   11:绘制指定单元及其节点号(结束单元号),
%                   12:全部天然裂缝(=1,天然裂缝;=2,天然裂缝+编号),单元接触状态(13),裂缝编号(14,=1),15(仅绘制给定材料的单元)
% 第2行,网格变形图: Deformation(1),Node(2),El(3),Gauss points(4)
%                   5:Crack Shape. =1,fracture area; =2,fracture volume
%                   6:Scaling Factor.
%                   7:Boundary Conditions. =1 or =2 or =3. 
%                   8:变形前的网格. =1,变形前的全部网格；=2,变形前的网格边界
%                   9:网格线=0,外边界线框;=1,全部网格;=2(=12,+编号),增强单元网格+外边界线框;=3,表面网格+外边界线框;=4(=14,+编号),增强单元网格+表面网格;
%                            =5(=15,+编号),裂尖增强单元(8个节点都是裂尖增强)网格+外边界线框);
%                            =6(=16,+编号),Heaviside增强单元(8个节点都是裂尖增强)网格+外边界线框).
%                   10:=1,裂缝边界离散点的局部坐标系;=2,裂尖增强单元的基准线(baseline);=3,裂尖增强单元的基准线(baseline)+基准局部坐标系;=4,裂尖增强节点对应的单元(参考单元);
%                      =5,裂缝边界离散点及编号
%                   11:=1,裂缝面前缘节点主应力矢量;
%                   12:杂项控制: =1, 绘制裂尖应力计算点,CFCP=2时有效;          =2, 绘制计算球内的高斯点,CFCP=2时有效;
%                                =3, 绘制增强节点的实际位移矢量;               =4, 绘制增强节点的增强自由度位移矢量.
%                   13:单元接触状态
%                   14:增强节点(=1,增强节点)
%                   15:Fracture zone
% 第3行,应力(变)云图: Nodal Stress(Strain) Contour(=1:应力节点值;  =2:x slice;  =3:y slice;  =4:z slice;  =11,节点应变值),
%                   分量(2:=1,Only Mises stress;=2,仅x分量;=3,仅y分量;=4,仅z分量),
%                   主应力(3:1,主应力;2,主应力+方向;3,仅最大主应力),x or y or z slice location(4),
%                   Blank(5),Scaling Factor(6),
%                   undeformed or Deformed(8),mesh(9),
%                   绘图坐标系(10: 1,笛卡尔坐标系; 2,柱坐标系),
%                   Blank(11),Blank(12),Blank(13),Blank(14),Slice location(15)
% 第4行,位移云图:   Nodal Displacement Contour(=1:节点值云图;=2:x slice;=3:y slice;=4:z slice),(2=1,only x;=2,only y;=3,only z;=4,only vector sum, 对Slice不起作用),Blank(3),Blank(4),
%                   Plot crack(5=1),Scaling Factor(6),
%                   Blank(7),Blank(8),mesh(9),
%                   绘图坐标系(10: 1,笛卡尔坐标系; 2,柱坐标系),
%                   Blank(11),Blank(12),Blank(13),Blank(14),Slice location(15)
% 第5行,裂缝云图:   Plot Crack Contour(=1,aperture contour;
%                                      =21,xx渗透率云图,22,yy渗透率云图,23,zz渗透率云图,24,xy渗透率云图,25,yz渗透率云图,26,xz渗透率云图),
%                   2: =1,results of fluid elements and nodes; =2, results of discrete crack surface
%                   3-裂缝面网格线,
%                   4-对于渗透率绘图,=0为国际单位制,=1为mDarcy,
%                   5:编号绘制:=1,绘制裂缝面离散点编号(or流体节点编号); =2,绘制裂缝面离散单元编号(or流体单元编号)
%                   6:Scaling Factor,7-Blank,8-Blank,
%                   9:网格线(=0,外边界线框;=1,全部网格;=2,增强单元网格+外边界线框;=3,表面网格+外边界线框;=4,增强单元网格+表面网格;=5,增强单元网格),Blank(10),
%                   Blank(11),12:全部天然裂缝(=1,天然裂缝;=2,天然裂缝+编号),Blank(13),Blank(14),Blank(15)
% 第6行,曲线绘制相关：开关(1=1,则绘制曲线),
%                   绘制内容(2=1,绘制裂缝前缘Vertex的最大主应力曲线;
%                             =2,绘制裂缝前缘Vertex的等效应力强度因子Keq曲线；
%                             =3,绘制裂缝前缘Vertex的I型应力强度因子KI曲线；
%                             =4,绘制裂缝前缘Vertex的II型应力强度因子KII曲线;
%                             =5,绘制裂缝前缘Vertex的I型应力强度因子KIII曲线;
%                             =6,绘制裂缝前缘Vertex的I型+II型+III型应力强度因子KIII曲线;
%                             =7,绘制HF分析的压力-时间曲线;
%                             =8,绘制HF分析的压力-扩展步曲线);
%                             =9,绘制扩展步-裂缝体积曲线;
%                             =10,绘制压裂实验仿真时间-压力曲线.
%                   裂缝号(3),4:=1绘制光滑处理前的数据,5:井筒号(=-999,绘制全部),Blank(6-15)
% 第7行,其他绘制: 1：=1，绘制单元渗透率(单位为mDarcy)切片云图（注意是单元结果云图）; 2 (=1,x轴切片,=2,y轴切片,=3,z轴切片)；3(切片坐标); 4(=1绘制单元网格线);
%                 5(=1，坐标轴); 6(=1，log图); 7-15(blank)
%      
%             
% Location=  -750;    %Key_PLOT(7,:)绘制Location. 2022-11-29.
% Location=  1;    %Key_PLOT(7,:)绘制Location. 2022-11-29.
% Location=  401.1;    %Key_PLOT(7,:)绘制Location. 2022-11-29.
Location=  172.5;    %Key_PLOT(7,:)绘制Location. 2022-11-29.
%
%                         1    2         3    4    5   6              7    8   9  10      11    12  13  14          15
Key_PLOT(1,:)         = [ 1,   0,        0,   0,   1,  0,             0,   0,  0  ,0     ,0     ,1  ,0  ,1  ,       0];  
Key_PLOT(2,:)         = [ 1,   0,        0,   0,   1,  Defor_Factor,  0,   0,  0  ,1     ,0     ,0  ,0  ,0  ,       0];  
Key_PLOT(3,:)         = [ 0,   1,        0,   0,   0,  Defor_Factor,  0,   1,  0  ,1     ,0     ,0  ,0  ,0  ,Location];  
Key_PLOT(4,:)         = [ 0,   4,        0,   0,   1,  Defor_Factor,  0,   0,  0  ,1     ,0     ,0  ,0  ,0  ,Location];
Key_PLOT(5,:)         = [ 1,   1,        0,   1,   0,  Defor_Factor,  0,   0,  0  ,0     ,0     ,1  ,0  ,0  ,       0];   
Key_PLOT(6,:)         = [ 0,   7,        1,   1,-999,  0,             0,   0,  0  ,0     ,0     ,0  ,0  ,0  ,       0];   
Key_PLOT(7,:)         = [ 0,   3, Location,   1,   0,  1,             0,   0,  0  ,0     ,0     ,0  ,0  ,0  ,       0];   

Key_Plot_Elements =  1;   %绘制指定的单元(Elements_to_be_Ploted). 2022-07-31.
Key_Plot_EleGauss =  1;   %绘制指定单元的Gauss积分点(仅适用于固定积分算法. *Key_Integral_Sol  = 2 ). 2022-07-31.
% Elements_to_be_Ploted = [1,2,3];   %Key_PLOT(1,10)=-999时激活.
% Elements_to_be_Ploted = [1,2,3];   %Key_PLOT(1,10)=-999时激活.
% Elements_to_be_Ploted = [1238,1463,1688,1913,2138];   %Key_PLOT(1,10)=-999时激活.
% Elements_to_be_Ploted = [1732,1736];   %Key_PLOT(1,10)=-999时激活.
% Elements_to_be_Ploted = [1415,1416,1417,1418,1419,1420,1430,1431,1432,1433,1434,1435,1445,1446,1447,1448,1449,1450,1460,1461,1462,1463,1464,1465,...
                         % 1475,1476,1477,1478,1479,1480,1490,1491,1492,1493,1494,1495,1640,1641,1642,1643,1644,1645,1655,1656,1657,1658,1659,1660,...
						 % 1670,1671,1672,1673,1674,1675,1685,1686,1687,1688,1689,1690,1700,1701,1702,1703,1704,1705,1715,1716,1717,1718,1719,1720];   %Key_PLOT(1,10)=-999时激活.
% Elements_to_be_Ploted = [];
% Elements_to_be_Ploted = [5110];
% Elements_to_be_Ploted = [1463,1464,1478,1479,1688,1689,1703,1704];
% Elements_to_be_Ploted = [5802,5835];
% Elements_to_be_Ploted = [9334,7234,6609,5984,      10361,9736,9111,8486,7861,7236,6611,6986,5361];
% Elements_to_be_Ploted = [7234,   7236];
% Elements_to_be_Ploted = [7861];
% Elements_to_be_Ploted = [1703,1899,1913,1927];
% Elements_to_be_Ploted = [1703];

% Elements_to_be_Ploted = [1237,2139];
% EleGauss_to_be_Ploted = [1237,2139];

% Elements_to_be_Ploted = [983,1206,1251,2125,2170,2393];
% EleGauss_to_be_Ploted = [983,1206,1251,2125,2170,2393];


% Elements_to_be_Ploted = [1237,2139,1223,2153];
% EleGauss_to_be_Ploted = [1237,2139,1223,2153];

% Elements_to_be_Ploted = [2163];
% EleGauss_to_be_Ploted = [2163];

% Elements_to_be_Ploted = [1925,1924,1926];
% EleGauss_to_be_Ploted = [1925,1924,1926];

% Elements_to_be_Ploted = [1614,1615];
% EleGauss_to_be_Ploted = [1614,1615];

% Elements_to_be_Ploted = [2166,2406,2408];
% EleGauss_to_be_Ploted = [2166,2406,2408];

% Elements_to_be_Ploted = [635,2735];
% Elements_to_be_Ploted = [5788,5793,4111,4116,8614,8619,10759,10764];
% Elements_to_be_Ploted = [3918,3917,3916];
% Elements_to_be_Ploted = [3486,3487,3488,3489,3490,3491];

% Elements_to_be_Ploted = [7234,7235,7236,7237,7238,7862,7863];
% EleGauss_to_be_Ploted = [7234,7235,7236,7237,7238,7862,7863];

% Elements_to_be_Ploted = [2138,2139,2153,2154,1913,1914,1928,1929];
% EleGauss_to_be_Ploted = [2138,2139,2153,2154,1913,1914,1928,1929];

% Elements_to_be_Ploted = [1388,1389,1403,1404,1613,1614,1628,1629];
% EleGauss_to_be_Ploted = [1388,1389,1403,1404,1613,1614,1628,1629];

% Elements_to_be_Ploted = [1521];
% EleGauss_to_be_Ploted = [1521];

% Elements_to_be_Ploted = [212,213,221,222,131,132,140,141];
% EleGauss_to_be_Ploted = [212,213,221,222,131,132,140,141];

% Elements_to_be_Ploted = [464,465,374  ,375 ,   383    ,384    ,455  ,  456];
% EleGauss_to_be_Ploted = [464,465,374  ,375 ,   383    ,384    ,455  ,  456];

% Elements_to_be_Ploted = [1459  , 1460,   1474 ,  1475  , 1234 ,  1235 ,  1249   ,1250];
% EleGauss_to_be_Ploted = [1459  , 1460,   1474 ,  1475  , 1234 ,  1235 ,  1249   ,1250];

% Elements_to_be_Ploted = [881];
% EleGauss_to_be_Ploted = [881];

% Elements_to_be_Ploted = [8386  ,      8387  ,      8388 ,       8411   ,     8412    ,    8413, ...
						 % 8436,        8437  ,      8438 ,       8461 ,       8462  ,      8463,...
						% 8486  ,      8487  ,      8414  ,      8415 ,       8416  ,      8439,...
						% 8440  ,      8441  ,      8464  ,      8465 ,       8466  ,      8488,...
						% 8489  ,      8490  ,      8491,        8514 ,       8515  ,      8485,...
						% 8510  ,      8511  ,      8512  ,      8513,        8539 ,       8540];

% EleGauss_to_be_Ploted = [1656, 1657,1658, 1659, 1660];
% Elements_to_be_Ploted = [36705];

% Elements_to_be_Ploted = [36158];
% EleGauss_to_be_Ploted = [191643];
%###################################################################################
%                                                                                  #
%                                                                                  #
%                        End of user defined part.                                 #    
%                                                                                  #
%                                                                                  #
%###################################################################################

PhiPsi_Post_1_Go_3D

%退出程序.
% exit
