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
%--------------------- PhiPsi_Post_Plot2 ---------------------------
%-------------------------------------------------------------------

%---------------- Start and define global variables ----------------
clear all; close all; clc; format compact;  format long;
global Key_Dynamic Version Num_Gauss_Points 
global Filename Work_Dirctory Full_Pathname num_Crack Defor_Factor
global Num_Processor Key_Parallel Max_Memory POST_Substep
global tip_Order split_Order vertex_Order junction_Order    
global Key_PLOT Key_POST_HF Num_Crack_HF_Curves num_Na_Crack
global Plot_Aperture_Curves Plot_Pressure_Curves Plot_Velocity_Curves Num_Step_to_Plot
global Key_TipEnrich Plot_Quantity_Curves Plot_Concentr_Curves Itera_Num
global Na_Crack_X Na_Crack_Y num_Na_Crack Itera_HF_Time Key_HF_Analysis
global num_Hole Hole_Coor   
global num_Circ_Inclusion num_Poly_Inclusion Key_Time_String
global Key_Animation Key_Ani_Ave
global Time_Delay


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
% Display welcome information.
Welcome   
       
tic;
Tclock=clock;
Tclock(1);

disp([' >> Start time is ',num2str(Tclock(2)),'/',num2str(Tclock(3)),'/',num2str(Tclock(1))...
     ,' ',num2str(Tclock(4)),':',num2str(Tclock(5)),':',num2str(round(Tclock(6))),'.'])
disp(' ') 

% Make the "patch" method supported by "getframe", added in version 4.8.10
% See more: http://www.mathworks.com/support/bugreports/384622
opengl('software')      

%----------------------- Pre-Processing ----------------------------
disp(' >> Reading input file....') 
Filename=[];

% -------------------------------------
%   Set color and font
% -------------------------------------                           
PhiPsi_Color_and_Font_Settings

% -------------------------------------
%   Start Post-processor.      
% -------------------------------------   
Key_PLOT   = zeros(6,15);                                   % Initialize the Key_PLOT

%###########################################################################################################
%##########################            User defined part        ############################################
%###########################################################################################################
% Filename='exa_3D_crack';Work_Dirctory='x:\PhiPsi_work\exa_3D_crack';
% Filename='exa_3D_hollow_cylinder';Work_Dirctory='x:\PhiPsi_work\exa_3D_hollow_cylinder';
% Filename='3D_Block_11x11x11';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_11x11x11';
% Filename='3D_Block_Tension';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Tension';
% Filename='3D_Block_21x21x21_Geo_Insitu';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_21x21x21_Geo_Insitu'; %21x21x21
% Filename='exa_3D_block_tension';Work_Dirctory='X:\PhiPsi_Work\exa_3D_block_tension';
% Filename='3D_Block_Tension_Fine';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_Tension_Fine';
% Filename='3D_beam_sifs';Work_Dirctory='X:\PhiPsi_Work\3D_beam_sifs';
% Filename='3D_Beam_ThreePointBending';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Beam_ThreePointBending';
% Filename='paper8_3d_edge_crack';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\paper8_3d_edge_crack';
% Filename='exa_3d_edge_crack';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\exa_3d_edge_crack'; %网站算例
% Filename='3D_Block_21x21x21_Fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\3D_Block_21x21x21_Fixed'; %21x21x21
% Filename='3D_Block_25x25x25_Fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\3D_Block_25x25x25_Fixed'; %25x25x25
% Filename='3D_Block_25x25x25_Fixed_NonUni2';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\3D_Block_25x25x25_Fixed_NonUni2'; %25x25x25, 非均匀2：完全非均匀网格
% Filename='3D_Block_35x35x35_Fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\3D_Block_35x35x35_Fixed'; %35x35x35
% Filename='3D_Block_35x35x35_Irregu_Fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\3D_Block_35x35x35_Irregu_Fixed'; %35x35x35-不规则-两种材料. 2022-10-03.
% Filename='3D_Block_45x45x45_Fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\3D_Block_45x45x45_Fixed'; %45x45x45
% Filename='3D_Block_55x55x55_Fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\3D_Block_55x55x55_Fixed'; %55x55x55
% Filename='3D_Block_75x55x55_Fixed';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_75x55x55_Fixed'; %75x75x55
% Filename='3D_Block_75x75x55_Fixed';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_75x75x55_Fixed'; %75x75x55
% Filename='3D_Block_250x250x250_Fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\3D_Block_250x250x250_Fixed'; %250x250x250，相比25x25x25体积放大1000倍. 2022-09-10.
% Filename='3D_Block_50x1000x1000_Fixed';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\3D_Block_50x1000x1000_Fixed'; %5万单元. 2022-10-04.
% Filename='XA_3D_HF_Three_Layers';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\XA_3D_HF_Three_Layers'; %2023-02-12
% Filename='XA_3D_HF_Three_Layers';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\XA_3D_HF_Three_Layers2'; %2023-02-12
% Filename='XA_3D_HF_Three_Layers';Work_Dirctory='E:\PhiPsi问题待排除-多层压裂-计算到65步终止-2023-02-13\PhiPsi_Project\PhiPsi_Work\XA_3D_HF_Three_Layers'; %2023-02-12
% Filename='XA_3D_HF_Three_Layers_2kx1kx500';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\XA_3D_HF_Three_Layers_2kx1kx500'; %2023-02-19
% Filename='XA_3D_HF_3_Layers_2kx1kx500_1m';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\XA_3D_HF_3_Layers_2kx1kx500_1m'; %2023-02-20. 110万单元.
% Filename='1kx1kx1k_single_HF_Xsite';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\1kx1kx1k_single_HF_Xsite'; %用于和Xsite对比. 116.万单元. 2022-11-09.
% Filename='3D_Block_200x200x345_Fixed';Work_Dirctory='x:\PhiPsi_Project\PhiPsi_work\3D_Block_200x200x345_Fixed'; %2023-04-16. 88.32万单元.
Filename='3D_True_Triaxial';Work_Dirctory='X:\PhiPsi_Project\PhiPsi_work\3D_True_Triaxial'; %真三轴实验，2023-01-20

Num_Step_to_Plot      = -999             ;%后处理结果计算步号或者破裂步号(对于水力压裂分析),若为-999,则全部绘制
Defor_Factor       = 1.0e0                ;%变形放大系数

%-------------------------------
Key_Animation= [0 0 0 1];                % 1:displacement(2,Gauss);2:stress(2,Gauss);3:deformation(裂缝面演化(网格));4:裂缝面云图
Key_Ani_Ave  = [1 1 1 1];                % Key for uniform maximum and minimum
Time_Delay   = 0.15;   %0.2;              % Delay time of each frame of the gif animation,default: 0.025
Key_Time_String = 1;                     % 时间的单位: =1,则为s;=2,min;=3,hour;=4,day;=5,month;=6,year
%-------------------------------
% 第1行,有限元网格: Plot Mesh(=1),
%                   2: Node(=1,only nodes;=2,plot nodes + node number)
%                   3:Element(=1,plot solid element;=2,plot only element number;=3,plot element number + element center),Gauss points(4),
%                   5: 裂缝及裂缝坐标点(=1,绘制裂缝面;=2,绘制裂缝面及坐标点;=3,绘制裂缝面及坐标点编号),
%                   6: 计算点(流体节点)及其编号
%                      =1,流体单元和节点;     =2,流体单元和节点+节点编号;    =3,流体单元和节点+单元编号;  =4,流体单元及编号+计算点及编号
%                   7: 流体单元计算点以及Gauss点相关信息
%                      =1, 节点集增强节点载荷(矢量);= 2, 计算点净水压(矢量); = 3, 计算点流量(矢量);  =4, 计算点开度(矢量);
%                      =12,流体单元净水压(云图);    =13,流体单元流量(云图);  =14,流体单元开度(云图),
%                      =21,流体单元Gauss点(单点积分)外法线向量;=22,流体单元Gauss点(单点积分)局部坐标系;=25,流体单元Gauss点(单点积分)接触力
%                   8:增强节点(=1,增强节点),9:网格线(=1),
%                   10:绘制指定单元及其节点号(起始单元号),11:绘制指定单元及其节点号(结束单元号),
%                   天然裂缝(12),单元接触状态(13),裂缝编号(14),15(Blank)
% 第2行,网格变形图: Deformation(1),Node(2),El(3),Gauss points(4)
%                   5:Crack Shape. =1,fracture area(含离散裂缝网格); =2,fracture volume; =3,仅绘制裂缝边界(不含离散裂缝网格)同时填充颜色.
%                   6:Scaling Factor.
%                   7:Boundary Conditions. =1 or =2 or =3. 
%                   8:变形前的网格. =1,变形前的全部网格；=2,变形前的网格边界
%                   9:网格线(=0,外边界线框;=1,全部网格;=2(=12,+编号),增强单元网格+外边界线框;=3,表面网格+外边界线框;=4(=14,+编号),增强单元网格+表面网格;
%                            =5(=15,+编号),裂尖增强单元(8个节点都是裂尖增强)网格+外边界线框),
%                   10:=1,裂缝边界离散点的局部坐标系;=2,裂尖增强单元的基准线(baseline);=3,裂尖增强单元的基准线(baseline)+基准局部坐标系;=4,裂尖增强节点对应的单元(参考单元);
%                      =5,裂缝边界离散点及编号
%                   11:=1,裂缝面前缘节点主应力矢量;
%                   12:=1,绘制裂尖应力计算点,CFCP=2时有效;=2,绘制计算球内的高斯点,CFCP=2时有效;=3,绘制增强节点的位移矢量,
%                   13:单元接触状态
%                   14:增强节点
%                   15:Fracture zone
% 第3行,应力云图:   Nodal Stress Contour(=1:节点值;=2:x slice;=3:y slice;=4:z slice),(2:=1,Only Mises stress;=2,仅x应力;=3,仅y应力;=4,仅z应力),
%                   主应力(3:1,主应力;2,主应力+方向;3,仅最大主应力),x or y or z slice location(4),
%                   Crack(5:=1,line;=2,shape;3,不显示增强单元),Scaling Factor(6),
%                   undeformed or Deformed(8),mesh(9),Blank(10)Blank(11),Blank(12),Blank(13),Blank(14),Slice location(15)
% 第4行,位移云图:   Nodal Displacement Contour(=1:节点值云图;=2:x slice;=3:y slice;=4:z slice),(2=1,only x;=2,only y;=3,only z),Blank(3),Blank(4),
%                   Crack(5:=1,plot crack;3,不显示增强单元),Scaling Factor(6),
%                   Blank(7),Blank(8),mesh(9),Blank(10),Blank(11),Blank(12),Blank(13),Blank(14),Slice location(15)
% 第5行,裂缝云图:   Plot Crack Contour(=1,aperture contour;=2,),
%                   2:=1,results of fluid elements and nodes; =2, results of discrete crack surface
%                   3:=1,Plot Fracture mesh
%                   4-Blank,
%                   5:编号绘制:=1,绘制裂缝面离散点编号(or流体节点编号); =2,绘制裂缝面离散单元编号(or流体单元编号); =3, nodes+elements
%                   6:Scaling Factor,7-Blank,8-Blank,
%                   9:网格线(=0,外边界线框;=1,全部网格;=2,增强单元网格+外边界线框;=3,表面网格+外边界线框;=4,增强单元网格+表面网格),Blank(10),
%                   Blank(11),12:全部天然裂缝(=1,天然裂缝;=2,天然裂缝+编号),Blank(13),Blank(14),Blank(15)
%                         1   2   3    4   5   6              7    8   9  10      11    12  13  14      15
Key_PLOT(1,:)         = [ 0,  0,  0,   0,  0,  1,             0,   1,  0  ,1     ,0     ,0  ,0  ,0  ,  8.0];  
Key_PLOT(2,:)         = [ 1,  0,  0,   0,  1,  Defor_Factor,  0,   0,  0  ,5     ,1     ,2  ,0  ,0  ,    0];  
Key_PLOT(3,:)         = [ 0,  1,  0,   0,  0,  Defor_Factor,  0,   1,  0  ,0     ,0     ,0  ,0  ,0  ,  0.5];  
Key_PLOT(4,:)         = [ 1,  3,  0,   0,  0,  Defor_Factor,  0,   0,  0  ,1     ,0     ,0  ,0  ,0  ,    5];
Key_PLOT(5,:)         = [ 1,  1,  0,   0,  0,  Defor_Factor,  0,   0,  0  ,0     ,0     ,1  ,0  ,0  ,    0];   
%##########################################################################################################
%##########################            End of user defined part        #####################################
%###########################################################################################################

PhiPsi_Post_2_Go_3D

%退出程序.
% exit
