%************************
%     初始设置
%************************
clear all; close all; clc; format compact;  format long;

%Ref: \theory_documents\042.2 Taubin_1995_Curve and surface smoothing without shrinkage.pdf

%************************
%     输入参数
%************************
%------------------------------------------
%Example 1
%------------------------------------------
% num_point = 9;
% point(1,1:3) = [ 1, 1,1];
% point(2,1:3) = [ 2, 6,0];
% point(3,1:3) = [ 3, 1,0];
% point(4,1:3) = [ 4, 1,0];
% point(5,1:3) = [ 5, 0,1];
% point(6,1:3) = [ 5.5,-1,0];
% point(7,1:3) = [ 3,-2,0];
% point(8,1:3) = [2,-2,0];
% point(9,1:3) = [1,-1,0];
% point(num_point+1,1:3) = point(1,1:3);

%------------------------------------------
%Example 2: 正方形
%------------------------------------------
% num_point = 4;
% point(1,1:3) = [ 0, 0,0];
% point(2,1:3) = [ 1, 0,0];
% point(3,1:3) = [ 1, 1,0];
% point(4,1:3) = [ 0, 1,0];
% point(num_point+1,1:3) = point(1,1:3);

%------------------------------------------
%Example 3: 标准圆
%------------------------------------------
% r= 1;
% num_point = 20;
% for i_point =1:num_point
    % c_theta = 2.0*pi/num_point*i_point;
    % point(i_point,1:3) = [r*cos(c_theta),r*sin(c_theta),0];
% end
% point(num_point+1,1:3) = point(1,1:3);

%------------------------------------------
%Example 5: 椭圆
%------------------------------------------
% num_point = 20;
% a = 2;
% b = 1;
% for i_point =1:num_point
    % c_theta = 2.0*pi/num_point*i_point;
    % point(i_point,1:3) = [a*cos(c_theta),b*sin(c_theta),0];
% end
% point(num_point+1,1:3) = point(1,1:3);

%------------------------------------------
%Example 6: 椭圆+x,y,z随机量 
%------------------------------------------
% num_point = 20;
% a = 200;
% b = 100;
% random = 30;
% for i_point =1:num_point
    % c_theta = 2.0*pi/num_point*i_point;
    % point(i_point,1:3) = [(a+random*(rand(1)*2-1))*cos(c_theta),(b+random*(rand(1)*2-1))*sin(c_theta),random*(rand(1)*2-1)];
% end
% point(num_point+1,1:3) = point(1,1:3);

%------------------------------------------
%Example 7: 椭圆+x,y,z随机量 
%------------------------------------------
% num_point = 200;
% a = 2;
% b = 1;
% random = 0.05;
% for i_point =1:num_point
    % c_theta = 2.0*pi/num_point*i_point;
    % point(i_point,1:3) = [(a+random*(rand(1)*2-1))*cos(c_theta),(b+random*(rand(1)*2-1))*sin(c_theta),random*(rand(1)*2-1)];
% end
% point(num_point+1,1:3) = point(1,1:3);

%------------------------------------------
%Example 8: 心形
%------------------------------------------
% num_point = 20;
% a = 2;
% for i_point =1:num_point
    % c_theta = 2.0*pi/num_point*i_point;
    % point(i_point,1:3) = [a*(2*cos(c_theta)-cos(2*c_theta)),a*(2*sin(c_theta)-sin(2*c_theta)),0];
% end
% point(num_point+1,1:3) = point(1,1:3);
%

%--------------------------------------------------------------------------------------------------------------------
%Example 9: 随机点：数据来自https://stackoverflow.com/questions/31464345/fitting-a-closed-curve-to-a-set-of-points
%--------------------------------------------------------------------------------------------------------------------
num_point = 47;
point(1,1:3) = [ 6.55525 ,  3.05472,0];
point(2,1:3) = [ 6.17284 ,  2.802609,0];
point(3,1:3) = [ 5.53946 ,  2.649209,0];
point(4,1:3) = [ 4.93053 ,  2.444444,0];
point(5,1:3) = [ 4.32544 ,  2.318749,0];
point(6,1:3) = [ 3.90982 ,  2.2875,0  ];
point(7,1:3) =    [ 3.51294 ,  2.221875,0];
point(8,1:3) =    [ 3.09107 ,  2.29375,0 ];
point(9,1:3) =    [ 2.64013 ,  2.4375,0  ];
point(10,1:3) =    [ 2.275444,  2.653124,0];
point(11,1:3) =    [ 2.137945,  3.26562,0 ];
point(12,1:3) =    [ 2.15982 ,  3.84375,0 ];
point(13,1:3) =    [ 2.20982 ,  4.31562,0 ];
point(14,1:3) =    [ 2.334704,  4.87873,0 ];
point(15,1:3) =    [ 2.314264,  5.5047,0  ];
point(16,1:3) =    [ 2.311709,  5.9135,0  ];
point(17,1:3) =    [ 2.29638 ,  6.42961,0 ];
point(18,1:3) =    [ 2.619374,  6.75021,0 ];
point(19,1:3) =    [ 3.32448 ,  6.66353,0 ];
point(20,1:3) =    [ 3.31582 ,  5.68866,0 ];
point(21,1:3) =    [ 3.35159 ,  5.17255,0 ];
point(22,1:3) =    [ 3.48482 ,  4.73125,0 ];
point(23,1:3) =    [ 3.70669 ,  4.51875,0 ];
point(24,1:3) =    [ 4.23639 ,  4.58968,0 ];
point(25,1:3) =    [ 4.39592 ,  4.94615,0 ];
point(26,1:3) =    [ 4.33527 ,  5.33862,0 ];
point(27,1:3) =    [ 3.95968 ,  5.61967,0 ];
point(28,1:3) =    [ 3.56366 ,  5.73976,0 ];
point(29,1:3) =    [ 3.78818 ,  6.55292,0 ];
point(30,1:3) =    [ 4.27712 ,  6.8283,0  ];
point(31,1:3) =    [ 4.89532 ,  6.78615,0 ];
point(32,1:3) =    [ 5.35334 ,  6.72433,0 ];
point(33,1:3) =    [ 5.71583 ,  6.54449,0 ];
point(34,1:3) =    [ 6.13452 ,  6.46019,0 ];
point(35,1:3) =    [ 6.54478 ,  6.26068,0 ];
point(36,1:3) =    [ 6.7873  ,  5.74615,0 ];
point(37,1:3) =    [ 6.64086 ,  5.25269,0 ];
point(38,1:3) =    [ 6.45649 ,  4.86206,0 ];
point(39,1:3) =    [ 5.41586 ,  6.46519,0 ];
point(40,1:3) =    [ 5.44711 ,  4.26519,0 ];
point(41,1:3) =    [ 5.04087 ,  4.10581,0 ];
point(42,1:3) =    [ 4.70013 ,  3.67405,0 ];
point(43,1:3) =    [ 4.83482 ,  3.4375,0  ];
point(44,1:3) =    [ 5.34086 ,  3.43394,0 ];
point(45,1:3) =    [ 5.76392 ,  3.55156,0 ];
point(46,1:3) =    [ 6.37056 ,  3.8778,0  ];
point(47,1:3) =    [ 6.53116 ,  3.47228,0 ];
point(num_point+1,1:3) = point(1,1:3);


%************************
%    控制参数
%************************
% n_Loops =  1;
n_Loops =  20;

% lambda =  0.5;   %lambda的绝对值要小于mu
% mu     =  0.0;

% lambda =  0.3;   %lambda的绝对值要小于mu
% mu     = -0.35;
% mu     = 0.0;

% lambda =  0.5;   %lambda的绝对值要小于mu
% mu     = -0.53;

% 默认1
% lambda =  0.33;   %lambda的绝对值要小于mu
% mu     = -0.34;

% 默认2
lambda =  0.330;   %lambda的绝对值要小于mu
mu     = -0.331;

%************************
%    控制参数相关变量
%************************
k_pb=0.1;
mu_calculated = 1./(k_pb-(1/lambda))
k=2;
f=(1-k.*lambda).*(1-mu.*k);

%************************
%     绘图设置
%************************
c_figure = figure('units','normalized','position',[0.1,0.1,0.8,0.8],'Visible','on','Renderer','OpenGL');   %采用OpenGL渲染器   
hold on;
axis off; 
axis equal;
% view(17,29)  %视角
offset = 0.01;
for i_point =1:num_point
	plot3(point(i_point,1),point(i_point,2),point(i_point,3),'k.','MarkerSize',20.0)  
	ts = text(point(i_point,1)+offset,point(i_point,2)+offset,point(i_point,3)+offset,num2str(i_point),'Color','blue','FontSize',15,'FontName','Consolas','FontAngle','italic');
end
plot3(point(:,1)',point(:,2)',point(:,3)','LineWidth',2.5,'Color','red')	


% old_point = point;
% Center(1) = sum(old_point(1:num_point,1))/num_point;
% Center(2) = sum(old_point(1:num_point,2))/num_point;
% Center(3) = sum(old_point(1:num_point,3))/num_point;
	
%************************
%        循环
%************************
myMap = rand(n_Loops, 3);
for i_Loop=1:n_Loops

    disp(['Loop_counter: ',num2str(i_Loop)])
	 
    % Step 1
    for i_point=1:num_point
	    c_point = point(i_point,1:3);
		nex_point = point(i_point+1,1:3);
	    if i_point>=2
			pre_point = point(i_point-1,1:3);
		else
			pre_point = point(num_point,1:3);
		end
		
		%%%%%%%%%%%%%%%%%%%%%%%%%
		%   OPTTION 1: 普通
		%%%%%%%%%%%%%%%%%%%%%%%%%
        % L1_vector(i_point,1:3) = 0.5*(nex_point - c_point) + 0.5*(pre_point - c_point);

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%   OPTTION 2：考虑权重修改L算子
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		pre_L = sqrt((pre_point(1)-c_point(1))^2+(pre_point(2)-c_point(2))^2+(pre_point(3)-c_point(3))^2);
		nex_L = sqrt((nex_point(1)-c_point(1))^2+(nex_point(2)-c_point(2))^2+(nex_point(3)-c_point(3))^2);
		weight_pre = 1/pre_L;
		weight_nex = 1/nex_L;
		L1_vector(i_point,1:3) = (weight_pre*pre_point+weight_nex*nex_point)/(weight_pre+weight_nex)- c_point;
		
		
    end
	
	% point(:,1:3) = point(:,1:3)+lambda*L_vector(:,1:3);

	for i_point=1:num_point
	    point(i_point,1:3) = point(i_point,1:3)+lambda*L1_vector(i_point,1:3);
	end
	point(num_point+1,1:3) = point(1,1:3);
	
	% Step 2
    for i_point=1:num_point
	    c_point = point(i_point,1:3);
		nex_point = point(i_point+1,1:3);
	    if i_point>=2
			pre_point = point(i_point-1,1:3);
		else
			pre_point = point(num_point,1:3);
		end
		
		%%%%%%%%%%%%%%%%%%%%%%%%%
		%   OPTTION 1: 普通
		%%%%%%%%%%%%%%%%%%%%%%%%%
        % L2_vector(i_point,1:3) = 0.5*(nex_point - c_point) + 0.5*(pre_point - c_point);

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%   OPTTION 2：考虑权重修改L算子
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		pre_L = sqrt((pre_point(1)-c_point(1))^2+(pre_point(2)-c_point(2))^2+(pre_point(3)-c_point(3))^2);
		nex_L = sqrt((nex_point(1)-c_point(1))^2+(nex_point(2)-c_point(2))^2+(nex_point(3)-c_point(3))^2);
		weight_pre = 1/pre_L;
		weight_nex = 1/nex_L;
		L2_vector(i_point,1:3) = (weight_pre*pre_point+weight_nex*nex_point)/(weight_pre+weight_nex)- c_point;
		
    end
	
	for i_point=1:num_point
	    point(i_point,1:3) = point(i_point,1:3)+mu*L2_vector(i_point,1:3);
	end
	point(num_point+1,1:3) = point(1,1:3);
	
	
	%绘图:绘制循环过程中的结果(每5步绘制一次)
	if mod(i_Loop,5)==0
		for i_point =1:num_point
			plot3(point(i_point,1),point(i_point,2),point(i_point,3),'k.','MarkerSize',10.0)  
			% ts = text(point(i_point,1)+offset,point(i_point,2)+offset,point(i_point,3)+offset,num2str(i_point),'Color','blue','FontSize',15,'FontName','Consolas','FontAngle','italic');
		end
		point(num_point+1,1:3) = point(1,1:3);
		% plot3(point(:,1)',point(:,2)',point(:,3)','LineWidth',0.1,'Color',myMap(i_Loop,:),'linestyle','--')	
		plot3(point(:,1)',point(:,2)',point(:,3)','LineWidth',0.1,'Color',[0.5,0.5,0.5],'linestyle','--')	
	end
	
	%最后一步的绘图
	if i_Loop==n_Loops
		for i_point =1:num_point
			plot3(point(i_point,1),point(i_point,2),point(i_point,3),'k.','MarkerSize',20.0)  
			ts = text(point(i_point,1)+offset,point(i_point,2)+offset,point(i_point,3)+offset,num2str(i_point),'Color','blue','FontSize',15,'FontName','Consolas','FontAngle','italic');
		end
		point(num_point+1,1:3) = point(1,1:3);
		plot3(point(:,1)',point(:,2)',point(:,3)','LineWidth',3,'Color',[127/255,255/255,0/255])	
	end
end

title(['Smoothed curve after ',num2str(n_Loops), ' loops.'],'FontSize',16) 		
	
fcw(gca);



