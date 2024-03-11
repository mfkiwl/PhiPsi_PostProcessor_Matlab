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

function Plot_Node_Disp_3D_Slice(isub,Crack_X,Crack_Y,Crack_Z,Post_Enriched_Nodes,POS)
% This function plots the 3D slice of nodal dispalcement contour.

global Node_Coor Elem_Node Key_POST_HF
global Num_Node Num_Elem
global Min_X_Coor Max_X_Coor Min_Y_Coor Max_Y_Coor Min_Z_Coor Max_Z_Coor
global Key_PLOT aveg_area_ele
global Size_Font Elem_Fontcolor Elem_Fontsize Node_Fontcolor Node_Fontsize
global num_Crack num_of_Material
global Color_Crack Width_Crack Full_Pathname
global Color_Backgro_Mesh_1 Color_Backgro_Mesh_2 Color_Backgro_Mesh_3 Color_Backgro_Mesh_4
global Color_Backgro_Mesh_5 Color_Backgro_Mesh_6 Color_Backgro_Mesh_7
global Color_Backgro_Mesh_8 Color_Backgro_Mesh_9 Color_Backgro_Mesh_10
global Elem_Material Num_Step_to_Plot DISP Num_Foc_z
global Num_Foc_x Num_Foc_y Foc_x Foc_y Num_Foc_z Foc_z FORCE_Matrix
global Key_Data_Format
global Crack_node_X Crack_node_Y Crack_node_Z
global Crack_Ele_1 Crack_Ele_2 Crack_Ele_3
global Model_Outline Model_OutArea
global Num_Bou_x Num_Bou_y Num_Bou_z Bou_x Bou_y Bou_z
global Crack_node_local_X Crack_node_local_Y Crack_node_local_Z
global Crack_Node_in_Ele
global Crack3D_Vector_S1_X Crack3D_Vector_S1_Y Crack3D_Vector_S1_Z
global Crack3D_CalP_X Crack3D_CalP_Y Crack3D_CalP_Z
global Crack3D_Vertex_Vector_X_X 
global Crack3D_Vertex_Vector_X_Y  
global Crack3D_Vertex_Vector_X_Z 
global Crack3D_Vertex_Vector_Y_X  
global Crack3D_Vertex_Vector_Y_Y  
global Crack3D_Vertex_Vector_Y_Z  
global Crack3D_Vertex_Vector_Z_X  
global Crack3D_Vertex_Vector_Z_Y  
global Crack3D_Vertex_Vector_Z_Z 	
global Crack3D_Meshed_Outline
global Num_Accuracy_Contour Key_Contour_Metd
global volume_ele max_volume_ele min_volume_ele aveg_volume_ele 
global Num_Contourlevel Key_Flipped_Gray Color_Contourlevel
global Color_Mesh
global Title_Font Key_Figure_Control_Widget
global Ele_Permeability_3D

disp('    > Plotting slice of nodal dispalcement contour....') 

% 绘图算法设置为1
if Key_Contour_Metd==2
    disp('      *************************************************************') 
	disp('      WARNING :: Key_Contour_Metd is set to 1 for 3D slice ploting!') 
	disp('                 Key_Contour_Metd = 1                              ') 
    disp('      *************************************************************') 
    Key_Contour_Metd = 1;
end

scale         = Key_PLOT(4, 6);
Slice_Type    = Key_PLOT(4, 1);
Location_Plot = Key_PLOT(4,15);

% Get the new coordinates of all nodes.
New_Node_Coor(1:Num_Node,1) = Node_Coor(1:Num_Node,1) + scale*DISP(1:Num_Node,2);
New_Node_Coor(1:Num_Node,2) = Node_Coor(1:Num_Node,2) + scale*DISP(1:Num_Node,3);
New_Node_Coor(1:Num_Node,3) = Node_Coor(1:Num_Node,3) + scale*DISP(1:Num_Node,4);


% Get the maximum and minimum value of the new coordinates of all nodes.
Min_X_Coor_New = min(min(New_Node_Coor(1:Num_Node,1)));
Max_X_Coor_New = max(max(New_Node_Coor(1:Num_Node,1)));
Min_Y_Coor_New = min(min(New_Node_Coor(1:Num_Node,2)));
Max_Y_Coor_New = max(max(New_Node_Coor(1:Num_Node,2)));
Min_Z_Coor_New = min(min(New_Node_Coor(1:Num_Node,3)));
Max_Z_Coor_New = max(max(New_Node_Coor(1:Num_Node,3)));


% Get resample coors.
delta = (aveg_volume_ele^(1/3))/Num_Accuracy_Contour;
gx    = Min_X_Coor:delta:Max_X_Coor; 
gy    = Min_Y_Coor:delta:Max_Y_Coor;
gz    = Min_Z_Coor:delta:Max_Z_Coor;

% 绘图设置.
Start_Dis_Type =1;     
End_Dis_Type   =2; 

% 根据切片类型不同进行设置待绘制的坐标
if Slice_Type==2   %x轴切片
   g1 =gy;
   g2 =gz;
   Ploted_Nodes = find(Node_Coor(1:Num_Node,1)==Location_Plot);%根据坐标找到节点号
   num_of_Ploted_Nodes = size(Ploted_Nodes,1);
   Ploted_Node_Coor_1(1:num_of_Ploted_Nodes) = New_Node_Coor(Ploted_Nodes,2);
   Ploted_Node_Coor_2(1:num_of_Ploted_Nodes) = New_Node_Coor(Ploted_Nodes,3);
   DISP_1 = DISP(Ploted_Nodes,3);
   DISP_2 = DISP(Ploted_Nodes,4);   
elseif Slice_Type==3   %y轴切片
   g1 =gx;
   g2 =gz;
   Ploted_Nodes = find(Node_Coor(1:Num_Node,2)==Location_Plot);%根据坐标找到节点号
   num_of_Ploted_Nodes = size(Ploted_Nodes,1);
   Ploted_Node_Coor_1(1:num_of_Ploted_Nodes) = New_Node_Coor(Ploted_Nodes,1);
   Ploted_Node_Coor_2(1:num_of_Ploted_Nodes) = New_Node_Coor(Ploted_Nodes,3);
   DISP_1 = DISP(Ploted_Nodes,2);
   DISP_2 = DISP(Ploted_Nodes,4);
elseif Slice_Type==4   %z轴切片
   g1 =gx;
   g2 =gy;
   Ploted_Nodes = find(Node_Coor(1:Num_Node,3)==Location_Plot);%根据坐标找到节点号
   num_of_Ploted_Nodes = size(Ploted_Nodes,1);
   Ploted_Node_Coor_1(1:num_of_Ploted_Nodes) = New_Node_Coor(Ploted_Nodes,1);
   Ploted_Node_Coor_2(1:num_of_Ploted_Nodes) = New_Node_Coor(Ploted_Nodes,2);
   DISP_1 = DISP(Ploted_Nodes,2);
   DISP_2 = DISP(Ploted_Nodes,3);   
   % DISP_1 = DISP(Ploted_Nodes,2)*1000.0;
   % DISP_2 = DISP(Ploted_Nodes,3)*1000.0;      
end
% Plot the contours.
for i_Plot = Start_Dis_Type:End_Dis_Type
    % New figure.
    Tools_New_Figure
	hold on;
	if i_Plot == 1
		if  Slice_Type==2       %x轴切片
		    disp('      Resample Uy....')  
		elseif  Slice_Type==3   %y轴切片
		    disp('      Resample Ux....') 
	    elseif Slice_Type==4    %z轴切片
		    disp('      Resample Ux....') 			
		end
		if 	Key_Contour_Metd==1
			[X,Y,Ux] = griddata(Ploted_Node_Coor_1,Ploted_Node_Coor_2,DISP_1,...
				       unique(Ploted_Node_Coor_1),unique(Ploted_Node_Coor_2)');
		elseif 	Key_Contour_Metd==2
			[Ux,X,Y] = Tools_gridfit(Ploted_Node_Coor_1,Ploted_Node_Coor_2,DISP_1,g1,g2);
		end	
		
		if  Slice_Type==2
		    disp('      Contouring Uy....')  
		elseif  Slice_Type==3
		    disp('      Contouring Ux....') 
	    elseif Slice_Type==4 
		    disp('      Contouring Ux....') 			
		end	
		
		contourf(X,Y,Ux,Num_Contourlevel,'LineStyle','none')
		
		
		if  Slice_Type==2
		    title(['Displacement y of x-axis slice (x= ',num2str(Location_Plot), ' m)'],'FontName',Title_Font,'FontSize',Size_Font) 
		elseif  Slice_Type==3
		    title(['Displacement x of y-axis slice (y= ',num2str(Location_Plot), ' m)'],'FontName',Title_Font,'FontSize',Size_Font)
	    elseif Slice_Type==4 
		    title(['Displacement x of z-axis slice (z= ',num2str(Location_Plot), ' m)'],'FontName',Title_Font,'FontSize',Size_Font)			
		end
		
		clear Ux
	elseif i_Plot == 2
		if  Slice_Type==2
		    disp('      Resample Uz....')  
		elseif  Slice_Type==3
		    disp('      Resample Uz....') 
	    elseif Slice_Type==4 
		    disp('      Resample Uy....') 			
		end
		if 	Key_Contour_Metd==1
	    	[X,Y,Uy] = griddata(Ploted_Node_Coor_1,Ploted_Node_Coor_2,DISP_2,...
				       unique(Ploted_Node_Coor_1),unique(Ploted_Node_Coor_2)');	
		elseif 	Key_Contour_Metd==2
			[Uy,X,Y] = Tools_gridfit(Ploted_Node_Coor_1,Ploted_Node_Coor_2,DISP_2,g1,g2);
		end			
		
		if  Slice_Type==2
		    disp('      Contouring Uz....')  
		elseif  Slice_Type==3
		    disp('      Contouring Uz....')  
	    elseif Slice_Type==4 
		    disp('      Contouring Uy....') 			
		end	 
		
		contourf(X,Y,Uy,Num_Contourlevel,'LineStyle','none')
		

		if  Slice_Type==2
		    title(['Displacement z of x-axis slice (x= ',num2str(Location_Plot), ' m)'],'FontName',Title_Font,'FontSize',Size_Font) 
		elseif  Slice_Type==3
		    title(['Displacement z of y-axis slice (y= ',num2str(Location_Plot), ' m)'],'FontName',Title_Font,'FontSize',Size_Font)
	    elseif Slice_Type==4 
		    title(['Displacement y of z-axis slice (z= ',num2str(Location_Plot), ' m)'],'FontName',Title_Font,'FontSize',Size_Font)			
		end
		clear Uy
	end
	
	% Set colormap.
	if Key_Flipped_Gray==0
		colormap(Color_Contourlevel)
	elseif Key_Flipped_Gray==1
		colormap(flipud(gray))
	end
	
    % Plot the deformed mesh.
	if Key_PLOT(4,9)==1
	    disp('      Plot deformed mesh....') 	
        for iElem = 1:Num_Elem
			NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
				  Elem_Node(iElem,3) Elem_Node(iElem,4) ...
				  Elem_Node(iElem,5) Elem_Node(iElem,6) ...
				  Elem_Node(iElem,7) Elem_Node(iElem,8)];                             % Nodes for current element
				  
			%6个面循环
			for i_face =1:6
			    if i_face==1
					NN_4 = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
							Elem_Node(iElem,3) Elem_Node(iElem,4) Elem_Node(iElem,1)];    
				elseif i_face==2
					NN_4 = [Elem_Node(iElem,5) Elem_Node(iElem,6) ...
							Elem_Node(iElem,7) Elem_Node(iElem,8) Elem_Node(iElem,5)];     
				elseif i_face==3
					NN_4 = [Elem_Node(iElem,1) Elem_Node(iElem,2) ...
							Elem_Node(iElem,6) Elem_Node(iElem,5) Elem_Node(iElem,1)];  
				elseif i_face==4
					NN_4 = [Elem_Node(iElem,2) Elem_Node(iElem,6) ...
							Elem_Node(iElem,7) Elem_Node(iElem,3) Elem_Node(iElem,2)];  
				elseif i_face==5
					NN_4 = [Elem_Node(iElem,3) Elem_Node(iElem,4) ...
							Elem_Node(iElem,8) Elem_Node(iElem,7) Elem_Node(iElem,3)];     
				elseif i_face==6
					NN_4 = [Elem_Node(iElem,1) Elem_Node(iElem,4) ...
							Elem_Node(iElem,8) Elem_Node(iElem,5) Elem_Node(iElem,1)];   
				end
				
				%假如当前面的4个节点在Slice上		
				if ismember(NN_4,Ploted_Nodes)
					if Slice_Type==2       %x轴切片
						xi = New_Node_Coor(NN_4',2);                                         
						yi = New_Node_Coor(NN_4',3);      
					elseif Slice_Type==3   %y轴切片
						xi = New_Node_Coor(NN_4',1);                                         
						yi = New_Node_Coor(NN_4',3);   
					elseif Slice_Type==4   %z轴切片
						xi = New_Node_Coor(NN_4',1);                                         
						yi = New_Node_Coor(NN_4',2);   
					end                                           % Deformed y-coordinates of nodes
					plot(xi,yi,'LineWidth',0.5,'Color',Color_Mesh)
				end	
            end			
		end 	
	end
    % Plot crack.
	if Key_PLOT(4,5)==1
	  file_type = ['_fracture_coordinates.txt'];
	  String = [Full_Pathname,file_type];	
	  if exist(String, 'file') ==2 
		delete(String)
	  end
	  if isempty(Crack_X)==0
	    disp('      Plot the crack....') 	
		    for i_Crack = 1:num_Crack(isub)
			    nnode = size(Crack_node_X{i_Crack},2);
	            % for i_nnode=1:nnode
					% c_node_x = [Crack_node_X{i_Crack}(i_nnode)];
					% c_node_y = [Crack_node_Y{i_Crack}(i_nnode)];
					% c_node_z = [Crack_node_Z{i_Crack}(i_nnode)];
				% end
				
				nele = size(Crack_Ele_1{i_Crack},2);
				for i_nele=1:nele
					cr_node_1 = Crack_Ele_1{i_Crack}(i_nele);
					cr_node_2 = Crack_Ele_2{i_Crack}(i_nele);
					cr_node_3 = Crack_Ele_3{i_Crack}(i_nele);
				    cr_el_line1_P1 = [Crack_node_X{i_Crack}(cr_node_1),Crack_node_Y{i_Crack}(cr_node_1),Crack_node_Z{i_Crack}(cr_node_1)];
					cr_el_line1_P2 = [Crack_node_X{i_Crack}(cr_node_2),Crack_node_Y{i_Crack}(cr_node_2),Crack_node_Z{i_Crack}(cr_node_2)];
				    cr_el_line2_P1 = [Crack_node_X{i_Crack}(cr_node_2),Crack_node_Y{i_Crack}(cr_node_2),Crack_node_Z{i_Crack}(cr_node_2)];
					cr_el_line2_P2 = [Crack_node_X{i_Crack}(cr_node_3),Crack_node_Y{i_Crack}(cr_node_3),Crack_node_Z{i_Crack}(cr_node_3)];
				    cr_el_line3_P1 = [Crack_node_X{i_Crack}(cr_node_3),Crack_node_Y{i_Crack}(cr_node_3),Crack_node_Z{i_Crack}(cr_node_3)];
					cr_el_line3_P2 = [Crack_node_X{i_Crack}(cr_node_1),Crack_node_Y{i_Crack}(cr_node_1),Crack_node_Z{i_Crack}(cr_node_1)];		
					%ooooooooooooooooooooooooooooooooooooooo
					%                x轴切片
					%ooooooooooooooooooooooooooooooooooooooo
					if  Slice_Type==2 
					    num_intersection = 0;
						%裂缝面三角形单元的第1条边
						if  (Location_Plot >=cr_el_line1_P1(1) & Location_Plot <= cr_el_line1_P2(1)) || (Location_Plot >=cr_el_line1_P2(1) & Location_Plot <= cr_el_line1_P1(1))
							num_intersection = num_intersection +1;
							%........................................................................
							%计算交点坐标
							%利用空间直线的两点式方程：(x-x1)/(x2-x1)=(y-y1)/(y2-y1)=(z-z1)/(z2-z1)
							%........................................................................
							x1 = cr_el_line1_P1(1); y1 = cr_el_line1_P1(2); z1 = cr_el_line1_P1(3);
							x2 = cr_el_line1_P2(1); y2 = cr_el_line1_P2(2); z2 = cr_el_line1_P2(3);
							inter_y = (Location_Plot-x1)/(x2-x1)*(y2-y1)+y1;
							inter_z = (Location_Plot-x1)/(x2-x1)*(z2-z1)+z1;
							intersection_Points(num_intersection,1:2) = [inter_y,inter_z];
							
						end		
						%裂缝面三角形单元的第2条边
						if  (Location_Plot >=cr_el_line2_P1(1) & Location_Plot <= cr_el_line2_P2(1)) || (Location_Plot >=cr_el_line2_P2(1) & Location_Plot <= cr_el_line2_P1(1))
							num_intersection = num_intersection +1;
							x1 = cr_el_line2_P1(1); y1 = cr_el_line2_P1(2); z1 = cr_el_line2_P1(3);
							x2 = cr_el_line2_P2(1); y2 = cr_el_line2_P2(2); z2 = cr_el_line2_P2(3);
							inter_y = (Location_Plot-x1)/(x2-x1)*(y2-y1)+y1;
							inter_z = (Location_Plot-x1)/(x2-x1)*(z2-z1)+z1;
							intersection_Points(num_intersection,1:2) = [inter_y,inter_z];
						end	
						%裂缝面三角形单元的第3条边
						if  (Location_Plot >=cr_el_line3_P1(1) & Location_Plot <= cr_el_line3_P2(1)) || (Location_Plot >=cr_el_line3_P2(1) & Location_Plot <= cr_el_line3_P1(1))
							num_intersection = num_intersection +1;
							x1 = cr_el_line3_P1(1); y1 = cr_el_line3_P1(2); z1 = cr_el_line3_P1(3);
							x2 = cr_el_line3_P2(1); y2 = cr_el_line3_P2(2); z2 = cr_el_line3_P2(3);
							inter_y = (Location_Plot-x1)/(x2-x1)*(y2-y1)+y1;
							inter_z = (Location_Plot-x1)/(x2-x1)*(z2-z1)+z1;
							intersection_Points(num_intersection,1:2) = [inter_y,inter_z];
						end	
					%ooooooooooooooooooooooooooooooooooooooo
					%                y轴切片
					%ooooooooooooooooooooooooooooooooooooooo
					elseif  Slice_Type==3 %y轴切片
					    num_intersection = 0;
						%裂缝面三角形单元的第1条边
						if  (Location_Plot >=cr_el_line1_P1(2) & Location_Plot <= cr_el_line1_P2(2)) || (Location_Plot >=cr_el_line1_P2(2) & Location_Plot <= cr_el_line1_P1(2))
							num_intersection = num_intersection +1;
							%........................................................................
							%计算交点坐标
							%利用空间直线的两点式方程：(x-x1)/(x2-x1)=(y-y1)/(y2-y1)=(z-z1)/(z2-z1)
							%........................................................................
							x1 = cr_el_line1_P1(1); y1 = cr_el_line1_P1(2); z1 = cr_el_line1_P1(3);
							x2 = cr_el_line1_P2(1); y2 = cr_el_line1_P2(2); z2 = cr_el_line1_P2(3);
							inter_x = (Location_Plot-y1)/(y2-y1)*(x2-x1)+x1;
							inter_z = (Location_Plot-y1)/(y2-y1)*(z2-z1)+z1;
							intersection_Points(num_intersection,1:2) = [inter_x,inter_z];
						end		
						%裂缝面三角形单元的第2条边
						if  (Location_Plot >=cr_el_line2_P1(2) & Location_Plot <= cr_el_line2_P2(2)) || (Location_Plot >=cr_el_line2_P2(2) & Location_Plot <= cr_el_line2_P1(2))
							num_intersection = num_intersection +1;
							x1 = cr_el_line2_P1(1); y1 = cr_el_line2_P1(2); z1 = cr_el_line2_P1(3);
							x2 = cr_el_line2_P2(1); y2 = cr_el_line2_P2(2); z2 = cr_el_line2_P2(3);
							inter_x = (Location_Plot-y1)/(y2-y1)*(x2-x1)+x1;
							inter_z = (Location_Plot-y1)/(y2-y1)*(z2-z1)+z1;
							intersection_Points(num_intersection,1:2) = [inter_x,inter_z];					
						end	
						%裂缝面三角形单元的第3条边
						if  (Location_Plot >=cr_el_line3_P1(2) & Location_Plot <= cr_el_line3_P2(2)) || (Location_Plot >=cr_el_line3_P2(2) & Location_Plot <= cr_el_line3_P1(2))
							num_intersection = num_intersection +1;
							x1 = cr_el_line3_P1(1); y1 = cr_el_line3_P1(2); z1 = cr_el_line3_P1(3);
							x2 = cr_el_line3_P2(1); y2 = cr_el_line3_P2(2); z2 = cr_el_line3_P2(3);
							inter_x = (Location_Plot-y1)/(y2-y1)*(x2-x1)+x1;
							inter_z = (Location_Plot-y1)/(y2-y1)*(z2-z1)+z1;
							intersection_Points(num_intersection,1:2) = [inter_x,inter_z];					
						end	
					%ooooooooooooooooooooooooooooooooooooooo
					%                z轴切片
					%ooooooooooooooooooooooooooooooooooooooo						
					elseif Slice_Type==4  %z轴切片
					    num_intersection = 0;
						%裂缝面三角形单元的第1条边
						if  (Location_Plot >=cr_el_line1_P1(3) & Location_Plot <= cr_el_line1_P2(3)) || (Location_Plot >=cr_el_line1_P2(3) & Location_Plot <= cr_el_line1_P1(3))
							num_intersection = num_intersection +1;
							%........................................................................
							%计算交点坐标
							%利用空间直线的两点式方程：(x-x1)/(x2-x1)=(y-y1)/(y2-y1)=(z-z1)/(z2-z1)
							%........................................................................
							x1 = cr_el_line1_P1(1); y1 = cr_el_line1_P1(2); z1 = cr_el_line1_P1(3);
							x2 = cr_el_line1_P2(1); y2 = cr_el_line1_P2(2); z2 = cr_el_line1_P2(3);
							inter_x = (Location_Plot-z1)/(z2-z1)*(x2-x1)+x1;
							inter_y = (Location_Plot-z1)/(z2-z1)*(y2-y1)+y1;
							intersection_Points(num_intersection,1:2) = [inter_x,inter_y];
							
						end		
						%裂缝面三角形单元的第2条边
						if  (Location_Plot >=cr_el_line2_P1(3) & Location_Plot <= cr_el_line2_P2(3)) || (Location_Plot >=cr_el_line2_P2(3) & Location_Plot <= cr_el_line2_P1(3))
							num_intersection = num_intersection +1;
							x1 = cr_el_line2_P1(1); y1 = cr_el_line2_P1(2); z1 = cr_el_line2_P1(3);
							x2 = cr_el_line2_P2(1); y2 = cr_el_line2_P2(2); z2 = cr_el_line2_P2(3);
							inter_x = (Location_Plot-z1)/(z2-z1)*(x2-x1)+x1;
							inter_y = (Location_Plot-z1)/(z2-z1)*(y2-y1)+y1;
							intersection_Points(num_intersection,1:2) = [inter_x,inter_y];
						end	
						%裂缝面三角形单元的第3条边
						if  (Location_Plot >=cr_el_line3_P1(3) & Location_Plot <= cr_el_line3_P2(3)) || (Location_Plot >=cr_el_line3_P2(3) & Location_Plot <= cr_el_line3_P1(3))
							num_intersection = num_intersection +1;
							x1 = cr_el_line3_P1(1); y1 = cr_el_line3_P1(2); z1 = cr_el_line3_P1(3);
							x2 = cr_el_line3_P2(1); y2 = cr_el_line3_P2(2); z2 = cr_el_line3_P2(3);
							inter_x = (Location_Plot-z1)/(z2-z1)*(x2-x1)+x1;
							inter_y = (Location_Plot-z1)/(z2-z1)*(y2-y1)+y1;
							intersection_Points(num_intersection,1:2) = [inter_x,inter_y];
						end						
					end
					%%/////////////////////////////////////////////////////
					%  如果有>两个交点，绘制交点(暂时不支持变形后的裂缝)
					%//////////////////////////////////////////////////////
					if num_intersection>=2
						% plot(intersection_Points(1:2,1),intersection_Points(1:2,2),'LineWidth',2.0,'Color','black')
						plot(intersection_Points(1:num_intersection,1),intersection_Points(1:num_intersection,2),'LineWidth',2.0,'Color','black') % 2024-02-28.
						%保存裂缝坐标
						if intersection_Points(1:2,2)>0
							c_File = fopen(String, 'a+');
							fprintf(c_File,'%10e          %10e\n',intersection_Points(1,1:2));
							fprintf(c_File,'%10e          %10e\n',intersection_Points(1,1),-intersection_Points(1,2));
							% fprintf(c_File,'%10e          %10e\n',intersection_Points(2,1:2));
							fclose(c_File);		
                        end						
					end						
				end					
			end
	  end
    end
	
	% The range of the plot.
	min_x = Min_X_Coor_New; max_x = Max_X_Coor_New;
    min_y = Min_Y_Coor_New; max_y = Max_Y_Coor_New;	
	min_z = Min_Z_Coor_New; max_z = Max_Z_Coor_New;
	
	if Slice_Type==2       %x轴切片
		axis([min_y max_y min_z max_z]);	      
	elseif Slice_Type==3   %y轴切片
        axis([min_x max_x min_z max_z]);	 
	elseif Slice_Type==4   %z轴切片
		axis([min_x max_x min_y max_y]);	
	end     	
	
	axis equal; 

	
    set(gca,'XTick',[],'YTick',[],'XColor','w','YColor','w')

	%colorbar('FontAngle','italic','FontName',Title_Font,'FontSize',Size_Font);
	
	cbar = colorbar;
	set(get(cbar,'title'),'string','m');
	% set(get(cbar,'title'),'string','mm');
	
	% colorbar('FontName',Title_Font,'FontSize',Size_Font);
	
	% 2024-02-25.
    if i_Plot == 1	   
	  % clim(1)=min(DISP_1)*1000.0;
	  % clim(2)=max(DISP_1)*1000.0;
	  clim(1)=min(DISP_1);
	  clim(2)=max(DISP_1);	  
    elseif i_Plot == 2
	  % clim(1)=min(DISP_2)*1000.0;
	  % clim(2)=max(DISP_2)*1000.0;
	  clim(1)=min(DISP_2);
	  clim(2)=max(DISP_2);	  
    end
	% clim
	ylim(cbar,[clim(1) clim(2)]);
	
	numpts = 10 ;    % Number of points to be displayed on colorbar
	kssv = linspace(clim(1),clim(2),numpts);
	set(cbar,'YtickMode','manual','YTick',kssv); % Set the tickmode to manual
	for i = 1:numpts
		imep = num2str(kssv(i),'%+3.2E');
		vasu(i) = {imep} ;
	end
	set(cbar,'YTickLabel',vasu(1:numpts),'FontName',Title_Font,'FontSize',Size_Font);

	% Active Figure control widget (2021-08-01)
	% Ref: https://ww2.mathworks.cn/matlabcentral/fileexchange/38019-figure-control-widget
	% Press q to exit.
	% Press r (or double-click) to reset to the initial.
	if Key_Figure_Control_Widget==1
		fcw(gca);
	end


	% Save pictures.	
	if  Slice_Type==2
		if i_Plot == 1
			Save_Picture(c_figure,Full_Pathname,'dpyn')
		elseif i_Plot ==2
			Save_Picture(c_figure,Full_Pathname,'dpzn')
		end
	elseif  Slice_Type==3
		if i_Plot == 1
			Save_Picture(c_figure,Full_Pathname,'dpxn')
		elseif i_Plot ==2
			Save_Picture(c_figure,Full_Pathname,'dpzn')
		end
	elseif Slice_Type==4 
		if i_Plot == 1
			Save_Picture(c_figure,Full_Pathname,'dpxn')
		elseif i_Plot ==2
			Save_Picture(c_figure,Full_Pathname,'dpyn')
		end		
	end	

end


