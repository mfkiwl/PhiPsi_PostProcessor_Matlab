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

function Plot_3D_curves(POST_Substep)
% ����3D����.
% Firstly written on 2022-04-24.

global Key_PLOT Full_Pathname Num_Node Num_Foc_x Num_Foc_y Foc_x Foc_y
global num_Crack Key_Dynamic Real_Iteras Real_Sub Key_Contour_Metd
global Output_Freq num_Output_Sub Key_Crush Num_Crack_HF_Curves Size_Font 
global Plot_Aperture_Curves Plot_Pressure_Curves Num_Step_to_Plot 
global Plot_Velocity_Curves Plot_Quantity_Curves Plot_Concentr_Curves
global Title_Font Key_Figure_Control_Widget
global wbpt_Matrix avol_Matrix



%================================================
%  �����ѷ�ǰԵVertex����������Ӧ��S1����
%================================================
if Key_PLOT(6,2)==1
    disp('    > Plot S1 curve of crack vertices....')  
    c_crack = Key_PLOT(6,3);
	if c_crack > num_Crack
	    disp('      Error :: can not find this crack, error in Plot_3D_curves.m!')
		Error_Message
	end
	if c_crack <= 0
	    disp('      Error :: illegal crack number, error in Plot_3D_curves.m!')
		Error_Message
	end	
	%��ȡ�ļ����⻬���������ݣ�
	if exist([Full_Pathname,'.cvps_',num2str(POST_Substep)], 'file') ==2 
	    namefile= [Full_Pathname,'.cvps_',num2str(POST_Substep)];
	    data=fopen(namefile,'r'); 
		data2 =fopen(namefile,'r'); 
		yes_empty_file = fgetl(data2); 
		if yes_empty_file~= -1   %����ļ���Ϊ��
			lineNum = 0;
			while ~feof(data)
				lineNum = lineNum+1;
				TemData = fgetl(data);      
				if lineNum ==c_crack
					Plot_Vertex_S1 = str2num(TemData);
				end
			end
		else
	        Plot_Vertex_S1=[];
		    disp('      Warnning :: cvps file is empty, in Plot_3D_curves.m!')
		    % Error_Message		
		end
		fclose(data); 
		fclose(data2); 
	else
	    Plot_Vertex_S1=[];
		disp('      Error :: can not find cvps file, error in Plot_3D_curves.m!')
		Error_Message
	end 
	if Key_PLOT(6,4)==1 %���ƹ⻬����ǰ������
		%��ȡ�ļ����⻬����ǰ�����ݣ� NEWFTU-2022052001.
		if exist([Full_Pathname,'.cvpr_',num2str(POST_Substep)], 'file') ==2 
			namefile= [Full_Pathname,'.cvpr_',num2str(POST_Substep)];
			data=fopen(namefile,'r'); 
			data2 =fopen(namefile,'r'); 
			yes_empty_file = fgetl(data2); 
			if yes_empty_file~= -1   %����ļ���Ϊ��
				lineNum = 0;
				while ~feof(data)
					lineNum = lineNum+1;
					TemData = fgetl(data);      
					if lineNum ==c_crack
						Plot_Vertex_S1_Real = str2num(TemData);
					end
				end
			else
				Plot_Vertex_S1=[];
				disp('      Warnning :: cvpr file is empty, in Plot_3D_curves.m!')
				% Error_Message		
			end
			fclose(data); 
			fclose(data2); 
		else
			Plot_Vertex_S1_Real=[];
			disp('      Error :: can not find cvpr file, error in Plot_3D_curves.m!')
			Error_Message
		end	
	end
	%��������
	num_Vertex = size(Plot_Vertex_S1,2);
	c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
	hold on;
	string_title = ['S1 of vertex of crack ',num2str(c_crack),' (MPa)'];
	title(string_title,'FontName',Title_Font,'FontSize',Size_Font)
	Plot_x = 1:1:num_Vertex;
	plot(Plot_x ,Plot_Vertex_S1/1.0E6,'black-o','LineWidth',1)
	if Key_PLOT(6,4)==1 %���ƹ⻬����ǰ������
	    plot(Plot_x ,Plot_Vertex_S1_Real/1.0E6,'red--o','LineWidth',0.5)
		legend_FontSize = legend('After Smoothing','Before Smoothing');%ͼ������
	end
	grid on %������	
	set(gca,'box','on','fontname','Times Newroman','Fontsize',9);
	xlabel('Vertex number','FontName',Title_Font,'FontSize',Size_Font) 
	ylabel('S1','FontName',Title_Font,'FontSize',Size_Font) 	
	Save_Picture(c_figure,Full_Pathname,'curve_s1')
end


%================================================
%�����ѷ�ǰԵVertex����ĵ�ЧӦ��ǿ������Keq����
%================================================
if Key_PLOT(6,2)==2
    disp('    > Plot Keq curve of crack vertices....')  
    c_crack = Key_PLOT(6,3);
	if c_crack > num_Crack
	    disp('      Error :: can not find this crack, error in Plot_3D_curves.m!')
		Error_Message
	end
	if c_crack <= 0
	    disp('      Error :: illegal crack number, error in Plot_3D_curves.m!')
		Error_Message
	end	
	%��ȡ�ļ�
	if exist([Full_Pathname,'.cvke_',num2str(POST_Substep)], 'file') ==2 
	    namefile= [Full_Pathname,'.cvke_',num2str(POST_Substep)];
	    data=fopen(namefile,'r'); 
		data2 =fopen(namefile,'r'); 
		yes_empty_file = fgetl(data2); 
		if yes_empty_file~= -1   %����ļ���Ϊ��
			lineNum = 0;
			while ~feof(data)
				lineNum = lineNum+1;
				TemData = fgetl(data);      
				if lineNum ==c_crack
					Plot_Vertex_Keq = str2num(TemData);
				end
			end
		else
	        Plot_Vertex_Keq=[];
		    disp('      Warnning :: cvke file is empty, in Plot_3D_curves.m!')
		    % Error_Message		
		end
		fclose(data); 
		fclose(data2); 
	else
	    Plot_Vertex_Keq=[];
		disp('      Error :: can not find cvke file, error in Plot_3D_curves.m!')
		Error_Message
	end	
	%��������
	num_Vertex = size(Plot_Vertex_Keq,2);
	c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
	hold on;
	string_title = ['Keq of vertex of crack ',num2str(c_crack),' (MPa*m^1^/^2)'];
	title(string_title,'FontName',Title_Font,'FontSize',Size_Font)
	Plot_x = 1:1:num_Vertex;
	plot(Plot_x ,Plot_Vertex_Keq/1.0E6,'black-o','LineWidth',1)
	grid on %������	
	set(gca,'box','on','fontname','Times Newroman','Fontsize',9);
	xlabel('Vertex number','FontName',Title_Font,'FontSize',Size_Font) 
	ylabel('Keq','FontName',Title_Font,'FontSize',Size_Font) 	
	Save_Picture(c_figure,Full_Pathname,'curve_keq')
end

%================================================
%�����ѷ�ǰԵVertex�����I��ǿ����������
%================================================
if Key_PLOT(6,2)==3
    disp('    > Plot KI curve of crack vertices....')  
    c_crack = Key_PLOT(6,3);
	if c_crack > num_Crack
	    disp('      Error :: can not find this crack, error in Plot_3D_curves.m!')
		Error_Message
	end
	if c_crack <= 0
	    disp('      Error :: illegal crack number, error in Plot_3D_curves.m!')
		Error_Message
	end	
	%��ȡ�ļ�
	if exist([Full_Pathname,'.cvk1_',num2str(POST_Substep)], 'file') ==2 
	    namefile= [Full_Pathname,'.cvk1_',num2str(POST_Substep)];
	    data  =fopen(namefile,'r'); 
		data2 =fopen(namefile,'r'); 
		yes_empty_file = fgetl(data2); 		
		if yes_empty_file~= -1   %����ļ���Ϊ��		
			lineNum = 0;
			while ~feof(data)
				lineNum = lineNum+1;
				TemData = fgetl(data);      
				if lineNum ==c_crack
					Plot_Vertex_K1 = str2num(TemData);
				end
			end
		else
	        Plot_Vertex_K1=[];
		    disp('      Warnning :: cvk1 file is empty, in Plot_3D_curves.m!')
		    % Error_Message				
		end
		fclose(data); 
		fclose(data2); 
	else
	    Plot_Vertex_K1=[];
		disp('      Error :: can not find cvk1 file, error in Plot_3D_curves.m!')
		Error_Message
	end
	if Key_PLOT(6,4)==1 %���ƹ⻬����ǰ������	
		%��ȡ�ļ�
		if exist([Full_Pathname,'.ckr1_',num2str(POST_Substep)], 'file') ==2 
			namefile= [Full_Pathname,'.ckr1_',num2str(POST_Substep)];
			data  =fopen(namefile,'r'); 
			data2 =fopen(namefile,'r'); 
			yes_empty_file = fgetl(data2); 		
			if yes_empty_file~= -1   %����ļ���Ϊ��		
				lineNum = 0;
				while ~feof(data)
					lineNum = lineNum+1;
					TemData = fgetl(data);      
					if lineNum ==c_crack
						Plot_Vertex_K1_Real = str2num(TemData);
					end
				end
			else
				Plot_Vertex_K1_Real=[];
				disp('      Warnning :: ckr1 file is empty, in Plot_3D_curves.m!')
				% Error_Message				
			end
			fclose(data); 
			fclose(data2); 
		else
			Plot_Vertex_K1_Real=[];
			disp('      Error :: can not find ckr1 file, error in Plot_3D_curves.m!')
			Error_Message
		end	
	end
	
	%��������
	num_Vertex = size(Plot_Vertex_K1,2);
	c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
	hold on;
	string_title = ['KI of vertex of crack ',num2str(c_crack),' (MPa*m^1^/^2)'];
	title(string_title,'FontName',Title_Font,'FontSize',Size_Font)
	Plot_x = 1:1:num_Vertex;
	plot(Plot_x ,Plot_Vertex_K1/1.0E6,'black-o','LineWidth',1)
	if Key_PLOT(6,4)==1 %���ƹ⻬����ǰ������	
		plot(Plot_x ,Plot_Vertex_K1_Real/1.0E6,'black--o','LineWidth',0.5)
		legend_FontSize = legend('After Smoothing','Before Smoothing');%ͼ������
	end
	grid on %������
	set(gca,'box','on','fontname','Times Newroman','Fontsize',9);
	xlabel('Vertex number','FontName',Title_Font,'FontSize',Size_Font) 
	ylabel('KI','FontName',Title_Font,'FontSize',Size_Font) 
    Save_Picture(c_figure,Full_Pathname,'curve_k1')	
end

%================================================
%�����ѷ�ǰԵVertex�����II��ǿ����������
%================================================
if Key_PLOT(6,2)==4
    disp('    > Plot KII curve of crack vertices....')  
    c_crack = Key_PLOT(6,3);
	if c_crack > num_Crack
	    disp('      Error :: can not find this crack, error in Plot_3D_curves.m!')
		Error_Message
	end
	if c_crack <= 0
	    disp('      Error :: illegal crack number, error in Plot_3D_curves.m!')
		Error_Message
	end	
	%��ȡ�ļ�
	if exist([Full_Pathname,'.cvk2_',num2str(POST_Substep)], 'file') ==2 
	    namefile= [Full_Pathname,'.cvk2_',num2str(POST_Substep)];
	    data=fopen(namefile,'r'); 
		data2 =fopen(namefile,'r'); 
		yes_empty_file = fgetl(data2); 
		if yes_empty_file~= -1   %����ļ���Ϊ��				
			lineNum = 0;
			while ~feof(data)
				lineNum = lineNum+1;
				TemData = fgetl(data);      
				if lineNum ==c_crack
					Plot_Vertex_K2 = str2num(TemData);
				end
			end
		else
	        Plot_Vertex_K2=[];
		    disp('      Warnning :: cvk2 file is empty, in Plot_3D_curves.m!')
		    % Error_Message				
		end
		fclose(data); 
		fclose(data2); 
	else
	    Plot_Vertex_K2=[];
		disp('      Error :: can not find cvk2 file, error in Plot_3D_curves.m!')
		Error_Message
	end
	if Key_PLOT(6,4)==1 %���ƹ⻬����ǰ������	
		%��ȡ�ļ�
		if exist([Full_Pathname,'.ckr2_',num2str(POST_Substep)], 'file') ==2 
			namefile= [Full_Pathname,'.ckr2_',num2str(POST_Substep)];
			data  =fopen(namefile,'r'); 
			data2 =fopen(namefile,'r'); 
			yes_empty_file = fgetl(data2); 		
			if yes_empty_file~= -1   %����ļ���Ϊ��		
				lineNum = 0;
				while ~feof(data)
					lineNum = lineNum+1;
					TemData = fgetl(data);      
					if lineNum ==c_crack
						Plot_Vertex_K2_Real = str2num(TemData);
					end
				end
			else
				Plot_Vertex_K2_Real=[];
				disp('      Warnning :: ckr1 file is empty, in Plot_3D_curves.m!')
				% Error_Message				
			end
			fclose(data); 
			fclose(data2); 
		else
			Plot_Vertex_K2_Real=[];
			disp('      Error :: can not find ckr1 file, error in Plot_3D_curves.m!')
			Error_Message
		end	
	end	
	%��������
	num_Vertex = size(Plot_Vertex_K2,2);
	c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
	hold on;
	string_title = ['KII of vertex of crack ',num2str(c_crack),' (MPa*m^1^/^2)'];
	title(string_title,'FontName',Title_Font,'FontSize',Size_Font)
	Plot_x = 1:1:num_Vertex;
	plot(Plot_x ,Plot_Vertex_K2/1.0E6,'black-o','LineWidth',1)
	if Key_PLOT(6,4)==1 %���ƹ⻬����ǰ������	
		plot(Plot_x ,Plot_Vertex_K2_Real/1.0E6,'black--o','LineWidth',0.5)
		legend_FontSize = legend('After Smoothing','Before Smoothing');%ͼ������
	end	
    grid on %������
	set(gca,'box','on','fontname','Times Newroman','Fontsize',9);
	xlabel('Vertex number','FontName',Title_Font,'FontSize',Size_Font) 
	ylabel('KII','FontName',Title_Font,'FontSize',Size_Font) 	
	Save_Picture(c_figure,Full_Pathname,'curve_k2')
end


%================================================
%�����ѷ�ǰԵVertex�����III��ǿ����������
%================================================
if Key_PLOT(6,2)==5
    disp('    > Plot KIII curve of crack vertices....')  
    c_crack = Key_PLOT(6,3);
	if c_crack > num_Crack
	    disp('      Error :: can not find this crack, error in Plot_3D_curves.m!')
		Error_Message
	end
	if c_crack <= 0
	    disp('      Error :: illegal crack number, error in Plot_3D_curves.m!')
		Error_Message
	end	
	%��ȡ�ļ�
	if exist([Full_Pathname,'.cvk3_',num2str(POST_Substep)], 'file') ==2 
	    namefile= [Full_Pathname,'.cvk3_',num2str(POST_Substep)];
	    data=fopen(namefile,'r'); 
	    data2=fopen(namefile,'r'); 
		yes_empty_file = fgetl(data2); 
		if yes_empty_file~= -1   %����ļ���Ϊ��		
			lineNum = 0;
			while ~feof(data)
				lineNum = lineNum+1;
				TemData = fgetl(data);      
				if lineNum ==c_crack
					Plot_Vertex_K3 = str2num(TemData);
				end
			end
		else
	        Plot_Vertex_K3=[];
		    disp('      Warnning :: cvk2 file is empty, in Plot_3D_curves.m!')
		    % Error_Message				
		end			
		fclose(data); 
		fclose(data2); 
	else
	    Plot_Vertex_K3=[];
		disp('      Error :: can not find cvps file, error in Plot_3D_curves.m!')
		Error_Message
	end
	if Key_PLOT(6,4)==1 %���ƹ⻬����ǰ������	
		%��ȡ�ļ�
		if exist([Full_Pathname,'.ckr3_',num2str(POST_Substep)], 'file') ==2 
			namefile= [Full_Pathname,'.ckr3_',num2str(POST_Substep)];
			data  =fopen(namefile,'r'); 
			data2 =fopen(namefile,'r'); 
			yes_empty_file = fgetl(data2); 		
			if yes_empty_file~= -1   %����ļ���Ϊ��		
				lineNum = 0;
				while ~feof(data)
					lineNum = lineNum+1;
					TemData = fgetl(data);      
					if lineNum ==c_crack
						Plot_Vertex_K3_Real = str2num(TemData);
					end
				end
			else
				Plot_Vertex_K3_Real=[];
				disp('      Warnning :: ckr1 file is empty, in Plot_3D_curves.m!')
				% Error_Message				
			end
			fclose(data); 
			fclose(data2); 
		else
			Plot_Vertex_K3_Real=[];
			disp('      Error :: can not find ckr1 file, error in Plot_3D_curves.m!')
			Error_Message
		end	
	end	
	%��������
	num_Vertex = size(Plot_Vertex_K3,2);
	c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
	hold on;
	string_title = ['KII of vertex of crack ',num2str(c_crack),' (MPa*m^1^/^2)'];
	title(string_title,'FontName',Title_Font,'FontSize',Size_Font)
	Plot_x = 1:1:num_Vertex;
	plot(Plot_x ,Plot_Vertex_K3/1.0E6,'black-o','LineWidth',1)
	if Key_PLOT(6,4)==1 %���ƹ⻬����ǰ������	
		plot(Plot_x ,Plot_Vertex_K3_Real/1.0E6,'black--o','LineWidth',0.5)
		legend_FontSize = legend('After Smoothing','Before Smoothing');%ͼ������
	end	
	grid on %������
	set(gca,'box','on','fontname','Times Newroman','Fontsize',9);
	xlabel('Vertex number','FontName',Title_Font,'FontSize',Size_Font) 
	ylabel('KIII','FontName',Title_Font,'FontSize',Size_Font) 
    Save_Picture(c_figure,Full_Pathname,'curve_k3')	
end

%================================================
%�����ѷ�ǰԵVertex�����I+II+III��ǿ����������
%================================================
if Key_PLOT(6,2)==6
    disp('    > Plot KI, KII, and KIII curves of crack vertices....')  
    c_crack = Key_PLOT(6,3);
	if c_crack > num_Crack
	    disp('      Error :: can not find this crack, error in Plot_3D_curves.m!')
		Error_Message
	end
	if c_crack <= 0
	    disp('      Error :: illegal crack number, error in Plot_3D_curves.m!')
		Error_Message
	end	
	%��ȡ�ļ�(�⻬������)
	if exist([Full_Pathname,'.cvk1_',num2str(POST_Substep)], 'file') ==2 
	    namefile= [Full_Pathname,'.cvk1_',num2str(POST_Substep)];
	    data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);      
			if lineNum ==c_crack
				Plot_Vertex_K1 = str2num(TemData);
			end
		end
		fclose(data); 
	else
	    Plot_Vertex_K1=[];
		disp('      Error :: can not find cvps file, error in Plot_3D_curves.m!')
		Error_Message
	end	
	if exist([Full_Pathname,'.cvk2_',num2str(POST_Substep)], 'file') ==2 
	    namefile= [Full_Pathname,'.cvk2_',num2str(POST_Substep)];
	    data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);      
			if lineNum ==c_crack
				Plot_Vertex_K2 = str2num(TemData);
			end
		end
		fclose(data); 
	else
	    Plot_Vertex_K2=[];
		disp('      Error :: can not find cvps file, error in Plot_3D_curves.m!')
		Error_Message
	end	
	if exist([Full_Pathname,'.cvk3_',num2str(POST_Substep)], 'file') ==2 
	    namefile= [Full_Pathname,'.cvk3_',num2str(POST_Substep)];
	    data=fopen(namefile,'r'); 
		lineNum = 0;
		while ~feof(data)
			lineNum = lineNum+1;
			TemData = fgetl(data);      
			if lineNum ==c_crack
				Plot_Vertex_K3 = str2num(TemData);
			end
		end
		fclose(data); 
	else
	    Plot_Vertex_K3=[];
		disp('      Error :: can not find ckr1 file, error in Plot_3D_curves.m!')
		Error_Message
	end
	%��ȡ�ļ�(�⻬����ǰ��)
	if Key_PLOT(6,4)==1 %���ƹ⻬����ǰ������
		if exist([Full_Pathname,'.ckr1_',num2str(POST_Substep)], 'file') ==2 
			namefile= [Full_Pathname,'.ckr1_',num2str(POST_Substep)];
			data=fopen(namefile,'r'); 
			lineNum = 0;
			while ~feof(data)
				lineNum = lineNum+1;
				TemData = fgetl(data);      
				if lineNum ==c_crack
					Plot_Vertex_K1_Real = str2num(TemData);
				end
			end
			fclose(data); 
		else
			Plot_Vertex_K1=[];
			disp('      Error :: can not find ckr2 file, error in Plot_3D_curves.m!')
			Error_Message
		end	
		if exist([Full_Pathname,'.ckr2_',num2str(POST_Substep)], 'file') ==2 
			namefile= [Full_Pathname,'.ckr2_',num2str(POST_Substep)];
			data=fopen(namefile,'r'); 
			lineNum = 0;
			while ~feof(data)
				lineNum = lineNum+1;
				TemData = fgetl(data);      
				if lineNum ==c_crack
					Plot_Vertex_K2_Real = str2num(TemData);
				end
			end
			fclose(data); 
		else
			Plot_Vertex_K2=[];
			disp('      Error :: can not find ckr3 file, error in Plot_3D_curves.m!')
			Error_Message
		end	
		if exist([Full_Pathname,'.ckr3_',num2str(POST_Substep)], 'file') ==2 
			namefile= [Full_Pathname,'.ckr3_',num2str(POST_Substep)];
			data=fopen(namefile,'r'); 
			lineNum = 0;
			while ~feof(data)
				lineNum = lineNum+1;
				TemData = fgetl(data);      
				if lineNum ==c_crack
					Plot_Vertex_K3_Real = str2num(TemData);
				end
			end
			fclose(data); 
		else
			Plot_Vertex_K3=[];
			disp('      Error :: can not find cvps file, error in Plot_3D_curves.m!')
			Error_Message
		end	
	end
	%��������
	num_Vertex = size(Plot_Vertex_K3,2);
	c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
	hold on;
	string_title = ['K of vertex of crack ',num2str(c_crack),' (MPa*m^1^/^2)'];
	title(string_title,'FontName',Title_Font,'FontSize',Size_Font)
	Plot_x = 1:1:num_Vertex;
	plot(Plot_x ,Plot_Vertex_K1/1.0E6,'black-o','LineWidth',1)
	plot(Plot_x ,Plot_Vertex_K2/1.0E6,'blue-o','LineWidth',1)
	plot(Plot_x ,Plot_Vertex_K3/1.0E6,'red-o','LineWidth',1)
	if Key_PLOT(6,4)==1 %���ƹ⻬����ǰ������
		plot(Plot_x ,Plot_Vertex_K1_Real/1.0E6,'black--o','LineWidth',0.5)
		plot(Plot_x ,Plot_Vertex_K2_Real/1.0E6,'blue--o','LineWidth',0.5)
		plot(Plot_x ,Plot_Vertex_K3_Real/1.0E6,'red--o','LineWidth',0.5)	
	end
	xlabel('Vertex number','FontName',Title_Font,'FontSize',Size_Font) 
	ylabel('K','FontName',Title_Font,'FontSize',Size_Font) 	
	grid on %������
	set(gca,'box','on','fontname','Times Newroman','Fontsize',9);
	if Key_PLOT(6,4)==1 %���ƹ⻬����ǰ������
	    legend_FontSize = legend('KI-after smoothing','KII-after smoothing','KIII-after smoothing','KI-before smoothing','KII-before smoothing','KIII-before smoothing');%ͼ������
	else 
	    legend_FontSize = legend('KI','KII','KIII');%ͼ������
	end
    set(legend_FontSize,'FontSize',10.5) 
	Save_Picture(c_figure,Full_Pathname,'curve_k1_k2_k3')
end

%==================================================================
%����HF������ѹ��-ʱ������/ѹ��-��չ������. NEWFTU-2022070301.
%==================================================================
if Key_PLOT(6,2)==7 | Key_PLOT(6,2)==8
    if Key_PLOT(6,2)==7
        disp('    > Plot Pressure-Time curves...')  
	elseif Key_PLOT(6,2)==8 %NEWFTU-2022101101.
	    disp('    > Plot Pressure curves...')  
	end
	if isempty(wbpt_Matrix) ==0 
	    %����׼��.
	    c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
		hold on;
		if Key_PLOT(6,2)==7
			string_title = ['Pressure-time curves of wellbore(s)']; 
		elseif Key_PLOT(6,2)==8 %NEWFTU-2022101101.
			string_title = ['Pressure curves of wellbore(s)']; 
		end
	    
	    title(string_title,'FontName',Title_Font,'FontSize',Size_Font)
		if Key_PLOT(6,2)==7
		    if Key_PLOT(6,6)==1
			    xlabel('Time (mins)','FontName',Title_Font,'FontSize',Size_Font)
			elseif(Key_PLOT(6,6)==2)
			    xlabel('Time (s)','FontName',Title_Font,'FontSize',Size_Font)
			else
			    disp('    > ERROR-2023091901 :: illegal Key_PLOT(6,2)! Key_PLOT(6,2) should be 1(s) or 2(mins)!') 
			end
		elseif Key_PLOT(6,2)==8 %NEWFTU-2022101101.
			xlabel('Propagation step ','FontName',Title_Font,'FontSize',Size_Font) 
		end		
	    
	    ylabel('Pressure (MPa)','FontName',Title_Font,'FontSize',Size_Font) 	
	    grid on %������
	    set(gca,'box','on','fontname','Times Newroman','Fontsize',9);

		
	    %����ܵľ�Ͳ��num_WB.
	    num_WB = max(wbpt_Matrix(:,1));
		
		num_Lines_of_wbpt = size(wbpt_Matrix,1);
	    
		
		Wb_num_to_Plot = Key_PLOT(6,5); %�����Ƶľ�Ͳ��.
		if Wb_num_to_Plot > num_WB
		    disp('      Error :: illegal Wb_num_to_Plot, check Key_PLOT(6,5), error in Plot_3D_curves.m!')
		    Error_Message
		end
		
		if Key_PLOT(6,5)==-999   %����ȫ����Ͳ.
			i_WB_Start = 1;
			i_WB_End   = num_WB;
		else                    %����ָ����Ͳ.
			i_WB_Start = Wb_num_to_Plot;
			i_WB_End   = Wb_num_to_Plot;
		end
		
		%������Ͳѭ��.
		stage_Count = 0;
        myMap = rand(1000,3);  %�����ɫ		
		for i_WB = i_WB_Start:i_WB_End 
			%��õ�ǰ��Ͳ�ķֶ���Stages_num_of_cWB.
			c_count = 0;
			c_WB_Matrx = [];
			for c_i= 1:num_Lines_of_wbpt
				if wbpt_Matrix(c_i,1)==i_WB
					c_count = c_count +1;
					c_WB_Matrx(c_count,1:4) = wbpt_Matrix(c_i,2:5);
				end
			end
			Stages_num_of_cWB = max(c_WB_Matrx(1:c_count,1));
			%==========================================
			%   �����ֶ�ѭ����ÿһ���ֶλ�һ������.
			%==========================================
			for i_Stage = 1:Stages_num_of_cWB
				%��õ�ǰ�ֶε����Ѳ���
				c_count = 0;
				c_Stage_Matrx = [];
				for c_i= 1:size(c_WB_Matrx,1)
					if c_WB_Matrx(c_i,1)==i_Stage
						c_count = c_count +1;
						c_Stage_Matrx(c_count,1:3) = c_WB_Matrx(c_i,2:4);
					end
				end	
				Prop_num_of_cStage = max(c_Stage_Matrx(1:c_count,1));
				%��������.
				stage_Count = stage_Count +1;
				
				if Key_PLOT(6,2)==7
					if Key_PLOT(6,6)==1
					    plot(c_Stage_Matrx(1:c_count,2)/60.0,c_Stage_Matrx(1:c_count,3)/1.0E6,'-o','color',myMap(stage_Count,1:3),'LineWidth',0.5)  %����Ϊ��λ
					elseif(Key_PLOT(6,6)==2)
					    plot(c_Stage_Matrx(1:c_count,2),c_Stage_Matrx(1:c_count,3)/1.0E6,'-o','color',myMap(stage_Count,1:3),'LineWidth',0.5)     %sΪ��λ
					end
				elseif Key_PLOT(6,2)==8 %NEWFTU-2022101101.
					plot([1:c_count],c_Stage_Matrx(1:c_count,3)/1.0E6,'-o','color',myMap(stage_Count,1:3),'LineWidth',0.5)  %����Ϊ��λ
				end
				
				Legend_Str(stage_Count) = {['Stage ',num2str(i_Stage),' of Wellbore ',num2str(i_WB)]};		
			end
		end
		c_legend = legend(Legend_Str(1:stage_Count));%ͼ������
		set(c_legend,'FontSize',10.5) 			
	else    
		disp('      Error :: can not find wbpt file, error in Plot_3D_curves.m!')
		Error_Message
	end
    Save_Picture(c_figure,Full_Pathname,'curve_pressure_time')
end

%==================================================
%������չ��-�ѷ��������; NEWFTU-2022110901.
%==================================================
if Key_PLOT(6,2)==9 
	disp('    > Plot crack volume curve...')  
	if isempty(avol_Matrix) ==0 
	    %����׼��.
	    c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
		hold on;
		string_title = ['Crack volume curves']; 
	    
	    title(string_title,'FontName',Title_Font,'FontSize',Size_Font)

		xlabel('Propagation step ','FontName',Title_Font,'FontSize',Size_Font) 	
	    
	    ylabel('Volume (m3)','FontName',Title_Font,'FontSize',Size_Font) 	
		
	    grid on %������
	    set(gca,'box','on','fontname','Times Newroman','Fontsize',9);
        
		%������Ѳ�����Ŀ.
	    num_Steps = size(avol_Matrix,1);
		
		% ȫ�����.
		plot([1:size(avol_Matrix)],avol_Matrix(:,2),'-o','color','black','LineWidth',0.5)  
		Legend_Str(1) = {['Sum of all crack volume']};
		
		% �������.
		plot([1:size(avol_Matrix)],avol_Matrix(:,3),'-o','color','red','LineWidth',0.5)  
		Legend_Str(2) = {['Sum of positive crack volume']};
		
		c_legend = legend(Legend_Str(1:2));%ͼ������
		set(c_legend,'FontSize',10.5,'Location','southeast') 			
	else    
		disp('      Error :: can not find avol file, error in Plot_3D_curves.m!')
		Error_Message
	end
    Save_Picture(c_figure,Full_Pathname,'curve_pressure_time')
end

%==================================================================
%����HF������ѹ��-ʱ������/ѹ��-��չ������. NEWFTU-2022070301.
%==================================================================
if Key_PLOT(6,2)==10
    disp('    > Plot pressure-time curve for the hydraulic fracturing simulation...')  
    namefile= [Full_Pathname,'.injp'];
	if exist(namefile, 'file') ==2 
	    disp('    > Plot pressure-time curve for the hydraulic fracturing simulation...')  
		Time_Pressure_Matrix   = load([namefile]);
		% ����׼��.
		c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on');
		hold on;
		disp('    > Plot Pressure-Time curves...')  
		string_title = ['Pressure-time curve']; 
		title(string_title,'FontName',Title_Font,'FontSize',Size_Font)
		xlabel('Time (s)','FontName',Title_Font,'FontSize',Size_Font) 
		ylabel('Pressure (MPa)','FontName',Title_Font,'FontSize',Size_Font) 	
		grid on %������
		set(gca,'box','on','fontname','Times Newroman','Fontsize',9);
		plot(Time_Pressure_Matrix(:,1),Time_Pressure_Matrix(:,2)/1.0e6,'-o','color','blue','LineWidth',0.5)  %����Ϊ��λ
		Save_Picture(c_figure,Full_Pathname,'curve_pressure_time') 
	else
		disp('      Error :: can not find injp file, error in Plot_3D_curves.m!')
		Error_Message
	end
end