function varargout = PhiPsi_2D_Post_1_GUI(varargin)
% PHIPSI_2D_POST_1_GUI MATLAB code for PhiPsi_2D_Post_1_GUI.fig
%      PHIPSI_2D_POST_1_GUI, by itself, creates a new PHIPSI_2D_POST_1_GUI or raises the existing
%      singleton*.
%
%      H = PHIPSI_2D_POST_1_GUI returns the handle to a new PHIPSI_2D_POST_1_GUI or the handle to
%      the existing singleton*.
%
%      PHIPSI_2D_POST_1_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHIPSI_2D_POST_1_GUI.M with the given input arguments.
%
%      PHIPSI_2D_POST_1_GUI('Property','Value',...) creates a new PHIPSI_2D_POST_1_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PhiPsi_2D_Post_1_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PhiPsi_2D_Post_1_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PhiPsi_2D_Post_1_GUI

% Last Modified by GUIDE v2.5 01-Oct-2018 12:54:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PhiPsi_2D_Post_1_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @PhiPsi_2D_Post_1_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PhiPsi_2D_Post_1_GUI is made visible.
function PhiPsi_2D_Post_1_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PhiPsi_2D_Post_1_GUI (see VARARGIN)

% Choose default command line output for PhiPsi_2D_Post_1_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PhiPsi_2D_Post_1_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PhiPsi_2D_Post_1_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%..........................................................
%                 ִ�з���
%..........................................................
function pushbutton_apply_Callback(hObject, eventdata, handles)
%global Input_Num_Step_to_Plot Input_Defor_Factor Input_KIKII_Crack Input_KIKII_Tip
global Input_Key_PLOT
% ��1��,����Ԫ����: Mesh(1),Node(2),El(3),Gauss points(4),
%                   5: �ѷ켰�ѷ������(=1,�����ѷ�;=2,�����ѷ켰�����),
%                   6: ����㼰����(=1,�����;=2,�����ͱ��),
%                   7: �ѷ�ڵ�(�����)���(=1,�ڵ㼯��ǿ�ڵ��غ�;=2,����㾻ˮѹ;=3,���������;=4,����㿪��;=5,�����ճ����x�������;=6,�����ճ����y�������),
%                   ��ǿ�ڵ�(8),������(9),
%                   ֧�ż�(10),��ԪӦ��״̬�Ƿ��1-��3>Tol(11),��Ȼ�ѷ�(12),��Ԫ�Ӵ�״̬/ճ���ѷ�״̬(13),�ѷ���(14),Fracture zone(15)
% ��2��,�������ͼ: Deformation(1),Node(2),El(3),Gauss points(4),Crack(5:1,line;2,shape),Scaling Factor(6),
%                   FsBs(7=1or2or3),Deformed and undefor(8),Blank(9),֧�ż�(10),Blank(11),Blank(12),
%                   ��Ԫ�Ӵ�״̬/ճ���ѷ�״̬(13),��ǿ�ڵ�(14),Fracture zone(15)
% ��3��,Ӧ����ͼ:   Stress Contour(1,2:Gauss points),(2:1,Only Mises stress;2,��xӦ��;3,��yӦ��;4,����Ӧ��),
%                   ��Ӧ��(3:1,��Ӧ��;2,��Ӧ��+����;3,�������Ӧ��),�������(4:1,��Ч����Ӧ��),
%                   Crack(5:1,line;2,shape),Scaling Factor(6),FsBs(7=1or2or3),
%                   undeformed or Deformed(8),mesh(9),֧�ż�(10),Blank(11),Blank(12),Blank(13),Blank(14),Fracture zone(15)
% ��4��,λ����ͼ:   Plot Dis-Contour(1,2:Gauss points),Blank(2),Blank(3),Blank(4),Crack(5:1,line;2,shape),Scaling Factor(6),
%                   FsBs(7=1or2or3),undeformed or Deformed(8),mesh(9),֧�ż�(10),Blank(11),Blank(12),Blank(13),Blank(14),Fracture zone(15)
% ��5��,��������ͼ: Plot Contour(1:��ֵ,2:����ʸ��,3:��ֵ+����ʸ��,4:Gauss�㳡ֵ),Blank(2),Blank(3),Blank(4),Crack(5:1,line;2,shape),Scaling Factor(6),
%                   Blank(7),undeformed or Deformed(8),mesh(9),֧�ż�(10),Blank(11),Blank(12),Blank(13),Blank(14),Fracture zone(15)
% ��6��,���Ӷ���ѧ: Plot Contour(1:λ��,2:�ٶ�,3:),�켣(2=1,ȫ���켣;>1,�켣�����ļ��㲽��),Blank(3),Blank(4),Blank(4),Scaling Factor(6),
%                   Blank(7),Original location(8),plot a mesh(9),Blank(10),Blank(11),Blank(12),Blank(13),Blank(14),Box(15)
% ��7��,���߻���:   ���߻�����أ�����(1=1,���������),Ӧ��ǿ����������(2=1),Ӧ��ǿ�����������ѷ��(3),Ӧ��ǿ�����������Ѽ��(4),
%                   ˮ��ѹ�ѷ������������ߵĻ���(5=1,�����ѷ��ھ�ˮѹ�ֲ�����;5=2,�ѷ쿪�ȷֲ����ߣ�=3���ѷ����򿪶ȷֲ����ߣ�
%                                                  5=4���ѷ�����ڵ����ٷֲ�����;5=5���ѷ�����ڵ������ֲ�����;5=6���ѷ�����ڵ㴦֧�ż�Ũ�������ֲ�����),
%                   ˮ��ѹ�ѷ������������߶�Ӧ���ѷ��(6),ˮ��ѹ�ѷ���עˮ��ѹ���仯����(7=1,��ʱ��Ϊx������;7=2,��ѹ�Ѳ�Ϊx��),
%                   ����ҳ���������ʱ仯����(8=1),�Ƿ�����ۻ������仯����(9=1),�Ƿ����ĳһ�����ѹ���仯(10=1),
%                   �����غ�λ������(11=1,�غ�λ�����ߵĻ���*.fdcu�ļ�,),Blank(12-15)
%                         1   2   3   4   5   6              7   8   9  10  11  12  13  14   15
% Key_PLOT(1,:)         = [ 1,  0,  0,  0,  2,  0,             0,  1,  1  ,1  ,0  ,1  ,1  ,0  ,1];  
% Key_PLOT(2,:)         = [ 1,  0,  0,  0,  2,  Defor_Factor,  3,  0,  1  ,1  ,1  ,1  ,1  ,1  ,1];  
% Key_PLOT(3,:)         = [ 2,  1,  0,  0,  2,  Defor_Factor,  0,  1,  1  ,0  ,0  ,1  ,0  ,0  ,1];  
% Key_PLOT(4,:)         = [ 2,  0,  0,  0,  2,  Defor_Factor,  0,  1,  1  ,1  ,0  ,1  ,0  ,0  ,1];
% Key_PLOT(5,:)         = [ 0,  0,  0,  0,  1,  Defor_Factor,  0,  0,  1  ,1  ,0  ,1  ,0  ,0  ,1];   

%��ȡ�غɲ���
Input_Num_Step_to_Plot = str2num(get(handles.edit_step_number,'string'));
%��ȡ���ηŴ�ϵ��
Input_Defor_Factor = str2num(get(handles.edit_deform_factor,'string'));
Input_KIKII_Crack  = str2num(get(handles.edit_crack_num_to_plot_KIKII,'string'));
Input_KIKII_Tip    = str2num(get(handles.edit_tip_num_to_plot_KIKII,'string'));
%��ȡInput_Filename��Input_Work_Dirctory
Input_Filename= get(handles.edit_filename,'string');
Input_Work_Dirctory= get(handles.edit_directory,'string');

%��ȡ�������Ʋ���
set(handles.edit_log,'String',['Reading settings......'])

Input_Plot_Mesh = get(handles.checkbox_plot_mesh,'value');
Input_Plot_deformed_mesh = get(handles.checkbox_plot_deformed_mesh,'value');
Input_plot_stress = get(handles.checkbox_plot_stress,'value');
Input_plot_displacement = get(handles.checkbox_plot_displacement,'value');
Input_plot_elem_number = get(handles.checkbox_plot_elem_number,'value');
Input_plot_node_number = get(handles.checkbox_plot_node_number,'value');
Input_plot_crack_number = get(handles.checkbox_plot_crack_number,'value');
Input_mesh_line = get(handles.checkbox_mesh_line,'value');
Input_Gauss_Points = get(handles.checkbox_Gauss_Points,'value');
Input_undeformed_mesh = get(handles.checkbox_undeformed_mesh,'value');
Input_enriched_nodes = get(handles.checkbox_enriched_nodes,'value');
Input_boundary_conditions = get(handles.checkbox_boundary_conditions,'value');
Input_external_forces = get(handles.checkbox_external_forces,'value');
Input_checkbox_real_shaped_crack = get(handles.checkbox_real_shaped_crack,'value');
Input_plot_KIKII = get(handles.checkbox_plot_KIKII,'value');
Input_Contour_jet = get(handles.radiobutton_contour_jet,'value');
Input_Contour_summer = get(handles.radiobutton_contour_summer,'value');
Input_Contour_gray = get(handles.radiobutton_contour_gray,'value');
Input_Contour_autumn = get(handles.radiobutton_contour_autumn,'value');
%����Input_Key_PLOT
Input_Key_PLOT(1:7,1:15) =0;
%yyy(1:15) = [ 1,  0,  0,  0,  2,  0,             0,  1,  1  ,1  ,0  ,1  ,1  ,0  ,1];  
Input_Key_PLOT(1,1:15) = [Input_Plot_Mesh, Input_plot_node_number, Input_plot_elem_number, Input_Gauss_Points,  2, ...
                          0,0,  Input_enriched_nodes,  Input_Plot_Mesh  ,0  , ...
                          0,1,1,Input_plot_crack_number ,1];
if Input_checkbox_real_shaped_crack==1
    REAL_CRACK =2;
else
    REAL_CRACK =1;
end
if Input_boundary_conditions==1 && Input_external_forces==0
    BOUNDARY_FORCE =1;
elseif Input_boundary_conditions==0 && Input_external_forces==1
    BOUNDARY_FORCE =2;
elseif Input_boundary_conditions==1 && Input_external_forces==1
    BOUNDARY_FORCE =3;
else
    BOUNDARY_FORCE =0;
end
if Input_plot_stress==1;PLOT_STRESS =2;else;PLOT_STRESS =0;end
if Input_plot_displacement==1;PLOT_DISP =2;else;PLOT_DISP =0;end
if Input_Contour_jet==1;Input_Contour_String ='jet';end
if Input_Contour_summer==1;Input_Contour_String ='summer';end
if Input_Contour_gray==1;Input_Contour_String ='gray';end
if Input_Contour_autumn==1;Input_Contour_String ='autumn';end

Input_Key_PLOT(2,1:15) = [Input_Plot_deformed_mesh, Input_plot_node_number, Input_plot_elem_number, Input_Gauss_Points,  REAL_CRACK, ...
                          Input_Defor_Factor,BOUNDARY_FORCE, Input_undeformed_mesh, 1  ,1  , ...
                          0,1,1,Input_enriched_nodes ,1];  
Input_Key_PLOT(3,1:15) = [PLOT_STRESS, 6, 0, 0,  REAL_CRACK, ...
                          Input_Defor_Factor,BOUNDARY_FORCE, 1, Input_mesh_line  ,1  , ...
                          0,0,0,0 ,1];      
Input_Key_PLOT(4,1:15) = [PLOT_DISP, 0, 0, 0,  REAL_CRACK, ...
                          Input_Defor_Factor,BOUNDARY_FORCE, 1, Input_mesh_line  ,1  , ...
                          0,0,0,0 ,1];  
Input_Key_PLOT(7,1:15) = [1,Input_plot_KIKII,Input_KIKII_Crack, Input_KIKII_Tip,  0, ...
                          0,0, 0, 0  ,0  , ...
                          0,0,0,0 ,1];                        
% Input_Key_PLOT(2,1:15) = 0
% Input_Key_PLOT(3,1:15) = 0
% Input_Key_PLOT(4,1:15) = 0
% Input_Key_PLOT(5,1:15) = 0
% ִ�к������
set(handles.edit_log,'String',['Starting post-processing......'])
PhiPsi_2D_Post_1_Function(Input_Defor_Factor,Input_Num_Step_to_Plot,Input_Filename, ...
                          Input_Work_Dirctory,Input_Key_PLOT,Input_KIKII_Crack,Input_KIKII_Tip, ...
                          Input_Contour_String)
set(handles.edit_log,'String',['Post-processing done, see Command Window.log for more info.'])
% ��������
sound(sin(2*pi*20*(1:1500)/100));
warndlg('POST-PROCESSING DONE!','NOTICE')


                      
                      
% ------------------------------------------------------                      
% Cancel��ť
% ------------------------------------------------------   
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
exit


function edit_directory_Callback(hObject, eventdata, handles)
% hObject    handle to edit_directory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_directory as text
%        str2double(get(hObject,'String')) returns contents of edit_directory as a double


% --- Executes during object creation, after setting all properties.
function edit_directory_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_directory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_filename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filename as text
%        str2double(get(hObject,'String')) returns contents of edit_filename as a double


% --- Executes during object creation, after setting all properties.
function edit_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%..........................................................
%                 �򿪽���ļ���ز���
%..........................................................
function pushbutton_open_file_Callback(hObject, eventdata, handles)
%global Open_Filename, Open_Pathname,Open_Full_Name
% hObject    handle to pushbutton_open_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%��*.elem�ļ�
[Open_Filename, Open_Pathname] = uigetfile({'*.elem'},'Please select PhiPsi *.elem file');
% ���δ�ҵ��ļ�,��
if Open_Pathname==0
    set(handles.edit_directory,'String','Result files not selected yet, please open file or enter directly.');  % ��·����ʾ���ı�����
end
if Open_Filename==0
    set(handles.edit_filename,'String','Result files not selected yet, please open file or enter directly.');   % ��·����ʾ���ı�����
end

% ������ļ��Ƿ����
if Open_Filename~=0
    if Open_Pathname~=0
    Open_Full_Name = [Open_Pathname,Open_Filename];
    set(handles.edit_log,'String',['Result files located: ',Open_Full_Name])
    Open_Pathname(end:end)=[];
    Open_Filename(end-4:end)=[];
    set(handles.edit_directory,'String',Open_Pathname);   % ��·����ʾ���ı�����
    set(handles.edit_filename,'String',Open_Filename);    % ��·����ʾ���ı�����    
    end
end



%���ü��㲽��ز���
function edit_step_number_Callback(hObject, eventdata, handles)
%set(handles.edit_log,'String',['Step number setting done.'])


% --- Executes during object creation, after setting all properties.
function edit_step_number_CreateFcn(hObject, eventdata, handles)
global Input_Num_Step_to_Plot
% hObject    handle to edit_step_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_log_Callback(hObject, eventdata, handles)
% hObject    handle to edit_log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_log as text
%        str2double(get(hObject,'String')) returns contents of edit_log as a double


% --- Executes during object creation, after setting all properties.
function edit_log_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_deform_factor_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_deform_factor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_save_settings.
function pushbutton_save_settings_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_plot_mesh.
function checkbox_plot_mesh_Callback(hObject, eventdata, handles)


% --- Executes on button press in checkbox_plot_node_number.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_node_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_node_number


% --- Executes on button press in checkbox_crack_number.
function checkbox_crack_number_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_crack_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_crack_number


% --- Executes on button press in checkbox_plot_elem_number.
function checkbox_plot_elem_number_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_elem_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_elem_number


% --- Executes on button press in checkbox_plot_deformed_mesh.
function checkbox_plot_deformed_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_deformed_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_deformed_mesh


% --- Executes on button press in checkbox_plot_stress.
function checkbox_plot_stress_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_stress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_stress


% --- Executes on button press in checkbox_plot_displacement.
function checkbox_plot_displacement_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_displacement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_displacement


% --- Executes on button press in checkbox_enriched_nodes.
function checkbox_enriched_nodes_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_enriched_nodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_enriched_nodes


% --- Executes on button press in checkbox_undeformed_mesh.
function checkbox_undeformed_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_undeformed_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_undeformed_mesh


% --- Executes on button press in checkbox_Gauss_Points.
function checkbox_Gauss_Points_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Gauss_Points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Gauss_Points


% --- Executes on button press in checkbox_mesh_line.
function checkbox_mesh_line_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_mesh_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_mesh_line


% --- Executes on button press in checkbox_plot_node_number.
function checkbox_plot_node_number_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_node_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_node_number


% --- Executes on button press in checkbox_plot_crack_number.
function checkbox_plot_crack_number_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_crack_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_crack_number


% --- Executes on button press in checkbox_elem_number.
function checkbox_elem_number_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_elem_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_elem_number


% --- Executes on button press in checkbox_boundary_conditions.
function checkbox_boundary_conditions_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_boundary_conditions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_boundary_conditions


% --- Executes on button press in checkbox_external_forces.
function checkbox_external_forces_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_external_forces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_external_forces


% --- Executes on button press in checkbox_points_of_cracks.
function checkbox_points_of_cracks_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_points_of_cracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_points_of_cracks


% --- Executes on button press in checkbox20.
function checkbox20_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox20


% --- Executes on button press in checkbox22.
function checkbox22_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox22


% --- Executes on button press in checkbox21.
function checkbox21_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox21


% --- Executes on button press in checkbox19.
function checkbox19_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox19


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_real_shaped_crack.
function checkbox_real_shaped_crack_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_real_shaped_crack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_real_shaped_crack


% --------------------------------------------------------------------
function Untitled_help_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_About_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox25.
function checkbox25_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox25


% --- Executes on button press in checkbox26.
function checkbox26_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox26


% --- Executes on button press in checkbox27.
function checkbox27_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox27


% --- Executes on button press in checkbox28.
function checkbox28_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox28


% --- Executes on button press in checkbox_plot_KIKII.
function checkbox_plot_KIKII_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_KIKII (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_KIKII



function edit_crack_num_to_plot_KIKII_Callback(hObject, eventdata, handles)
% hObject    handle to edit_crack_num_to_plot_KIKII (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_crack_num_to_plot_KIKII as text
%        str2double(get(hObject,'String')) returns contents of edit_crack_num_to_plot_KIKII as a double


% --- Executes during object creation, after setting all properties.
function edit_crack_num_to_plot_KIKII_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_crack_num_to_plot_KIKII (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tip_num_to_plot_KIKII_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tip_num_to_plot_KIKII (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tip_num_to_plot_KIKII as text
%        str2double(get(hObject,'String')) returns contents of edit_tip_num_to_plot_KIKII as a double


% --- Executes during object creation, after setting all properties.
function edit_tip_num_to_plot_KIKII_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tip_num_to_plot_KIKII (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% �ر�ȫ��figures
function pushbutton_close_all_fig_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close_all_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%close all
delete(findall(0,'type','figure'))
