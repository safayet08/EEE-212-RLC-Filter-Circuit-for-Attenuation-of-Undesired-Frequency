function varargout = project(varargin)
    % PROJECT MATLAB code for project.fig
    %      PROJECT, by itself, creates a new PROJECT or raises the existing
    %      singleton*.
    %
    %      H = PROJECT returns the handle to a new PROJECT or the handle to
    %      the existing singleton*.
    %
    %      PROJECT('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in PROJECT.M with the given input arguments.
    %
    %      PROJECT('Property','Value',...) creates a new PROJECT or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before project_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to project_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help project

    % Last Modified by GUIDE v2.5 08-Sep-2019 16:12:29

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @project_OpeningFcn, ...
                       'gui_OutputFcn',  @project_OutputFcn, ...
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


% --- Executes just before project is made visible.
function project_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to project (see VARARGIN)

    % Choose default command line output for project
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes project wait for user response (see UIRESUME)
    % uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = project_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;


function WC2_Callback(hObject, eventdata, handles)
    % hObject    handle to WC2 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of WC2 as text
    %        str2double(get(hObject,'String')) returns contents of WC2 as a double


% --- Executes during object creation, after setting all properties.
function WC2_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to WC2 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on selection change in type.
function type_Callback(hObject, eventdata, handles)
    % hObject    handle to type (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns type contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from type


% --- Executes during object creation, after setting all properties.
function type_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to type (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: listbox controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end



function WC1_Callback(hObject, eventdata, handles)
    % hObject    handle to WC1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of WC1 as text
    %        str2double(get(hObject,'String')) returns contents of WC1 as a double

    

function baseF_Callback(hObject, eventdata, handles)
    % hObject    handle to baseF (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of baseF as text
    %        str2double(get(hObject,'String')) returns contents of baseF as a double


% --- Executes during object creation, after setting all properties.
function baseF_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to baseF (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

% --- Executes during object creation, after setting all properties.
function WC1_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to WC1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


    % --- Executes on button press in pushbutton1.
    function pushbutton1_Callback(hObject, eventdata, handles)
    % hObject    handle to pushbutton1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    WC2 = str2double(get(handles.WC1, 'string'));
    WC1 = str2double(get(handles.WC2, 'string'));
    base = str2double(get(handles.baseF, 'string'));
    listBoxItems = handles.type.String;
    listBoxSelectedIndexes = handles.type.Value;
    selection = listBoxItems{listBoxSelectedIndexes};
    
    axes(handles.ckt); hold off;
    axes(handles.DB); hold off;
    axes(handles.TF); hold off;
    axes(handles.Response); hold off;
    
    if strcmp(selection, 'Band Pass Filter')
        band_pass(WC1, WC2, hObject, eventdata, handles, base);
    elseif strcmp(selection, 'Band Reject Filter')
        band_stop(WC1, WC2, hObject, eventdata, handles, base);
    elseif strcmp(selection, 'High Pass Filter')
        high_pass(WC1, hObject, eventdata, handles, base);
    elseif strcmp(selection, 'Low Pass Filter')
        low_pass(WC2, hObject, eventdata, handles, base);
    end
    
   
    

function pos = normalize(x, idx) 
    pos = [];
    pos = idx(1);
    for i = 2: 1 : length(idx)
        if idx(i) - idx(i - 1) > 1
            pos = [pos idx(i)];
        end
    end

function [idx] = cutoff(r, value, tol)
    v = abs(r - (value));
    d = v .* ones(size(r));
    minDev = min(d);
    eps = tol;
    idx = find(abs(d - minDev) <= eps);


function band_pass(WC1, WC2, hObject, eventdata, handles, base) 
    L = 2.7 * 10 ^ -3;
    BW = abs(WC2 - WC1);
    R = BW * L;
    W0 = sqrt(WC1 * WC2);
    Q = W0 / BW;
    C = 1 / (W0^2 * L);
    
    R_label = sprintf('R = %f', R);
    L_label = sprintf('L = %f', L);
    C_label = sprintf('C = %f', C);
    Q_label = sprintf('Quality Factor = %f', Q);
    set(handles.RR, 'String', R_label);
    set(handles.LL, 'String', L_label);
    set(handles.CC, 'String', C_label);
    set(handles.QQ, 'String', Q_label);
    
    w = linspace(max(0, -BW + WC1), WC2 + 2*BW, 10000);
    TF =@(w, R, L, C) (1i*R*C .* w) ./ (1i*R*C .* w - w .* w .* L*C + 1); %Across R
    
    r = abs(TF(w, R, L, C));
    
    axes(handles.ckt);
    im = imread('BAND_PASS.jpg');
    imshow(im);
    
    axes(handles.DB);
    DB = 20 .* log10(r);
    plot(w, DB, 'LineWidth', 1.5); hold on; grid on;
    
    [p, idx0] = findpeaks(DB);
    plot(w(idx0), DB(idx0), 'ro', 'LineWidth', 3); hold on; %PEAK POINT PLOTTING
    
    idxC = cutoff(DB, DB(idx0(1)) - 3, 2 * 10 ^ -3);
    idxC = normalize(w, idxC);
    
    DBWC1 = w(idxC(1));
    DBWC2 = w(idxC(end));
    DBWC1_label = sprintf('Cutoff = %f', DBWC1);
    DBWC2_label = sprintf('Cutoff = %f', DBWC2);
    set(handles.DBWC1, 'String', DBWC1_label);
    set(handles.DBWC2, 'String', DBWC2_label);
    
    plot(w(idxC(1:end)), DB(idxC(1:end)), 'ro', 'LineWidth', 3); %CUTOOF POINT PLOTTING
    
    xlabel("---------------> Angular Frequency(\omega)")
    ylabel("Gain in DB")
    title("Plot of gain in DB Across R")
    
    axes(handles.TF);
    plot(w, r, 'LineWidth', 1.5); hold on; grid on; %PLOTTING THE TRANSFER FUNCTION
    title("Plot of the transfer function");
    xlabel("---------------> Angular Frequency(\omega)")
    
    [p, idx0] = findpeaks(r);
    plot(w(idx0), r(idx0), 'ro', 'LineWidth', 3); hold on; %PEAK POINT PLOTTING
    
    idxC = cutoff(r, r(idx0(1)) / sqrt(2), 2 * 10 ^ -3);
    idxC = normalize(w, idxC);
    
    TFWC1 = w(idxC(1));
    TFWC2 = w(idxC(end));
    TFWC1_label = sprintf('Cutoff = %f', TFWC1);
    TFWC2_label = sprintf('Cutoff = %f', TFWC2);
    set(handles.TFWC1, 'String', TFWC1_label);
    set(handles.TFWC2, 'String', TFWC2_label);
    
    plot(w(idxC(1:end)), r(idxC(1:end)), 'ro', 'LineWidth', 3); %CUTOOF POINT PLOTTING
    
    band_pass_simulation(R, L, C, TF, WC1, WC2, hObject, eventdata, handles, base);

function band_pass_simulation(R, L, C, TF, WC1, WC2, hObject, eventdata, handles, base) 
    sinusoid =@(w, t, phi) sin(w * t + phi);
    
%   For a series of wt, 2wt, 3wt, 4wt, ....., nwt; set base here
    upto = 10;
    
    t = linspace(0, 2.5 * (2 * pi) / (base), 10000);
    
    total_response = zeros(size(t));
    for i = 1 : upto
        w = base * i;
        v_out = TF(w, R, L, C);
        response = abs(v_out) * sinusoid(w, t, 0);
        total_response = total_response + response;
    end
    
    pass_band_response = zeros(size(t));
    for i = 1 : upto
        w = base * i;
        if w >= WC1 && w <= WC2
            v_out = TF(w, R, L, C);
            response = abs(v_out) * sinusoid(w, t, 0);
            pass_band_response = pass_band_response + response;
        end
    end
    
    axes(handles.Response);
    plot(t, total_response, 'LineWidth', 1.5); hold on; grid on;
    plot(t, pass_band_response, 'LineWidth', 1.5); hold on;
    legend("Total Response", "Pass Band Response");
    title("Plot of Total Response and Pass Band Response");
    xlabel("---------------> Time (t)");
    ylabel("---------------> Voltage across R (V_{out})");
   
function band_stop(WC1, WC2, hObject, eventdata, handles, base) 

    L = 1 * 10 ^ -6;
    BW = abs(WC2 - WC1);
    R = BW * L;
    W0 = sqrt(WC1 * WC2);
    Q = W0 / BW;
    C = 1 / (W0^2 * L);
    TF =@(w, R, L, C) (w .* w .* L*C + 1) ./ (1i*R*C .* w - w .* w .* L*C + 1); %ACROSS LC
    
    R_label = sprintf('R = %f', R);
    L_label = sprintf('L = %f', L);
    C_label = sprintf('C = %f', C);
    Q_label = sprintf('Quality Factor = %f', Q);
    set(handles.RR, 'String', R_label);
    set(handles.LL, 'String', L_label);
    set(handles.CC, 'String', C_label);
    set(handles.QQ, 'String', Q_label);
    
    w = linspace(max(-BW + WC1, 0), WC2 + BW, 10000);
    ff =@(w, R, L, C) (1i*R*C .* w) ./ (1i*R*C .* w - w .* w .* L*C + 1);
    
    r = abs(ff(w, R, L, C));
    
    axes(handles.ckt);
    im = imread('BAND_STOP.jpg');
    imshow(im);
     
    axes(handles.DB);
    DB = 20 .* log10(abs(r));
    plot(w, -DB, 'LineWidth', 1.5); hold on; grid on;
    
    [p, idx0] = findpeaks(DB);
    plot(w(idx0), -DB(idx0), 'ro', 'LineWidth', 3); hold on; %PEAK POINT PLOTTING
    
    idxC = cutoff(DB, DB(idx0(1)) - 3, 2 * 10 ^ -3);
    idxC = normalize(w, idxC);
    
    DBWC1 = w(idxC(1));
    DBWC2 = w(idxC(end));
    DBWC1_label = sprintf('Cutoff = %f', DBWC1);
    DBWC2_label = sprintf('Cutoff = %f', DBWC2);
    set(handles.DBWC1, 'String', DBWC1_label);
    set(handles.DBWC2, 'String', DBWC2_label);
    
    plot(w(idxC(1:end)), -DB(idxC(1:end)), 'ro', 'LineWidth', 3); %CUTOOF POINT PLOTTING
    
    xlabel("---------------> Angular Frequency(\omega)")
    ylabel("Gain in DB")
    title("Plot of gain in DB Across LC")
   
    axes(handles.TF);
    plot(w, -r, 'LineWidth', 1.5); hold on; %PLOTTING THE TRANSFER FUNCTION
    title("Plot of the transfer function");
    xlabel("---------------> Angular Frequency(\omega)")
    
    [p, idx0] = findpeaks(r); 

    plot(w(idx0), -r(idx0(1)), 'ro', 'LineWidth', 3); hold on; %PEAK POINT PLOTTING
    
    idxC = cutoff(r, r(idx0(1)) / sqrt(2), 2 * 10 ^ -3);
    q = normalize(w, idxC);
    idxC = q;
    
    TFWC1 = w(idxC(1));
    TFWC2 = w(idxC(end));
    TFWC1_label = sprintf('Cutoff = %f', TFWC1);
    TFWC2_label = sprintf('Cutoff = %f', TFWC2);
    set(handles.TFWC1, 'String', TFWC1_label);
    set(handles.TFWC2, 'String', TFWC2_label);
    
    plot(w(idxC(1:end)), -r(idxC(1:end)), 'ro', 'LineWidth', 3); %CUT OFF POINT PLOTTING
    grid on;
    
    band_stop_simulation(R, L, C, TF, WC1, WC2, hObject, eventdata, handles, base);

    
function band_stop_simulation(R, L, C, TF, WC1, WC2, hObject, eventdata, handles, base) 
    sinusoid =@(w, t, phi) sin(w * t + phi);
    
%   For a series of wt, 2wt, 3wt, 4wt, ....., nwt; set base here
    upto = 5;
    
    t = linspace(0, 3 * (2 * pi) / (base), 10000);
    
    total_response = zeros(size(t));
    for i = 1 : upto
        w = base * i;
        v_out = -TF(w, R, L, C);
        response = abs(v_out) * sinusoid(w, t, 0);
        total_response = total_response + response;
    end
    
    WC1
    WC2
    
    pass_band_response = zeros(size(t));
    for i = 1 : upto
        w = base * i;
        if w >= WC2 || w <= WC1
            w
            v_out = -TF(w, R, L, C);
            response = abs(v_out) * sinusoid(w, t, 0);
            pass_band_response = pass_band_response + response;
        end
    end
    
    axes(handles.Response);
    plot(t, total_response, 'LineWidth', 1.5); hold on; grid on;
    plot(t, pass_band_response, 'LineWidth', 1.5); hold on;
    legend("Total Response", "Pass Band Response");
    title("Plot of Total Response and Pass Band Response");
    xlabel("---------------> Time (t)");
    ylabel("---------------> Voltage across LC (V_{out})");
    
function high_pass(WC,  hObject, eventdata, handles, base) %RC circuit

    C = 10 ^ -6;
    R = 1 / (WC * C);
    
    Q = 1 / (WC * R * C);
    w = linspace(0, 9 * WC, 10000);
    TF =@(w, R, C) (1i * R * C .* w) ./ ((1i * R * C .* w) + 1) ; %Across R
    r = abs(TF(w, R, C));
    
    R_label = sprintf('R = %f', R);
    L_label = sprintf('L = X');
    C_label = sprintf('C = %f', C);
    Q_label = sprintf('Quality Factor = %f', Q);
    set(handles.RR, 'String', R_label);
    set(handles.LL, 'String', L_label);
    set(handles.CC, 'String', C_label);
    set(handles.QQ, 'String', Q_label);
    
    
    axes(handles.ckt);
    im = imread('RC_HIGH_PASS.jpg');
    imshow(im);
    
    axes(handles.DB);
    DB = 20 .* log10(r);
    plot(w, DB, 'LineWidth', 1.5); hold on; grid on;
    
    title("Plot of gain in DB");
    xlabel("---------------> Angular Frequency(\omega)");
    ylabel("---------------> Gain in DB Across R");
    peak = max(DB);
    idxC = cutoff(DB, peak - 3, 2 * 10 ^ -3);
    idxC = normalize(w, idxC);
    plot(w(idxC(1:end)), DB(idxC(1:end)), 'ro', 'LineWidth', 3); %CUTOOF POINT PLOTTING
    
    DBWC1 = w(idxC(1));
    DBWC2 = w(idxC(end));
    DBWC1_label = sprintf('Cutoff = %f', DBWC1);
    DBWC2_label = sprintf('Cutoff = %f', DBWC2);
    set(handles.DBWC1, 'String', DBWC1_label);
    set(handles.DBWC2, 'String', DBWC2_label);
    
    axes(handles.TF);
    plot(w, r, 'LineWidth', 1.5); hold on; %PLOTTING THE TRANSFER FUNCTION
    grid on;
    title("Plot of the transfer function");
    xlabel("---------------> Angular Frequency(\omega)")
    ylabel("---------------> The Transfer Function Across R")
    
    idxC = cutoff(r, r(length(r)) / sqrt(2), 2 * 10 ^ -3);
    idxC = normalize(w, idxC);
    
    TFWC1 = w(idxC(1));
    TFWC2 = w(idxC(end));
    TFWC1_label = sprintf('Cutoff = %f', TFWC1);
    TFWC2_label = sprintf('Cutoff = %f', TFWC2);
    set(handles.TFWC1, 'String', TFWC1_label);
    set(handles.TFWC2, 'String', TFWC2_label);
    
    plot(w(idxC(1:end)), r(idxC(1:end)), 'ro', 'LineWidth', 3); %CUTOOF POINT PLOTTING
    hold on;
    
    high_pass_simulation(R, C, TF, WC,  hObject, eventdata, handles, base);
    
function high_pass_simulation(R, C, TF, WC,  hObject, eventdata, handles, base) 
    sinusoid =@(w, t, phi) sin(w * t + phi);
    
%   For a series of wt, 2wt, 3wt, 4wt, ....., nwt; set base here
    % base = 1000;
    upto = 5;
    
    t = linspace(0, 2 * (2 * pi) / (base), 10000);
    
    total_response = zeros(size(t));
    for i = 1 : upto
        w = base * i;
        v_out = TF(w, R, C);
        response = abs(v_out) * sinusoid(w, t, 0);
        total_response = total_response + response;
    end
    
    pass_band_response = zeros(size(t));
    for i = 1 : upto
        w = base * i;
        if w >= WC
            v_out = TF(w, R, C);
            response = abs(v_out) * sinusoid(w, t, 0); 
            pass_band_response = pass_band_response + response;
        end
    end
    
    
    axes(handles.Response);
    plot(t, total_response, 'LineWidth', 1.5); hold on; grid on;
    plot(t, pass_band_response, 'LineWidth', 1.5); hold on;
    legend("Total Response", "Pass Band Response");
    title("Plot of Total Response and Pass Band Response");
    xlabel("---------------> Time (t)");
    ylabel("---------------> Voltage across R (V_{out})");

function low_pass(WC, hObject, eventdata, handles, base) %RL circuit
    L = 2.7 * 10 ^ -3;
    R = WC * L;
    Q = WC * L / R;

    R_label = sprintf('R = %f', R);
    L_label = sprintf('L = %f', L);
    C_label = sprintf('C = X');
    Q_label = sprintf('Quality Factor = %f', Q);
    set(handles.RR, 'String', R_label);
    set(handles.LL, 'String', L_label);
    set(handles.CC, 'String', C_label);
    set(handles.QQ, 'String', Q_label);
    
    w = linspace(0, 9 * WC, 10000);
    TF =@(w, R, L) (R) ./ (1i*L.* w + R); %Across L
    r = abs(TF(w, R, L));
    
    axes(handles.ckt);
    im = imread('RL_HIGH_PASS.jpg');
    imshow(im);
    
        
    axes(handles.DB);
    DB = 20 .* log10(r);
    plot(w, DB, 'LineWidth', 1.5); hold on; grid on;
    
    title("Plot of gain in DB");
    xlabel("---------------> Angular Frequency(\omega)");
    ylabel("---------------> Gain in DB Across L");
    peak = max(DB);
    idxC = cutoff(DB, peak - 3, 2 * 10 ^ -3);
    idxC = normalize(w, idxC);
    plot(w(idxC(1:end)), DB(idxC(1:end)), 'ro', 'LineWidth', 3); %CUTOOF POINT PLOTTING
    
    DBWC1 = w(idxC(1));
    DBWC2 = w(idxC(end));
    DBWC1_label = sprintf('Cutoff = %f', DBWC1);
    DBWC2_label = sprintf('Cutoff = %f', DBWC2);
    set(handles.DBWC1, 'String', DBWC1_label);
    set(handles.DBWC2, 'String', DBWC2_label);
    
    axes(handles.TF);
    plot(w, r, 'LineWidth', 1.5); hold on; %PLOTTING THE TRANSFER FUNCTION
    grid on;
    title("Plot of the transfer function");
    xlabel("---------------> Angular Frequency(\omega)")
    ylabel("---------------> The Transfer Function Across L")

    idxC = cutoff(r, 1 / sqrt(2), 2 * 10 ^ -3);
    idxC = normalize(w, idxC);
    
    TFWC1 = w(idxC(1));
    TFWC2 = w(idxC(end));
    TFWC1_label = sprintf('Cutoff = %f', TFWC1);
    TFWC2_label = sprintf('Cutoff = %f', TFWC2);
    set(handles.TFWC1, 'String', TFWC1_label);
    set(handles.TFWC2, 'String', TFWC2_label);
    
    plot(w(idxC(1:end)), r(idxC(1:end)), 'ro', 'LineWidth', 3); %CUTOOF POINT PLOTTING
    hold on;
    low_pass_simulation(R, L, TF, WC, hObject, eventdata, handles, base);
    
function low_pass_simulation(R, L, TF, WC, hObject, eventdata, handles, base) 
    sinusoid =@(w, t, phi) sin(w * t + phi);
    
%   For a series of wt, 2wt, 3wt, 4wt, ....., nwt; set base here
    upto = 5;
    
    t = linspace(0, 6 * (2 * pi) / (base), 10000);
    
    total_response = zeros(size(t));
    for i = 1 : upto
        w = base * i;
        v_out = TF(w, R, L);
        response = abs(v_out) * sinusoid(w, t, 0);
        total_response = total_response + response;
    end
    
    pass_band_response = zeros(size(t));
    for i = 1 : upto
        w = base * i;
        if w <= WC
            w
            v_out = TF(w, R, L);
            response = abs(v_out) * sinusoid(w, t, 0);
            pass_band_response = pass_band_response + response;
        end
    end
    
    axes(handles.Response);
    plot(t, total_response, 'LineWidth', 1.5); hold on; grid on;
    plot(t, pass_band_response, 'LineWidth', 1.5); hold on;
    legend("Total Response", "Pass Band Response");
    title("Plot of Total Response and Pass Band Response");
    xlabel("---------------> Time (t)");
    ylabel("---------------> Voltage across R (V_{out})");
