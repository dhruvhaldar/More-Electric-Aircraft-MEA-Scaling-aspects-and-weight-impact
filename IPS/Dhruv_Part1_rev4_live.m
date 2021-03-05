%% 
% 
% Namn: Haldar, Dhruv
% 
% Datum: 20200804
%% *Airfoil Generation in MATLAB using XFOIL*
% |Drela M. (1989) XFOIL: An Analysis and Design System 
%for Low Reynolds Number 
% Airfoils. In: Mueller T.J. (eds) Low Reynolds Number Aerodynamics. 
%Lecture Notes in Engineering, vol 54.
%Springer, Berlin, Heidelberg. 
%https://doi.org/10.1007/978-3-642-84010-4_1|
% |Drela, M., & Youngren, H. (2001). 
% _XFOIL 6.9 User Primer_. Retrieved from 
% MIT: https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt

%% Clear Workspace and Program Memory
clear;
clc;
% Download and Unzip XFoil and Gmsh
% The following automatically downloads the XFoil 
%file dhruvhaldar/XFoil-cross 
% and Gmsh file from dhruvhaldar/Gmsh-cross Github repository 
%into the user folder.If the XFoil
% and Gmsh directory already exist 
%the code proceeds to the next section.

sevenzip_path=pwd;
while isfile("Xfoil-cross.zip")==0 ...
&&isfile("xfoil.exe")==0&&isfolder("Xfoil-cross")==0
    disp("Downloading XFoil from dhruvhaldar ...
/XFoil-cross Github Repository");
    unzip('https://github.com/dhruvhaldar/XFoil-Cross ...
/raw/main/Xfoil-cross.zip');
    addpath(genpath('Xfoil-cross'));
    break    
end
while isfile("Gmsh-cross.zip")==0 ...
&&isfile("gmsh.exe")==0&&isfolder("Gmsh-cross")==0
    disp("Downloading Gmsh from dhruvhaldar ...
/Gmsh-cross Github Repository, Please wait...");
    unzip('https://github.com/dhruvhaldar/Gmsh-cross ...
/raw/main/Gmsh-cross.zip');
    addpath(genpath('Gmsh-cross'));
    break
end
% Download and Unzip UIUC Airfoil Database (For reference)

while isfolder("coord_seligFmt")==0
    disp("Downloading Airfoil files from UIUC Airfoil ...
Database, Please wait...");
    unzip('https://m-selig.ae.illinois.edu/ads ...
/archives/coord_seligFmt.zip');
    addpath(genpath('coord_seligFmt'));
    break
end
%% *Type of Airfoil Geometry Input*
%% 
% * *Generate standard airfoil for NACA 4 and 5-digit airfoil*
% * *Generate/Load custom airfoil 
% for NACA 6,7 and 8-digit airfoil series*
% *Generate NACA 4-digit or 5-digit airfoil type*
% We can simply change the NACA airfoil code by
% changing this line, 
%XFoil will 
% understand it automatically and obtain data
% from its stored source code file. 
% The airfoil data files are obtained at the 
%UIUC Database. In order to fetch/refresh 
% the directory we would need Python coding which 
%is beyond the scope of MATLAB. 
% Thus, it is possible to manually 
%download a dat file if we want to input to 
% XFoil.

NACA = '24012';
if isempty(NACA)
    NACA=inputdlg('Input NACA 4,5,6,7 and 8-digit Airfoil name');
end
strlen_NACA=strlength(NACA);
if strlen_NACA<=3 || strlen_NACA>=8
    error("Invalid NACA Input")
end
% Generation for NACA 6,7 and 8-digit Airfoil Series
% By default, XFoil database does not have the 
%airfoil data for 6,7 and 8 Airfoil 
% series as the original XFoil was developed 
%in the year 2000 when these types 
% of airfoils did not exist. Thus, the airfoil 
%data is to be fetched and loaded 
% from separately created repository.The XFoil 
%Input data file parameters are 
% changed accordingly.
% Angle of attack in degrees (°)

AoA = '3';       
% Number of Panel nodes for coarseness of airfoil geometry

nodes = '200'; 
% Reynold's Number

Re = 1e6;
% Mach Number

Mach = 0.2;
if Mach>1
    error("Local Supersonic flow not allowed in Program");
elseif Mach<0
    error ("Invalid input");
end
%% 
% Input Number of Panel Nodes

if isempty(Re)
    nodes=inputdlg('Input Number of Panel Nodes');
end
%% 
% Input Bunching Parameters of Panel

bunching=1;
if isempty(bunching)
   bunching=inputdlg('Input Bunching Parameters of Panel'); 
end
% Atmospheric Input Conditions

if isempty(Re)  
    v=inputdlg(sprintf(['Enter Velocity of Fluid  ...
in m/s=%g\n' ...
        'Enter Characteristic Length or Chord Width  ...
of Airfoil in m=%g\n' ...
        'Enter Kinematic Viscosity in m3/s=%g\n'],v,l,nu));
    v=str2double(v);
    l=str2double(l);
    nu=str2double(nu);
    Re = (v*l)/nu;
    Re_answer=msgbox(sprintf('Reynolds Number is=%g' ...
,Re),'Reynolds Number','modal');
end
% Mach Number Input

if isempty(Mach)
    Mach=inputdlg('Enter Mach number in m/s');
    if str2double(Mach) 
        error("Supersonic flow not allowed in Program");
    end
end
% Airfoil Interpolation
% XFoil also offers the linear interpolation of two 
%different airfoils, essentially 
% required for complex Airfoil geometries thus
%effectively creating a mixed-airfoil 
% shape.It allows interpolating or "blending" of 
%airfoils in various proportions 
% from cubic splines;
% Airfoil coordinates filename

AirfoilVarDat = 'Save_Airfoil.dat';
% Pressure coefficient filename
% Pressure distribution in Aerodynamics is given as $C_p$

CpVarDat = 'Save_Cp.dat';
% Polar coordinates filename

Polvardat= 'Save_Pol.dat';
% Leading edge computation filename

Leadedgedat='Save_Leadedge.dat';
%% Generate Airfoil from XFoil
% Create the XFoil input data file

fid = fopen('xfoil_input.dat','w+');
if fid==-1
    error ('Error Opening Xfoil Input File');
end
if strlen_NACA>3 && strlen_NACA<=5
    fprintf(fid,['NACA ' NACA '\n']);
    fprintf(fid,'PPAR\n');
    fprintf(fid,['N ' nodes '\n']);
    fprintf(fid,'\n\n');
elseif strlen_NACA>5 && strlen_NACA<=8
    for i=1:1
        if strcmpi(NACA,"ste87151")
            fprintf(fid,['LOAD ' 'nlf415.dat' '\n']);
        clear i s;
        break;
        elseif strcmpi(NACA,"ste87391")
            fprintf(fid,['LOAD ' 'nlf1015.dat' '\n']);
        clear i s;
        break;
        elseif strcmpi(NACA,"bacj")
            fprintf(fid,['LOAD ' 'bacj.dat' '\n']);
        clear i s;
        break;
        elseif strcmpi(NACA,"b737")
            fprintf(fid,['LOAD ' 'b737a.dat' '\n']);
        clear i s;
        break;
        elseif strcmpi(NACA,"b707")
            fprintf(fid,['LOAD ' 'b707a.dat' '\n']);
        clear i s;
        break;
        end
     end
    fprintf(fid,['LOAD ' 'naca' NACA '.dat' '\n']);
    fprintf(fid,'PPAR\n');
    fprintf(fid,['N ' nodes '\n']);
    fprintf(fid,['P' bunching '\n']);
    fprintf(fid,'\n\n');
end
% Save the airfoil data points

fprintf(fid,['PSAV ' AirfoilVarDat '\n']);
% Compute Polar Coordinates
% Polar coodinates are used in the  
%calculation of Lift and Drag Coefficients.

fprintf(fid,'OPER\n');
fprintf(fid,'PACC 1 \n');
fprintf(fid,'\n\n');
fprintf(fid,['ALFA ' AoA '\n']);
fprintf(fid,['CPWR ' CpVarDat '\n']);
fprintf(fid,'PWRT\n');
fprintf(fid,[Polvardat '\n']);
while isfile(Polvardat)                                        
    fprintf(fid,'y \n');
    break
end

fprintf(fid,'\n');
fprintf(fid,'GDES\n');
fprintf(fid,'LERA\n');
fprintf(fid,'\n\n');
% Close File

fclose(fid);
%% Call XFoil
% Run XFoil using XFoil input data file 
% |Referenced from Louis Edelman (2020). 
%Xfoil Interface Updated
%(https://www.mathworks.com/ ...
%matlabcentral/fileexchange/49706-xfoil-interface-updated), 
% MATLAB Central File Exchange. Retrieved August 9, 2020.|
% Run XFoil

while ispc == 1 
    copyfile xfoil_input.dat Xfoil-cross\ f;
    cd Xfoil-cross;
    cmd_xfoil = 'xfoil < xfoil_input.dat';
    [status_xfoilwin,result_xfoilwin] = system(cmd_xfoil);
        while (status_xfoilwin<=0)
            disp(result_xfoilwin);
            disp("XFoil is unable to load input file");
            error(cmd_xfoil);
        end
    break
end
%% Read the generated Airfoil Data File
% Open Airfoil data file for reading

Airfoil_file_index = fopen(AirfoilVarDat);                                 
% Read data from Airfoil data file

datbuffer_airfoil = textscan(Airfoil_file_index, ...
'%f %f','CollectOutput',1,...     
                                 'Delimiter','','HeaderLines',0);
% Close Airfoil data file

fclose(Airfoil_file_index);   
fclose('all');
% Separate X and Y coordinate data points
% X-coordinate data buffer

XB = datbuffer_airfoil{1}(:,1);  
% Y-coordinate data buffer

YB = datbuffer_airfoil{1}(:,2);
% Read Pressure Coefficient .dat file
% Open $C_p$ file for reading

CP_file_index = fopen(CpVarDat);
% Ready data from $C_p$ file

while ispc
    datbuffer_cp = textscan ...
(CP_file_index,'%f %f %f','HeaderLines',3,...
                            'CollectOutput',1,...
                            'Delimiter','');
    break
end
while isunix
    datbuffer_cp = textscan ...
(CP_file_index,'%f %f','HeaderLines',1,...         
                            'CollectOutput',1,...
                            'Delimiter','');
    datbuffer_cp{1,1}(:,3) = datbuffer_cp{1,1}(:,2);                    
    datbuffer_cp{1,1}(:,2) = datbuffer_airfoil{1,1}(:,2);
    break
end
%  Close $C_p \;$file 

fclose(CP_file_index);                                                  
% Separate $C_p$ data from $C_p \;$file  
% X-coordinate

X_0  = datbuffer_cp{1,1}(:,1);   
%% 
% Y-coordinate

Y_0  = datbuffer_cp{1,1}(:,2);     
%% 
% Pressure Coefficient data

Cp_0 = datbuffer_cp{1,1}(:,3); 
% Read Polar Coordinate input File
% Polar Coordinate input Filename

Polvardat= 'Save_Pol.dat';
% Open Polar Coordinate input file for reading

Pol_file_index = fopen(Polvardat);
% Ready data from Polar Coordinate input file

poldatbuffer = textscan ...
(Pol_file_index,'%f %f %f %f %f %f %f' ...
,'CollectOutput',1,...
                                 'Delimiter','','HeaderLines',12);
% Close Polar Coordinate input file

fclose(Pol_file_index);
% Separate Polar data from Polar coordinate input file

CL = poldatbuffer{1,1}(2);
CD = poldatbuffer{1,1}(3);
CM = poldatbuffer{1,1}(5);
% Obtain xdata, ydata and pressure coefficients of Airfoil
% Split airfoil into (U)pper and (L)ower

XB_U = XB(YB >= 0);
XB_L = XB(YB < 0);
YB_U = YB(YB >= 0);
YB_L = YB(YB < 0);
% Split Xfoil results into (U)pper and (L)ower

Cp_U = Cp_0(YB >= 0);
Cp_L = Cp_0(YB < 0);
X_U  = X_0(YB >= 0);
X_L  = X_0(YB < 0);
% Plot Airfoil

figure(1);
cla; hold on; grid off;
set(gcf,'Color','White');
set(gca,'FontSize',12);
plot(XB_U,YB_U,'b.-');
plot(XB_L,YB_L,'r.-');
xlabel('X Coordinate');
ylabel('Y Coordinate');
axis equal;
%  Plot: Pressure coefficient

figure(2);
cla; hold on; grid on;
set(gcf,'Color','White');
set(gca,'FontSize',12);
plot(X_U,Cp_U,'bo-','LineWidth',2);
plot(X_L,Cp_L,'ro-','LineWidth',2);
xlabel('X Coordinate');
ylabel('Cp');
ylim('auto');
set(gca,'Ydir','reverse')
% Display Lift,Drag and Moment Coefficients from Polar computations

disp((sprintf('Lift Coefficient= %
g\nDrag Coefficient= %g\nMoment Coefficient= %g',CL,CD,CM)));
% Gmsh
% Change to Gmsh directory

while isfile("xfoil.exe")==1
    cd (sevenzip_path);
    cd Gmsh-cross;
    GmshPath=pwd;
    break
end
% Figure to Gmsh Geometry
% Input Parameters 

meshsize=0.5;
chordlen=1;
gridPtsU=size(XB_U,1);
gridPtsL=size(XB_L,1);
gridPts=gridPtsL+gridPtsU;
while gridPtsL<3||gridPtsU<3
    error("Enter specified number of grid points (minimum 3)");
end
saveFilename_geo=append('NACA',NACA,'_gmsh','.geo');   
%% 
% Multiply airfoil coordinates by the chord length

XB_L_chord = XB_L*chordlen;                                                          % Multiply lower X data by chord
YB_L_chord = YB_L*chordlen;                                                          % Multiply lower Y data by chord
XB_U_chord = XB_U*chordlen;                                                          % Multiply upper X data by chord 
YB_U_chord = YB_U*chordlen;                                                          % Multiply upper Y data by chord
%% 
% Creation of array for writing  ...
%consisting of upper and lower data points

pointU = zeros(gridPtsU-1,3);
pointL = zeros(gridPtsL-1,3);
indU = 1;
    for i = 1:1:gridPtsU
        pointU(indU,1) = XB_U_chord(i);
        pointU(indU,2) = YB_U_chord(i);
        pointU(indU,3) = 0;
        pointU(indU,4) = meshsize;
        indU = indU + 1;
    end
    indL = 1;
    for i = 1:1:gridPtsL
        pointL(indL,1) = XB_L_chord(i);
        pointL(indL,2) = YB_L_chord(i);
        pointL(indL,3) = 0;
        pointL(indL,4) = meshsize;
        indL = indL + 1;
    end
%% 
% Determine number of points and lines

numPtsU = length(pointU(:,1));
numPtsL = length(pointL(:,1));
numPts=numPtsL+numPtsU;
numLns = numPts;
%% 
% Write data to the .geo file

fid = fopen(saveFilename_geo,'w');
fprintf(fid,'/*Gmsh Airfoil Data ...
 file generated from MATLAB script\n');
fprintf(fid,'@author Dhruv Haldar*/\n');
fprintf(fid,'\n/*Airfoil Data Points */\n');
fprintf(fid,'\r/*Airfoil Upper Surface Data Points*/\n');
    for i = 1:1:numPtsU
        fprintf(fid,'Point(%i) =  ...
{%g, %g, %g, %g};\r\n' ...
,i,pointU(i,1),pointU(i,2),...
                                  pointU(i,3),pointU(i,4));
    end
fprintf(fid,'\r\n/*Airfoil Lower Surface Data Points*/\n');
    for i = 1:1:numPtsL
        j=i+numPtsU;
        fprintf(fid,'Point(%i) =  ...
{%g, %g, %g, %g};\r\n' ...
,j,pointL(i,1),pointL(i,2),...
                                   pointL(i,3),pointL(i,4));                                    
    end
fprintf(fid,'\r/*Define Airfoil Data lines*/\n');
     for i = 1:1:1   
        fprintf(fid,'Spline(%i) = {%i, %i};\r\n',j+1,1,numLns);
        break
    end
fprintf(fid,'\r\n/*Define Trailing edge*/\n');
    for i = 1:1:1
        fprintf(fid,'Spline(%i) = {%i,%i};\r\n',j+2,numLns,1);
    end
% Gmsh Enclosure
% A Bullet-shaped domain is constructed for the airfoil.
% 
% Determine number of lines and points

encPts=6;
encLns = encPts;
%% 
% Determine coordinates of point G 
&from the first written value in the upper 
% layer coordinate data array

pointG = zeros(1,4);
for k=1:1
   pointG(1,1)=pointU(1,1);
   pointG(1,2)=pointU(1,2);
   pointG(1,3)=pointU(1,3);
   pointG(1,4)=pointU(1,4);
   break
end
%% 
% Determination of Horizontal Distance for  ...
%Enclosure coordinate calculation

horizontal_distance=max(pointU(:,1))-min(pointU(:,1));
if round(horizontal_distance)==round(pointG(:,1))
else
    error('Horizontal distance for ...
 the enclosure is not the same as X coordinate of Point G');
end
if round(horizontal_distance)==round(chordlen)
else
    error('Horizontal distance for  ...
the enclosure is not the same as chord length');
end
%% 
% Determination of other enclosure coordinates w.r.t Point G

chordlen=pointG(1,1);
pointC=[20*pointG(1,1) pointG(1,2) 0 meshsize];
pointB=[pointC(1,1) 12.5*pointG(1,1) 0 meshsize];
pointA=[pointG(1,1) (25-12.5)*pointG(1,1) 0 meshsize];
pointD=[pointC(1,1) -(25-12.5)*pointG(1,1) 0 meshsize];
pointE=[pointG(1,1) -(25-12.5)*pointG(1,1) 0 meshsize];
pointF=[-12.5*pointG(1,1) pointG(1,2) 0 meshsize];
pointH=[pointF(1,1) pointA(1,2) 0 meshsize];
pointI=[pointH(1,1) -pointH(1,2) 0 meshsize];
%% 
% Write Enclosure data to file
% 
% Add enclosure points to file

fprintf(fid,'\r/*Enclosure coordinates w.r.t Point G*/\n');
j=i+499;
for i =1:1:1
    fprintf(fid,'/*PointC*/\n');
    fprintf(fid,'Point(%i) =  ...
{%g, %g, %g, %g};\r\n' ...
,j+1,pointC(i,1),pointC(i,2),...
                                   pointC(i,3),pointC(i,4));
    fprintf(fid,'/*PointB*/\n');
    fprintf(fid,'Point(%i) =  ...
{%g, %g, %g, %g};\r\n' ...
,j+2,pointB(i,1),pointB(i,2),...
                                    pointB(i,3),pointB(i,4));
    fprintf(fid,'/*PointA*/\n');
    fprintf(fid,'Point(%i) =  ...
{%g, %g, %g, %g};\r\n' ...
,j+3,pointA(i,1),pointA(i,2),...
                                   pointA(i,3),pointA(i,4));
    fprintf(fid,'/*PointD*/\n');
    fprintf(fid,'Point(%i) =  ...
{%g, %g, %g, %g};\r\n' ...
,j+4,pointD(i,1),pointD(i,2),...
                                   pointD(i,3),pointD(i,4));
    fprintf(fid,'/*PointE*/\n');
    fprintf(fid,'Point(%i) =  ...
{%g, %g, %g, %g};\r\n' ...
,j+5,pointE(i,1),pointE(i,2),...
                                   pointE(i,3),pointE(i,4));
    fprintf(fid,'/*PointF*/\n');
    fprintf(fid,'Point(%i) =  ...
{%g, %g, %g, %g};\r\n' ...
,j+6,pointF(i,1),pointF(i,2),...
                                   pointF(i,3),pointF(i,4));
    fprintf(fid,'/*PointH*/\n');
    fprintf(fid,'Point(%i) =  ...
{%g, %g, %g, %g};\r\n' ...
,j+7,pointH(i,1),pointH(i,2),...
                                   pointH(i,3),pointH(i,4));
    fprintf(fid,'/*PointI*/\n');
    fprintf(fid,'Point(%i) =  ...
{%g, %g, %g, %g};\r\n' ...
,j+8,pointI(i,1),pointI(i,2),...
                                   pointI(i,3),pointI(i,4));
end
%% 
% Join the enclosure points in the file

fprintf(fid,'\r/*Enclosure creation  ...
from enclosure coordinates*/\n');
k=j;
fprintf(fid,'/*CB*/\n');
fprintf(fid,'Line(%i) = {%i, %i};\r\n',k,k+1,k+2);
fprintf(fid,'/*BA*/\n');
fprintf(fid,'Line(%i) = {%i, %i};\r\n',k+1,k+2,k+3);
fprintf(fid,'/*DE*/\n');
fprintf(fid,'Line(%i) = {%i, %i};\r\n',k+2,k+4,k+5);
fprintf(fid,'/*DC*/\n');
fprintf(fid,'Line(%i) = {%i, %i};\r',k+3,k+4,k+1);
fprintf(fid,'/*AH*/\n');
fprintf(fid,'Line(%i) = {%i, %i};\r',k+4,k+3,k+7);
fprintf(fid,'/*HF*/\n');
fprintf(fid,'Line(%i) = {%i, %i};\r',k+5,k+7,k+6);
fprintf(fid,'/*FI*/\n');
fprintf(fid,'Line(%i) = {%i, %i};\r',k+6,k+6,k+8);
fprintf(fid,'/*IE*/\n');
fprintf(fid,'Line(%i) = {%i, %i};\r',k+7,k+8,k+5);
%% 
% Close the written file

fprintf(fid,'\n');
fclose(fid);
%% Call Gmsh
% |C. Geuzaine and J.-F. Remacle.|
%<https://gmsh.info/doc/preprints/gmsh_paper_preprint.pdf 
% |_Gmsh: a three-dimensional finite element 
%mesh generator with built-in pre- 
% and post-processing facilities_|>|. 
%International Journal for Numerical Methods 
% in Engineering 79(11), pp. 1309-1331, 2009.|
% Gmsh automatic meshing input script
% <http://onelab.info/pipermail
%/gmsh/2017/011027.html 
%http://onelab.info/pipermail/gmsh/2017/011027.html>
% 
% |Mesh options:|
% 
% |'-1, -2, -3'|
% 
% |Perform 1D, 2D or 3D mesh generation, then exit|
% 
% |Mesh save options:|
% 
% |'-save'|
% 
% |Save mesh, then exit|
% 
% |'-save_all'|
% 
% |Save all elements, then exit|
% 
% |'-save_parametric'|
% 
% |Save nodes with their parametric coordinates, then exit|
% 
% |'-save_topology'|
% 
% |Save model topology, then exit|

mesh_options=2;
mesh_save_options='save_all';
gmsh_run_script= sprintf  ...
('gmsh %s %d -%s',saveFilename_geo ...
,-mesh_options,mesh_save_options');
% Run Gmsh using Gmsh input script variable  

while ispc||isunix
[status_geo,result_geo] = system(gmsh_run_script);
    if (status_geo==-1)
        disp("Gmsh is unable to load input script");
        disp(result_geo);
        error(result_geo);
        break
    end
    saveFilename_mesh=append('NACA',NACA,'_gmsh','.msh');
    if isfile(saveFilename_mesh)==1
        movefile (saveFilename_mesh,sevenzip_path,'f');
        movefile (saveFilename_geo,sevenzip_path,'f');
       else
        error("Gmsh did not generate Mesh file")
    end
    cd(sevenzip_path);
    break
end