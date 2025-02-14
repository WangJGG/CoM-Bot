clear
%% select the folder with data in it
data_directory ='..\static\uploads';
cd(data_directory);
avifiles = dir('*.avi');
cd('..\..\ViconMarkerless2Opensim\uploads');
% get names of all of the any existing csvfiles
mp_path = fileparts(which('mediapipe_opensim.m'));
if count(py.sys.path,mp_path) == 0
    insert(py.sys.path,int32(0),mp_path)
end
run_mediapipe='Y';
s=pwd;
cd('..');
if strcmp(run_mediapipe,'Y')
    for F = 1:length(avifiles)
        disp(['Processing ' avifiles(F).name ' .....']);
        py.mediapipe_pose_estimation_markerless.mediapipe_video_analysis([avifiles(F).folder '\' avifiles(F).name]);
    end
end
cd(s);
% get names of all of the newly made csvfiles
csvfiles = dir('*.csv');
%% loop through csv files and open them and extract the landmark information and store to data structure
for F = 1:length(csvfiles)
    % load data fle
    data.data2d{F} = load_mediapipe_csv(csvfiles(F).name);
    % get frame_np time stamps for filtering and gap filling
    frame_no = table2array(data.data2d{F}(:,1));    
    time = table2array(data.data2d{F}(:,2));
    dt = median(diff(time));
    time = (1:frame_no(end))'*dt;    
    
    idx=ismember(1:numel(time),frame_no);
    
    for i = 3:4:size(data.data2d{1},2)
        MNAME = data.data2d{F}.Properties.VariableNames{i};
        MNAME = MNAME(1:end-2);        
        % get data from the table and interpolate any missing points from
        % the first to last frame
        XY = table2array(data.data2d{F}(:,i:i+1));
        % undistort points
        %XY = undistortPoints(XY,cam.cam_dist{F});
        % interpolate across any missing frames
        X = (interp1(frame_no,XY(:,1),1:frame_no(end),'linear'))';
        Y = (interp1(frame_no,XY(:,2),1:frame_no(end),'linear'))';
        C = (interp1(frame_no,table2array(data.data2d{F}(:,i+3)),1:frame_no(end),'linear'))';        
        X(~idx) = 0;
        Y(~idx) = 0;
        C(~idx) = 0;
        % low pass filter
        X_filt = matfiltfilt_low(dt,6,2,X); % 6Hz low pass filter
        Y_filt = matfiltfilt_low(dt,6,2,Y); % 6Hz low pass filter
        %store to marker structure
        data.markers2D.(MNAME).XY(:,1,F) = X_filt;
        data.markers2D.(MNAME).XY(:,2,F) = Y_filt;
        data.markers2D.(MNAME).C(:,F) = C;
    end
end

% savve the frame no, time, frame_rate to data structure
data.frame(:,1) = frame_no(1):1:frame_no(end);
data.time = time;
data.frame_rate = 1/dt;
% get information about the number of frames and markers for doing
% triangulation
nframes = size(data.data2d{1},1);
markers = fieldnames(data.markers2D);

%generate the videos to get 3d coordinates
s=pwd;
videopath = "..\static\uploads\front.avi";
%disp(videopath);
cd(s);
resultpath='result.csv';
run_pose_3d = 'Y';
cd('..');
if strcmp(run_pose_3d,'Y')
    py.mediapipe_pose_point.extract_pose_from_video(videopath,resultpath);
end
cd(s);
fid = fopen(resultpath, 'rt');
numMarkers = 33;
fileContent = textscan(fid, '%s', 'Delimiter', '\n');
allLines = fileContent{1};
numFrames = sum(cellfun('isempty', allLines));
dataMatrix = readmatrix(resultpath);
data.XYZ = zeros(numMarkers, 3, numFrames);
dataRowIndex = 1;

for frameIndex = 1:numFrames
    if dataRowIndex + numMarkers - 1 > size(dataMatrix, 1)
        break;
    end
    frameData = dataMatrix(dataRowIndex:dataRowIndex + numMarkers - 1, :);
    dataRowIndex = dataRowIndex + numMarkers;
    data.XYZ(:, :, frameIndex) = frameData;
end

fclose(fid);

Fname = ['squat.trc'];
mediapipepose2trc(data, Fname);

%% process opensim model - start by scaling generic model
import org.opensim.modeling.*

mass=48;
p_name='username';

generic_model = '..\Rajagopal2015_mediapipe_head.osim';
scale_settings_file = '..\scale_settings_mediapipe.xml';
scaled_model = [cd '\ModelScaling\' p_name '_modelScaled' '.osim'];

disp('Scale Tool Processing....')

%load scale tool and associated tools
ScTool = ScaleTool(scale_settings_file);
ScTool.setName(p_name)
ScTool.setPathToSubject('');
GMM = ScTool.getGenericModelMaker;
MS = ScTool.getModelScaler;
MP = ScTool.getMarkerPlacer;

%add model to generic setup file
GMM.setModelFileName(generic_model);
GMM.setName(generic_model);

%add output files
MS.setOutputScaleFileName([cd '\ModelScaling\' p_name '_scaleSetApplied' '.xml']);
MS.setOutputModelFileName([cd '\ModelScaling\' p_name '_modelScaled' '.osim']);
MP.setOutputMotionFileName([cd '\ModelScaling\' p_name '_static' '.mot']);
MP.setOutputModelFileName([cd '\ModelScaling\' p_name '_modelScaled' '.osim']);
MP.setOutputMarkerFileName([cd '\ModelScaling\' p_name '_markersScaled' '.xml']);
ScaledModel = [cd '\ModelScaling\' p_name '_modelScaled' '.osim'];

% setup the specific file parameters
% note that the name and path of the file is in the data_structure for static file
MS.setMarkerFileName(Fname);
MP.setMarkerFileName(Fname);

markerData = MarkerData(Fname);
% Get initial and intial time
InitialTime = markerData.getStartFrameTime();
FinalTime = markerData.getLastFrameTime();

%add time range
time_range = ArrayDouble;
time_range.append(FinalTime-0.1);
time_range.append(FinalTime);
MS.setTimeRange(time_range);
MP.setTimeRange(time_range);

%add mass
ScTool.setSubjectMass(mass);

%create the ModelScaling directory above c3d directory if it doesn't exist
if isempty(dir([cd '\ModelScaling']))
    mkdir(cd,'ModelScaling');
end
%write new .xml file in setup folder
ScTool.print([cd '\ModelScaling\' p_name '_setupScale.xml']);
% run scaling tool;
ScTool.run();
disp('DONE.');
%%
disp('IK Tool Processing....');
data_path = cd;

ik_settings_file = '..\IK_settings_mediapipe.xml';

ikTool = InverseKinematicsTool(ik_settings_file);

model = Model(scaled_model);

% Tell Tool to use the loaded model
ikTool.setModel(model);

%create the InverseKinematics directory above c3d directory if it
%doesn't exist
if isempty(dir([data_path '\InverseKinematics']))
    mkdir(data_path,'InverseKinematics');
end

% define the file names
marker_file = Fname;
[~,fname,~] = fileparts(Fname);
mot_file = [data_path '\InverseKinematics\' fname '.mot'];

setup_file = [data_path '\InverseKinematics\' fname '_iksetup.xml'];
% Get trc data to determine time range
markerData = MarkerData(Fname);
% Get initial and intial time
initial_time = markerData.getStartFrameTime();
final_time = markerData.getLastFrameTime();
% Setup the ikTool for this trial
ikTool.setMarkerDataFileName(marker_file);
ikTool.setStartTime(initial_time);
ikTool.setEndTime(final_time);
ikTool.setOutputMotionFileName(mot_file);
ikTool.setName(fname)
ikTool.setResultsDir([data_path '\InverseKinematics\']);

%write the XML setup file in same directory as MOT file
ikTool.print(setup_file);
disp(['Processing ' marker_file '....']);
% Run IK via API
ikTool.run();
disp(['DONE.']);

%% run body kinematics

body_settings_file = '..\body_kinematics_settings_mediapipe.xml';
           
%create the InverseKinematics directory above c3d directory if it
%doesn't exist
if isempty(dir([data_path '\BodyKinematics']))
    mkdir(data_path,'BodyKinematics');
end

hd = ['time' '\t' 'vertical_acc' '\t' 'horizontal_acc' '\t' 'vertical_force' '\t' 'horizontal_force' '\n'];
fm = ['%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\n'];

[~,fname,~] = fileparts(Fname);
mot_file = [data_path '\InverseKinematics\' fname '.mot'];
setup_file = [data_path '\BodyKinematics\' fname '_BKsetup.xml'];

tool = AnalyzeTool(body_settings_file);
tool.setModelFilename(scaled_model);
tool.setCoordinatesFileName(mot_file);
tool.setSolveForEquilibrium(false);
tool.setName(fname);
tool.setResultsDir([data_path '\BodyKinematics\' fname '\']);
motData = Storage(mot_file);
initial_time = motData.getFirstTime();
final_time = motData.getLastTime();
tool.setStartTime(initial_time);
tool.setFinalTime(final_time);
AS = tool.getAnalysisSet();
BK = AS.get(0);
BK.setStartTime(initial_time);
BK.setEndTime(final_time);

tool.print(setup_file);

disp(['Processing ' mot_file '....']);
tool.run();
disp('DONE.');

data1 = load_sto_file([data_path '\BodyKinematics\' fname '\' fname '_BodyKinematics_pos_global.sto']);
figure('Visible','off');
plot(data1.time, data1.center_of_mass_X,'b-', ...
    data1.time,data1.center_of_mass_Y,'r-', ...
    data1.time,data1.center_of_mass_Z,'g-');
legend('X','Y','Z');
xlabel('Time(s)'); 
ylabel('pos coondinates(m)');
title(['Pos of ' p_name ' Mass Center']);
saveas(gcf,['Mass_Center_of_' p_name '.png']);

disp('Successful');
output = 'Successful';
%%

filename1 = '.\BodyKinematics\squat\squat_BodyKinematics_pos_global.sto';
filename2 = '..\static_result\BodyKinematics\squat\squat_BodyKinematics_pos_global.sto';
opts1 = detectImportOptions(filename1, 'FileType', 'text', 'NumHeaderLines', 18);
opts2 = detectImportOptions(filename2, 'FileType', 'text', 'NumHeaderLines', 18);
opts1.SelectedVariableNames = {'time', 'center_of_mass_X', 'center_of_mass_Y', 'center_of_mass_Z'};
opts2.SelectedVariableNames = {'time', 'center_of_mass_X', 'center_of_mass_Y', 'center_of_mass_Z'};

com1 = readtable(filename1, opts1);
com2 = readtable(filename2, opts2);

x1 = com1.center_of_mass_X;
y1 = com1.center_of_mass_Y;
z1 = com1.center_of_mass_Z;

x2 = com2.center_of_mass_X;
y2 = com2.center_of_mass_Y;
z2 = com2.center_of_mass_Z;

folderName = 'images';

figure('Visible','off');
plot3(com1.center_of_mass_X, com1.center_of_mass_Y, com1.center_of_mass_Z, 'r', 'DisplayName', 'COM Group 1');
hold on;
plot3(com2.center_of_mass_X, com2.center_of_mass_Y, com2.center_of_mass_Z, 'b', 'DisplayName', 'COM Group 2');
xlabel('X Position');
ylabel('Y Position');
zlabel('Z Position');
title('3D Trajectory Comparison of COM Movement');
legend show;
grid on;
fileName = '3D_Trajectory.png';
fullFileName = fullfile(pwd, folderName, fileName);
saveas(gcf, fullFileName);

figure('Visible','off');
plot(com1.time, x1,'b-',com1.time,y1,'r-', com1.time,z1,'g-');
legend('X','Y','Z');
xlabel('Time(s)'); 
ylabel('pos coondinates(m)');
title(['Pos of ' p_name ' Mass Center']);
fileName = 'Mass_of_Center.png';
fullFileName = fullfile(pwd, folderName, fileName);
saveas(gcf, fullFileName);

v = VideoReader('..\..\static\uploads\front.avi');

s=pwd;
cd('..');
vw = VideoWriter('..\Result\result', 'MPEG-4');
open(vw);

cd(s);
filename4 = '.\front.csv';
opts4=detectImportOptions(filename4);
opts4.SelectedVariableNames = {'LEFT_HIP_X','LEFT_HIP_Y','RIGHT_HIP_X','RIGHT_HIP_Y'};
center1 = readtable(filename4, opts4);
filename5 = '..\static_result\front.csv';
opts5=detectImportOptions(filename5);
opts5.SelectedVariableNames = {'LEFT_HIP_X','LEFT_HIP_Y','RIGHT_HIP_X','RIGHT_HIP_Y'};
center2 = readtable(filename5, opts5);

frameIndex = 0;
while hasFrame(v)
    frame = readFrame(v);
    frameIndex = frameIndex + 1;

    if frameIndex <= height(com1)&&frameIndex <= height(com2)&&frameIndex <= height(center1)&&frameIndex <= height(center2)
        offsetx=(center1.LEFT_HIP_X(1)+center1.RIGHT_HIP_X(1))/2-(center2.LEFT_HIP_X(1)+center2.RIGHT_HIP_X(1))/2;
        offsety=(center1.LEFT_HIP_Y(1)+center1.RIGHT_HIP_Y(1))/2-(center2.LEFT_HIP_Y(1)+center2.RIGHT_HIP_Y(1))/2;
        x1 = com1.center_of_mass_X(frameIndex)*500+(center1.LEFT_HIP_X(frameIndex)+center1.RIGHT_HIP_X(frameIndex))/2;
        y1 = -com1.center_of_mass_Z(frameIndex)*500+(center1.LEFT_HIP_Y(frameIndex)+center1.RIGHT_HIP_Y(frameIndex))/2;
        x2 = com2.center_of_mass_Y(frameIndex)*10+(center2.LEFT_HIP_X(frameIndex)+center2.RIGHT_HIP_X(frameIndex))/2+offsetx;
        y2 = com2.center_of_mass_Z(frameIndex)*10+(center2.LEFT_HIP_Y(frameIndex)+center2.RIGHT_HIP_Y(frameIndex))/2+offsety;

        frame = insertShape(frame, 'FilledCircle', [x1, y1, 20], 'Color', 'green', 'LineWidth', 5);
        frame = insertShape(frame, 'Circle', [x2, y2, 40], 'Color', 'blue', 'LineWidth', 5);
    end
    
    writeVideo(vw, frame);
end

close(vw);