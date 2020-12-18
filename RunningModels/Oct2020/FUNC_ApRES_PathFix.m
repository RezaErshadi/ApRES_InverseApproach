function [DataPath,RadarData,InvPath,ps] = FUNC_ApRES_PathFix()
if ismac || isunix % path seperature
    ps = '/';
elseif ispc
    ps = '\';
end
CurDir = pwd; % current directory
SplitPath = split(CurDir,ps);
NumFolder = length(SplitPath);
for i = NumFolder:-1:1
    if SplitPath{i} == "ApRES"
        break
    else
        cd ..
    end
end
HomeDir = pwd;
% Add fiex classes/functions to the path
CalssesAndFunctionsPath = strcat(HomeDir,ps,'FixedClassesAndFunctions');
addpath(genpath(CalssesAndFunctionsPath));
% Add interactive classes/functions to the path
InteractiveClassesFunctions = strcat(CurDir,ps,'InteractiveClassesAndFunctions'); 
addpath(genpath(InteractiveClassesFunctions)); 
% Add Radar data to the apth
DataPath = strcat(HomeDir,ps,'Data'); 
addpath(genpath(DataPath)); 

RadarFolder = dir(strcat(DataPath,ps,'RadarData_ApRES'));
for i = 1:length(RadarFolder)
    RadarData(i,1) = string(RadarFolder(i).name);
    RadarData(i,2) = string(RadarFolder(i).folder);
end
RadarData(RadarData(1,:)==".",:) = [];
RadarData(RadarData(1,:)=="..",:) = [];

InvPath = strcat(HomeDir,ps,'InversionResults');

cd(CurDir)