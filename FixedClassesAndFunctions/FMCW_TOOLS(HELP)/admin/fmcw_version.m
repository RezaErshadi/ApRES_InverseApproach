function ver = fmcw_version(verpath)
 
% ver = fmcw_version
%
% Display release date of current software version
%
% Craig Stewart

if nargin ==0
    verpath = []; % defaul to current path
end
[fid,mes] = fopen([verpath 'version.txt']);
if fid == -1
    error(mes)
else
    ver = fgetl(fid);
    fclose(fid);
    disp(['FMCW code release ' ver])
end