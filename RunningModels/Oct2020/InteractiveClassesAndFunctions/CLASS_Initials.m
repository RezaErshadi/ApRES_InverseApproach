classdef CLASS_Initials
    methods(Static)
%  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
function [AxOut] = PossiblePrincipalAxis(HV,ao)
%     HV = CLASS_S2P.PCA_denoising(HV,1);
    [~,iB] = sort(HV,2); % sort HV and the first and second index of the smallest HV at each depth
    imin = iB(:,1:2); % The indices of the minimum of HV
    for i = 1:size(iB,1)
        A(i,1) = ao(iB(i,1));
        A(i,2) = ao(iB(i,2));
    end
    A = sort(A,2); % Angle of minimum HV
%  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
% two separate min angle
    psbl_Ax = [A-180 A A+180]; 
    a = psbl_Ax(1,4);
    for i = 1:size(HV,1)-1
        b = psbl_Ax(i+1,:);
        [~,ia] = min(abs(a(i,1)-b));
        a(i+1,1) = b(ia);
    end
%  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
    AxMin = a;
    AxMean = mean(a);
    AxStd = std(a);
%  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
    AxMean(2) = AxMean(1)-90;
    AxMin(:,2) = AxMin(:,1)-90;
%  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
    AxMin(AxMin>179) = AxMin(AxMin>179)-180;
    AxMin(AxMin<0) = AxMin(AxMin<0)+180;
    AxMean(AxMean>179) = AxMean(AxMean>179)-180;
    AxMean(AxMean<0) = AxMean(AxMean<0)+180;  
%  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
    if AxMean(2) < AxMean(1)
        AxMean(1,[1 2]) = AxMean(1,[2 1]);
        AxMin(:,[1 2]) = AxMin(:,[2 1]);
    end
    imin = round(AxMin+1,0);
%  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
    AxOut.AxMin = round(AxMin,0);
    AxOut.AxMean = round(AxMean,0);
    AxOut.AxStd = round(AxStd,0);
    AxOut.imin = imin;
end
%%  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
function dE0 = PossibledE(Z,Zmdl,AxOut,dEhor)
    psbl_dE = nan(size(Z));
    for i = 1:length(Z)
        psbl_dE(i,1) = dEhor(i,AxOut.imin(i,1));
    end
    dE1 = abs(psbl_dE);
    dE1(dE1>0.5) = nan;
    dE1 = fillmissing(dE1,'nearest');
    dE1 = smoothdata(dE1,'movmean');
    j1 = 1;
    for i = 1:size(Zmdl,1)
        [~,j2] = min(abs(Zmdl(i,1)-Z));
        dE0(i,1) = nanmean(dE1(j1:j2,1));
        j1 = j2+1;
    end
end
%%  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    end
end