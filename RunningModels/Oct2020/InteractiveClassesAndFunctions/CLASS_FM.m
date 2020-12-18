% Forward model functions
classdef CLASS_FM
    methods(Static)
    function MyForwardModel = BeginForwardModel(hb,dE,rr,Ax,dZ,s2nr)
        %------------------------------------------------------------
        f = 300e6; % frequency
        ao = 0:1:179;
        %------------------------------------------------------------
        Zmx = max(hb); % Maximum depth
        Z = [0:dZ:Zmx]'; % depth vector
        Z(1) = 1e-20; % fix the surface depth to not be zero
        Ax = Ax(1:length(hb),:);
        dE = dE(1:length(hb),:);
        rr = rr(1:length(hb),:);
        mp = [dE,rr,Ax]; % model parameters input for the forward model
        P1 = repmat([1 0 ; 0 1],1,1,length(ao));
        P2 = repmat([1 0 ; 0 1],1,1,length(ao));
        %% run the model to obtain 4 signals (HH,VV,HV,VH) in two different mode
        [Dta,P1,P2] = CLASS_FM.RunMainLayers(Z,dZ,hb,ao,mp,f,s2nr,P1,P2);
        %% output
        MyForwardModel.f = f;
        MyForwardModel.ao = ao;
        MyForwardModel.Z = Z;
        MyForwardModel.dZ = dZ;
        MyForwardModel.Zmx = Zmx;
        MyForwardModel.mp = mp;
        MyForwardModel.Dta = Dta;
        MyForwardModel.P1 = P1;
        MyForwardModel.P2 = P2;
    end
%%
    function [Dta,P1,P2] = RunMainLayers(Z,dZ,hb,ao,mp,f,s2nr,P1,P2)
        NL = length(hb); % number of main layers
        iZ0 = 1; % index of the depth at the surface
        for i = 1:NL
            if NL == 1
                iZ1 = length(Z);
                dhb = Z(end)-Z(1);
            else
                [~,iZ1] = min(abs(hb(i)-Z)); % index of current layer boundary
                dhb = mean(diff(hb));
            end
            cZ = Z(iZ0:iZ1); % Depth vector of the current layer 
            %----------------
            [Signal,P1,P2] = CLASS_FM.RunInternalLayers(mp(i,:),f,cZ,dZ,ao,P1,P2);
            SIGNAL(iZ0:iZ1,:,:) = Signal;
            %----------------
            iZ0 = iZ1 + 1;
        end
        if s2nr > 0
            for i = 1:size(SIGNAL,3)
                s = SIGNAL(:,:,i);
                SIGNAL(:,:,i) = CLASS_FM.AddNoise(s,s2nr); %add noise 
            end
%             dnsfactor = s2nr*100*5;
%             NoiseWin=max(hb)*0.05; CohWin=max(hb)*0.1; ConvWin=max(hb)*0.01;
%             NoiseWin=dnsfactor; CohWin=dnsfactor; ConvWin=dnsfactor;
            DenoisingFlag=[]; C_DepthWin=dZ; C_ConvWin=dZ; 
        else
            DenoisingFlag=[]; C_DepthWin=dZ; C_ConvWin=dZ; 
        end
        HH0 = SIGNAL(:,:,1);    VV0 = SIGNAL(:,:,2);   HV0 = SIGNAL(:,:,3);    VH0 = SIGNAL(:,:,4);
        HH = HH0; VV = VV0; HV = HV0; VH = VH0;
%         if size(HH0,2) == length(ao)
%             HH = HH0; VV = VV0; HV = HV0; VH = VH0;
%         elseif size(HH0,2) == 1
%             [HH,VH,HV,VV] = CLASS_S2P.AzimuthSynthesizer(HH0,VH0,HV0,VV0,ao,0);
%         else
%             disp("ERROR: check the 'ao' as internal layers function input")
%             return
%         end
        Dta = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"model");
    end 
%%
    function [Signal,P1out,P2out] = RunInternalLayers(mp,f,cZ,dZ,ao,P1inp,P2inp)
        cZ(cZ==0) = 1e-20;
        da = 0.034; % dielectric anisotropy  
        c = 299792458; % speed of light [m/s]
        epsilon0 = 8.85e-12; % dielectric permitivity in a vacuum [F/m]
        Px = 3.15;
        rx = 1e-12;
        
        Py = (mp(1).*da) + Px;
        rr = 10.^(mp(2)/20); % dB to ratio
%         rr = mp(2); % ratio
        ry = rx.*rr;
        Ax = mp(3);

        Cx = 1e-5;
        Cy = 1e-5;
        omega = 2 * pi * f; % angular frequency [rad/s]
        lambda0 = c / f; % wavelength in a vacuum [m]
        k0 = 2 * pi / lambda0; % [rad/m]
        mu0 = pi*4e-7; % magnetic permeability in a vacuum [Henry/m]
        j = sqrt(-1); % complex number
        
        % variables calculated from model parameters
        dpX1 = Px;
        dpX2 = Cx./(epsilon0.*omega);
        rdpX = dpX1 - (j.*dpX2);
        dpY1 = Py;
        dpY2 = Cy./(epsilon0.*omega);
        rdpY = dpY1 - (j.*dpY2);
        sL = length(cZ); % number of sublayer in the current main layer

        % ao (theta') is the angle between measurement frame and T antenna
        % Ax (alpha) is the angle between measrement frame and ice principal axes
        ti = ao - Ax; % ti (theta_i) is the  angle between ice principal axis and T antenna

        kx = sqrt((epsilon0 .* mu0 .* rdpX .* (omega.^2)) + (j .* mu0 .* Cx .* omega));
        ky = sqrt((epsilon0 .* mu0 .* rdpY .* (omega.^2)) + (j .* mu0 .* Cy .* omega));
        Tx = exp( (-j .* k0 .* dZ) + (j .* kx .* dZ) );
        Ty = exp( (-j .* k0 .* dZ) + (j .* ky .* dZ) );
        T = [Tx complex(0) ; complex(0) Ty];
        G = [rx 0 ; 0 ry];
        HH = nan(sL,size(ti,2));
        VV = nan(sL,size(ti,2));
        HV = nan(sL,size(ti,2));
        VH = nan(sL,size(ti,2));
        P1out = nan(2,2,size(ti,2));
        P2out = nan(2,2,size(ti,2));
        for m = 1:size(ti,2)
            if size(P1inp,3) == 1
                P1 = P1inp;
                P2 = P2inp;
            else
                P1 = P1inp(:,:,m);
                P2 = P2inp(:,:,m);
            end
            dg = ti(1,m);

            n=(1:sL)';
            TT = [T(1,1) T(2,2)].^n;
            D = (exp(j.*k0.*cZ)./(4.*pi.*cZ)).^2;
            csd = cosd(dg);
            snd = sind(dg);
            PRGRP = P1 * CLASS_FM.R(dg)*G*CLASS_FM.R(-dg) * P2;
            A1 = PRGRP(1,1);
            A2 = PRGRP(1,2);
            A3 = PRGRP(2,1);
            A4 = PRGRP(2,2);
            a1 = ((csd.^2) .* TT(:,1) ) + ((snd.^2) .* TT(:,2) );
            a2 = csd .* snd .* (TT(:,1) - TT(:,2));
            a4 = ((snd.^2) .* TT(:,1) ) + ((csd.^2) .* TT(:,2) );
            
            HH(:,m) = D.* ( a1.*(a1.*A1 + a2.*A3) + a2.*(a1.*A2 + a2.*A4));
            VH(:,m) = D.* ( a2.*(a1.*A1 + a2.*A3) + a4.*(a1.*A2 + a2.*A4));
            HV(:,m) = D.* ( a1.*(a2.*A1 + a4.*A3) + a2.*(a2.*A2 + a4.*A4));
            VV(:,m) = D.* ( a2.*(a2.*A1 + a4.*A3) + a4.*(a2.*A2 + a4.*A4));

            P1out(:,:,m) =  ((CLASS_FM.R(dg)*T*CLASS_FM.R(-dg))^(sL)) * P1;
            P2out(:,:,m) =  P2 * ((CLASS_FM.R(dg)*T*CLASS_FM.R(-dg))^(sL));
        end
        Signal(:,:,1) = HH;
        Signal(:,:,2) = VV;
        Signal(:,:,3) = HV;
        Signal(:,:,4) = VH;
    end
%%
    function R_mat = R(deg)
             R_mat = [cosd(deg) -sind(deg);
                             sind(deg)  cosd(deg)];
    end  
%%
    function [Nsgnl] = AddNoise(sgnl,pr)
        % Add random synthetic noise to the synthetic data
%         j = sqrt(-1);
%         re = real(sgnl);
%         im = imag(sgnl);
%         r_sigmas = pr .* re;
%         i_sigmas = pr .* im;
%         re_randomNoise = randn(size(sgnl)) .* r_sigmas;
%         im_randomNoise = randn(size(sgnl)) .* i_sigmas;
%         Nsgnl = (re + re_randomNoise) + j.*(im+im_randomNoise);
        absS = abs(sgnl);
        angS = angle(sgnl);
        abs_sigmas = pr .* absS;
        ang_sigmas = pr .* angS;
        abs_randomNoise = randn(size(sgnl)) .* abs_sigmas;
        ang_randomNoise = randn(size(sgnl)) .* ang_sigmas;
        Nsgnl = (absS + abs_randomNoise) .* exp(i.*(angS+ang_randomNoise));
    end
%%
function [rnd_ab] = RandBetween2Values(a,b,n)
    rnd_ab = (b-a).*rand(n) + a;
end
%%
end
end