% updated 06/26/2019
%% Inputs
% r_glb_base : coordinate of the base in the global PCS
% a_glb_base : coordinate of the base in the global PCS
% r_loc_trgd : coordinate of the trgd in the local  PCS
% a_loc_trgd : coordinate of the trgd in the local  PCS

%% Outputs
% r_glb_trgd : coordinate of the trgd in the global PCS
% a_glb_trgd : coordinate of the trgd in the global PCS

% RGI : point of interest in global PCS
% AGI : point of interest in global PCS

%% Example
% r_glb_base = 2.8;
% a_glb_base = pi-pi/6;
% r_glb_trgd = 4.4786;
% a_glb_trgd = 3.5348;
% r_loc_trgd = 3.555277766926234;
% a_loc_trgd = 4.21;

function [r_glb_trgd,a_glb_trgd]=fun_polar_loc2glob(r_glb_base,a_glb_base,r_loc_trgd,a_loc_trgd)

[x_loc_trgd,y_loc_trgd] = pol2cart(a_loc_trgd,r_loc_trgd);
[x_glb_base,y_glb_base] = pol2cart(a_glb_base,r_glb_base);

x_glb_trgd = x_loc_trgd + x_glb_base;
y_glb_trgd = y_loc_trgd + y_glb_base; 

[a_glb_trgd,r_glb_trgd] = cart2pol(x_glb_trgd,y_glb_trgd);

% a_glb_trgd      = wrapTo2Pi(a_glb_trgd);
ind             = r_glb_trgd<0;
a_glb_trgd(ind) = wrapTo2Pi(a_glb_trgd(ind) + pi);
r_glb_trgd      = abs(r_glb_trgd);

end

% RGI(:,:,1) = +(cos(AGL)./2 + 1./2).*(cos(ALI)./2 + 1./2).*((RGL.^2 + 2.*cos(AGL - ALI).*RGL.*RLI + RLI.^2)./(cos(AGL./2).^4.*cos(ALI./2).^4)).^(1./2);
% RGI(:,:,2) = -(cos(AGL)./2 + 1./2).*(cos(ALI)./2 + 1./2).*((RGL.^2 + 2.*cos(AGL - ALI).*RGL.*RLI + RLI.^2)./(cos(AGL./2).^4.*cos(ALI./2).^4)).^(1./2);
% 
% AGI(:,:,1) = +2.*atan((RGL + RLI - RGL.*(cos(AGL) + 1) - RLI.*(cos(ALI) + 1) + (cos(AGL)./2 + 1./2).*(cos(ALI)./2 + 1./2).*((RGL.^2 + 2.*cos(AGL - ALI).*RGL.*RLI + RLI.^2)./(cos(AGL./2).^4.*cos(ALI./2).^4)).^(1./2))./(RGL.*sin(AGL) + RLI.*sin(ALI)));
% AGI(:,:,2) = -2.*atan((RGL.*(cos(AGL) + 1) - RLI - RGL + RLI.*(cos(ALI) + 1) + (cos(AGL)./2 + 1./2).*(cos(ALI)./2 + 1./2).*((RGL.^2 + 2.*cos(AGL - ALI).*RGL.*RLI + RLI.^2)./(cos(AGL./2).^4.*cos(ALI./2).^4)).^(1./2))./(RGL.*sin(AGL) + RLI.*sin(ALI)));
% 
% AGI      = wrapTo2Pi(AGI);
% ind      = RGI<0;
% AGI(ind) = wrapTo2Pi(AGI(ind) + pi);
% RGI      = abs(RGI);
% 
% RGI      = RGI(:,:,1);
% AGI      = AGI(:,:,1);
