% updated 06/26/2019
%% Inputs
% r_glb_base : coordinate of the base in the global PCS
% a_glb_base : coordinate of the base in the global PCS
% r_glb_trgd : coordinate of the trgd in the global PCS
% a_glb_trgd : coordinate of the trgd in the global PCS

%% Outputs
% r_loc_trgd : coordinate of the trgd in the local  PCS
% a_loc_trgd : coordinate of the trgd in the local  PCS

%% Example
% r_glb_base = 2.8;
% a_glb_base = pi-pi/6;
% r_glb_trgd = 4.4786;
% a_glb_trgd = 3.5348;
% r_loc_trgd = 3.555277766926234;
% a_loc_trgd = 4.21;

function [r_loc_trgd,a_loc_trgd]=fun_polar_glob2loc(r_glb_base,a_glb_base,r_glb_trgd,a_glb_trgd)

[x_glb_trgd,y_glb_trgd] = pol2cart(a_glb_trgd,r_glb_trgd);
[x_glb_base,y_glb_base] = pol2cart(a_glb_base,r_glb_base);

x_loc_trgd = x_glb_trgd - x_glb_base;
y_loc_trgd = y_glb_trgd - y_glb_base;

[a_loc_trgd,r_loc_trgd] = cart2pol(x_loc_trgd,y_loc_trgd);



% a_loc_trgd      = wrapTo2Pi(a_loc_trgd);
ind             = r_loc_trgd<0;
a_loc_trgd(ind) = wrapTo2Pi(a_loc_trgd(ind) + pi);
r_loc_trgd      = abs(r_loc_trgd);

end

% RLI(:,:,1) = [-(cos(AGI)./2 + 1./2).*(cos(AGL)./2 + 1./2).*((RGI.^2 - 2.*cos(AGI - AGL).*RGI.*RGL + RGL.^2)./(cos(AGI./2).^4.*cos(AGL./2).^4)).^(1./2)]; 
% RLI(:,:,2) = [+(cos(AGI)./2 + 1./2).*(cos(AGL)./2 + 1./2).*((RGI.^2 - 2.*cos(AGI - AGL).*RGI.*RGL + RGL.^2)./(cos(AGI./2).^4.*cos(AGL./2).^4)).^(1./2)];
% 
% ALI(:,:,1) = [-2.*atan((RGL - RGI + RGI.*(cos(AGI) + 1) - RGL.*(cos(AGL) + 1) + (cos(AGI)./2 + 1./2).*(cos(AGL)./2 + 1./2).*((RGI.^2 - 2.*cos(AGI - AGL).*RGI.*RGL + RGL.^2)./(cos(AGI./2).^4.*cos(AGL./2).^4)).^(1./2))./(RGI.*sin(AGI) - RGL.*sin(AGL)))];
% ALI(:,:,2) = [+2.*atan((RGI - RGL - RGI.*(cos(AGI) + 1) + RGL.*(cos(AGL) + 1) + (cos(AGI)./2 + 1./2).*(cos(AGL)./2 + 1./2).*((RGI.^2 - 2.*cos(AGI - AGL).*RGI.*RGL + RGL.^2)./(cos(AGI./2).^4.*cos(AGL./2).^4)).^(1./2))./(RGI.*sin(AGI) - RGL.*sin(AGL)))];
% 
% ALI      = wrapTo2Pi(ALI);
% ind      = RLI<0;
% ALI(ind) = wrapTo2Pi(ALI(ind) + pi);
% RLI      = abs(RLI);
% 
% RLI      = RLI(:,:,1);
% ALI      = ALI(:,:,1);
