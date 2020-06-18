clc
clear
close all

%% parameters
d_star          = 3*.249e-9; 
cell_res        = 5e-6;
cell_max        = (sqrt(2)/2)*cell_res*2;
cell_maxd       = cell_max/2;
data_num        = 1000000;
flb             = -1.3387e+09;
fub             = +1.3387e+09;

%% functions 
fun_scale       = @(x,lb,ub) ((x-(ones(size(x,1),1)*min(x))) ...
                    ./((ones(size(x,1),1)*max(x))-(ones(size(x,1),1) ...
                    *min(x))))*(ub-lb)+lb; 
                
% a, r, and g are for target in pcs2 and t is for the base
fun_fg  = @(a,r,g,t) (cos(2*g).*(cos(3*a - 2*t) + ...
                      cos(a - 2*t)))./(2*r) ...
                      + (4*cos(2*t).*cos(a).^2.*cos(g).*sin(a).*sin(g))./r;

fun_fg_inv =@(a,r,g,t) ((cos(2*t).*(cos(3*a - 2*g) + cos(a - 2*g)))./2 + ...
                         4*cos(2*g).*cos(a).^2.*sin(a).*cos(t).*sin(t))./r;

%% restart random generator
rng('shuffle')    ;
s = rng;

%% data
data.num = data_num;

% base in pcs2
data.pcs2.bse.a = rand(data.num,1)*2*pi;
data.pcs2.bse.t = rand(data.num,1)*2*pi;
data.pcs2.bse.r = zeros(data.num,1);

% target in pcs2
data.pcs2.trg.a = rand(data.num,1)*2*pi;
data.pcs2.trg.t = rand(data.num,1)*2*pi;
% - generate forces using a URG
data.pcs2.trg.f = fun_scale(rand(data.num,1),flb,fub);
% - use inverse function to generate radials 
data.pcs2.trg.r = fun_fg_inv(data.pcs2.trg.a, ...
                             data.pcs2.trg.f, ... 
                             data.pcs2.trg.t, ...
                             data.pcs2.bse.t);
                         
%% transform target data with negetive r
ind_neg                    = data.pcs2.trg.r<0;
% absolute r
data.pcs2.trg.r(ind_neg,1) = abs(data.pcs2.trg.r(ind_neg,1));
% azimuth - pi
data.pcs2.trg.a(ind_neg,1) = data.pcs2.trg.a(ind_neg,1) - pi;
% wrap between 0 and 2pi
data.pcs2.trg.a            = wrapTo2Pi(data.pcs2.trg.a);

%% remove data that violates r>=d_star
ind_bad                    = data.pcs2.trg.r<d_star;

data.pcs2.bse.r(ind_bad,:) = [];
data.pcs2.bse.a(ind_bad,:) = [];
data.pcs2.bse.t(ind_bad,:) = [];

data.pcs2.trg.r(ind_bad,:) = [];
data.pcs2.trg.a(ind_bad,:) = [];
data.pcs2.trg.t(ind_bad,:) = [];
data.pcs2.trg.f(ind_bad,:) = [];

% update data.num
data.num = length(data.pcs2.bse.r);

%% transform data into pcs3
% base in pcs3
data.pcs3.bse.a = rand(data.num,1)*2*pi;
data.pcs3.bse.t = data.pcs2.bse.t; % the same in pcs2 and 3
data.pcs3.bse.r = fun_scale(rand(data.num,1),0,cell_maxd);


% target in pcs3
[data.pcs3.trg.r, data.pcs3.trg.a] = fun_polar_loc2glob(data.pcs3.bse.r, ...
                                                        data.pcs3.bse.a, ...
                                                        data.pcs2.trg.r, ...
                                                        data.pcs2.trg.a);


data.pcs3.trg.t = data.pcs2.trg.t;
data.pcs3.trg.f = data.pcs2.trg.f;

%% remove data outside of cell boundries 
ind_bad = data.pcs3.bse.r > cell_maxd | data.pcs3.trg.r > cell_maxd;

data.pcs2.bse.r(ind_bad,:) = [];
data.pcs2.bse.a(ind_bad,:) = [];
data.pcs2.bse.t(ind_bad,:) = [];

data.pcs2.trg.r(ind_bad,:) = [];
data.pcs2.trg.a(ind_bad,:) = [];
data.pcs2.trg.t(ind_bad,:) = [];
data.pcs2.trg.f(ind_bad,:) = [];

data.pcs3.bse.r(ind_bad,:) = [];
data.pcs3.bse.a(ind_bad,:) = [];
data.pcs3.bse.t(ind_bad,:) = [];

data.pcs3.trg.r(ind_bad,:) = [];
data.pcs3.trg.a(ind_bad,:) = [];
data.pcs3.trg.t(ind_bad,:) = [];
data.pcs3.trg.f(ind_bad,:) = [];

%% plot
figure;
hist(data.pcs3.bse.t,100)
figure;
hist(data.pcs3.trg.t,100)
figure;
hist(data.pcs3.trg.a,100)
figure;
hist(data.pcs3.trg.f,100)
figure;
polarplot(data.pcs3.bse.a, data.pcs3.bse.r,'b.','markersize',1)
hold on
polarplot(data.pcs3.trg.a, data.pcs3.trg.r,'r.','markersize',1)
rlim([0,cell_maxd])

%% inputs & outputs 
v = [data.pcs3.bse.t       , ...
     data.pcs2.trg.a       , ...
     log10(data.pcs2.trg.r), ...
     data.pcs2.trg.t];

 u = data.pcs3.trg.f;
 
% scaling
datain_raw = v;
dataou_raw = u;

datain_scl = [data.pcs3.bse.t/(2*pi), ...
              data.pcs2.trg.a/(2*pi), ...
              fun_scale(log10(data.pcs2.trg.r),0,1), ...
              data.pcs2.trg.t/(2*pi)];

dataou_scl = u/max(abs(u));
          
%% save 
% save in .mat
save('data_raw','datain_raw','dataou_raw','-v7.3')
save('data_scl','datain_scl','dataou_scl','-v7.3')

% save in dict in .txt for python
% raw
datain_list = fun_mat2list(datain_raw);
dataou_list = fun_mat2list(dataou_raw);
dict_val    = append("{'datain': ",string(datain_list),", 'dataou': ", string(dataou_list), "}");
fid         = fopen('data_inou_raw.txt','wt');
fprintf(fid, dict_val);
fclose(fid);

% scl
datain_list = fun_mat2list(datain_scl);
dataou_list = fun_mat2list(dataou_scl);
dict_val    = append("{'datain': ",string(datain_list),", 'dataou': ", string(dataou_list), "}");
fid         = fopen('data_inou_scl.txt','wt');
fprintf(fid, dict_val);
fclose(fid);
