% SINDy code with figure
close all; clear vars; clc;
% add path and load data
addpath('./utils');
load('NF.csv') % data cells of each condition 
load('NM.csv')
load('OF.csv')
load('OM.csv')
%% prep data:

N = {NF(:,1), NM(:,1), OF(:,1), OM(:,1)};  % number of bats
C = {NF(:,4), NM(:,4), OF(:,4), OM(:,4)};  % total calls
P = {NF(:,5), NM(:,5), OF(:,5), OM(:,5)};  % total power

% parameters for sindy
polyorder = 2;  % search space up to second order polynomials
usesine = 0;    % no trig functions
n = 3;          % 3D system 

% initialize
x = cell(1,4);
dx = cell(1,4);
Theta = cell(1,4);
norm_dx = cell(1,4);
norm_Theta = cell(1,4);

for i=1:4
    x{:,i} = [N{i} C{i} P{i}];
    dx{:,i} = [diff(x{:,i})];
    x{:,i} = x{:,i}(1:end-1,:);
    Theta{:,i} = poolData(x{:,i},n,polyorder,usesine);

    % normalize using std and mean
    norm_dx{:,i} = dx{:,i}; %initialize
    std_dx = std(dx{:,i});
    mean_dx = mean(dx{:,i});
    norm_Theta{:,i} = Theta{:,i}; %initialize
    std_the = std(Theta{:,i});
    mean_the= mean(Theta{:,i});
    for c = 1:3 % based on 3D system
        norm_dx{1,i}(:,c) = (norm_dx{1,i}(:,c)-mean_dx(c)) ./ std_dx(c);
    end
    for c= 2:10 % based on polyorder 
        norm_Theta{1,i}(:,c) = (norm_Theta{1,i}(:,c) - mean_the(c)) ./ std_the(c);
    end
end


%% split data (train, test, validate)
% initialize
Train_dx = cell(1,4);
Test_dx = cell(1,4);
Val_dx = cell(1,4);
Train_Theta = cell(1,4);
Test_Theta = cell(1,4);
Val_Theta = cell(1,4);

% choose sections of data that represent the system well
for i = 1:4
    Val_dx{1,i} = norm_dx{1,i}(1:450,:);      % 12.5%
    Test_dx{1,i} = norm_dx{1,i}(451:900,:);   % 12.5%
    Train_dx{1,i} = norm_dx{1,i}(901:end,:);  % 75%
    Val_Theta{1,i} = norm_Theta{1,i}(1:450,:);
    Test_Theta{1,i} = norm_Theta{1,i}(451:900,:);
    Train_Theta{1,i} = norm_Theta{1,i}(901:end,:);
end
 
%% test sindy for lambda
% initialize
mean_error = [];
l = [];

% what lambda minimizes MSE:
for i = 1:4
    lambda = [0:.001:1];
    for z = 1: length(lambda)
        XI = sparsifyDynamics(Test_Theta{1,i},Test_dx{1,i},lambda(z),n);
        calc_dx = Test_Theta{1,i}*XI;
        true_dx = Test_dx{1,i};
        mean_error(z) = mean(sqrt(sum((calc_dx(:,2:3) - true_dx(:,2:3)).^2,2))); % only using C and P
        temp = abs(XI) >10^(-10);
        s(z) = sum(sum(temp(:,2:3)));
    end
    [M,I] = min(mean_error);
    l(i) = lambda(I);
    % figure
    % plot(lambda, mean_error);
    % xlabel('lambda'), ylabel('MSE')
end
% we adjusted lambda values based on the error; maximize up to .003 from
% minimum error
%% train sindy using lambda found above

% initialize matrices
Xi = cell(1,4);
n=3; % 3D system
% Train SINDy
for i=1:4 
    Xi{1,i} = sparsifyDynamics(Train_Theta{1,i},Train_dx{1,i},l(i),n);
end

%% Validate SINDy
% initialize
Xdot={};
Rsqu = zeros(2,4);

for k=1:4
    % equations based on Xi (polyorder 2)
    Xdot{1,k} = Val_Theta{1,k}*Xi{1,k};

    % save R squared values
    m=fitlm(Xdot{1,k}(:,2),Val_dx{1,k}(:,2));
    Rsqu(1,k)=m.Rsquared.Ordinary;
    m=fitlm(Xdot{1,k}(:,3),Val_dx{1,k}(:,3));
    Rsqu(2,k)=m.Rsquared.Ordinary;


end

%%
% % colormaps of Xi for cdot and pdot
Cmat = [Xi{1,1}(:,2) Xi{1,2}(:,2) Xi{1,3}(:,2) Xi{1,4}(:,2)];
Pmat = [Xi{1,1}(:,3) Xi{1,2}(:,3) Xi{1,3}(:,3) Xi{1,4}(:,3)];

figure
tiledlayout(1,2,"TileSpacing","compact");
nexttile
imagesc(Cmat')
title('Total Calls')
set(gca,'yTick',1:4,'yticklabel',{'NF','NM','OF','OM'},'xTick',1:10,...
    'xticklabel',{'1','{\it N}','{\it C}','{\it P}','{\it N}^2','{\it NC}', '{\it NP}','{\it C}^2','{\it CP}','{\it P}^2'})
clim([-1.2,1]);
%colorbar('Ticks',[-1,-0.5,0,0.5,1]);
nexttile
imagesc(Pmat')
title('Acoustic Power')
set(gca,'yTick',1:4,'yticklabel',{'NF','NM','OF','OM'},'xTick',1:10,...
    'xticklabel',{'1','{\it N}','{\it C}','{\it P}','{\it N}^2','{\it NC}', '{\it NP}','{\it C}^2','{\it CP}','{\it P}^2'})
clim([-1.2,1]);
cb = colorbar('Ticks',[-1,-0.5,0,0.5,1]);
cb.Layout.Tile = 'east';

