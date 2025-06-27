% Transfer Entropy code with figure
clc; close all; clearvars;

% add path for JIDT and load data
javaaddpath('/infodynamics-dist-1.6/infodynamics.jar');
load('NF.csv') % data cells of each condition 
load('NM.csv')
load('OF.csv')
load('OM.csv')
data = {NF, NM, OF, OM};

% reshape for time step of t=1/30 sec
for n=1:4
    N = mean(reshape(data{1,n}(:,1),2,[]))';
    C = sum(reshape(data{1,n}(:,2),2,1800))'./N;
    P = sum(reshape(data{1,n}(:,3),2,1800))'./N;
    t30 = [N, C, P];
    t30(isnan(t30))=0;
    t30(isinf(t30))=0;
    data{1,n} = t30;
end

dataqu = {'NF','NM', 'OF', "OM"};
z = zeros(2,4);
col = parula(3);
ytitles = {'Calls per Bat', 'Acoustic Power per Bat'};

figure
for m =1:2 % calls, power
    subplot(1,2,m);
for j = 1:4 % for each condition (NF, NM, OF, OM)
    numObservations = size(data{1,1},1);
    sourceArray=data{1,j}(:,1); % column 1 = number of bats
    destArray = data{1,j}(:,m+1); % column 2 = calls; column 3 = power
 
    % Create a TE calculator and run it:
    teCalc=javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov');
    teCalc.setProperty('k', '6'); % Use Kraskov parameter K=6 for 6 nearest points
    teCalc.initialise(1); % Use history length 1 (Schreiber k=1)

    % Perform calculation with correlated source:
    teCalc.setObservations(sourceArray, destArray);
    result = teCalc.computeAverageLocalOfObservations();
    result(result<0)=0;
    fprintf('TE result %.4f nats \n', result);   

    % for permuted source
    resultvec = [];
    for a = 1: 100
        rand_source = [];
        for k = 1:(size(data{1,j},1)/30)
            rand_source(30*(k-1)+1:30*k)=data{1,j}(randperm(30)+(30*(k-1)));
        end
        sourceArray2=rand_source; % Uncorrelated source
        % Perform calculation with uncorrelated source:
        teCalc.initialise(); % Initialise leaving the parameters the same
        teCalc.setObservations(sourceArray2, destArray);
        result2 = teCalc.computeAverageLocalOfObservations();
        result2(result2<0)=0;
        resultvec(a) = result2;     
    end

    % z-score significance
    M = mean(resultvec);
    S = std(resultvec);
    z = (result-M)/S;
    if z>1.96 % 95% confidence
        fprintf('significant \n')
    end

    % plot swarm chart
    swarmchart(j*ones(1,100),resultvec,3, col(2,:),'filled', 'MarkerFaceAlpha',.5);
    hold on, box on
    plot(j, result,'-kp', 'markersize', 7, 'MarkerFaceColor','k');
    xticks([1,2,3,4]);
    xticklabels(dataqu);
        axis padded
        ylabel('Transfer Entropy');
        ylim([0 .3]);
        title(ytitles{m});
end
end

% size of image
x0=10;
y0=10;
width=850;
height=430;
set(gcf,'position',[x0,y0,width,height])