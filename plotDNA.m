% $Id: plotDNA.m,v 1.5 2014-01-10 17:36:29 schowell Exp $
% script to plot the results of calculating the Re and Rg of DNA
%
% data columns: moves  rg/lp      re/lp           a       r
clc;clc;
% close all;

d.files = {'140101_164355_3lp',...
    '140106_102302_6lp',...
    '140106_113435_10lp',...
    '140106_113500_1lp',...
    '140106_114806_15lp',...
    };

d.name = {'lp003',...
    'lp006',...
    'lp010',...
    'lp001',...
    'lp015',...
    };
    
d.L = [003,006,010,001,015];

d.nsteps = [100000,100000,100000,100000,100000];

expRg = zeros(length(d.L),3);
expRe = zeros(length(d.L),3);

d.path = './validate/';
d.ext = '_dnaMoves.o';

% rg = loadxy('extractedRg.csv');
% re = loadxy('extractedRe.csv');
xVal = [0.7, 1700];
llp = logspace(log(xVal(1))/log(10),log(xVal(2))/log(10));
Rg = sqrt(llp/3-1+2./llp-2./llp.^2.*(1-exp(-llp)));
rg = [llp; Rg]';
Re = sqrt(2*(llp-1+exp(-llp)));
re = [llp; Re]';

[d.L,j] = sort(d.L);
d.files = d.files(j);
d.name = d.name(j);

n.lag = 10;
n.equil = 1e6;

for i = 1:length(d.name)
    d.(d.name{i}).cmp = spline(rg(:,1),rg(:,2),d.L(i));
    d.(d.name{i}).all = loadxy([d.path,d.files{i},d.ext]);
    d.(d.name{i}).equilI = d.(d.name{i}).all(:,1)>n.equil;
    fprintf([d.name{i},', averaged over %d iterations of %d\n'],length(d.(d.name{i}).equilI),d.nsteps(i));
    
    d.(d.name{i}).rg = d.(d.name{i}).all(d.(d.name{i}).equilI,[1,2]);
    d.(d.name{i}).re = d.(d.name{i}).all(d.(d.name{i}).equilI,[1,3]);
    n.kept = length(d.(d.name{i}).rg);
    
    d.(d.name{i}).rgA = mean(d.(d.name{i}).rg(:,2));
    d.(d.name{i}).reA = mean(d.(d.name{i}).re(:,2));
    d.(d.name{i}).rgStD = std(d.(d.name{i}).rg(:,2));
    d.(d.name{i}).reStD = std(d.(d.name{i}).re(:,2));
    d.(d.name{i}).rgStE = d.(d.name{i}).rgStD/sqrt(n.kept);
    d.(d.name{i}).reStE = d.(d.name{i}).reStD/sqrt(n.kept);
    expRg(i,1) = d.L(i);
    expRg(i,2) = d.(d.name{i}).rgA;
    expRg(i,3) = d.(d.name{i}).rgStE;
    expRe(i,1) = d.L(i);
    expRe(i,2) = d.(d.name{i}).reA;
    expRe(i,3) = d.(d.name{i}).reStE;
    
    d.(d.name{i}).a = d.(d.name{i}).all(:,4);
    d.(d.name{i}).iter = sum(sum(d.(d.name{i}).all(:,[4,5])));
    d.(d.name{i}).accep = sum(d.(d.name{i}).a)/d.(d.name{i}).iter;
    
    d.(d.name{i}).mav = d.(d.name{i}).rg;
    
    n.steps = length(d.(d.name{i}).rg);
    if n.lag > n.steps
        d.(d.name{i}).mav(:,2) = tsmovavg(d.(d.name{i}).rg(:,2),'s',ceil(n.steps/4),1);
    else
        d.(d.name{i}).mav(:,2) = tsmovavg(d.(d.name{i}).rg(:,2),'s',n.lag,1);
    end
    mySubplot;
    hold all;
%     xyplot(rg)
    xyplot(d.(d.name{i}).rg,'s')
    plot(d.(d.name{i}).rg(:,1),d.(d.name{i}).cmp*ones(length(d.(d.name{i}).rg),1),'-','linewidth',2);
    xyplot(d.(d.name{i}).mav,'.','linewidth',2)
    title(['Rg of ',d.name{i},'  (a=',num2str(d.(d.name{i}).accep),')']);
    xlabel('iterations')
    ylabel('Rg/Lp')
end
legend('Rg','actual Rg','Moving Av','location','NorthEastOutside')
legend boxoff

figure;
hold all;
% xVal = [0.7, 1700];
% x=logspace(log(xVal(1))/log(10),log(xVal(2))/log(10));
% plot(x,spline(rg(:,1),rg(:,2),x));
% xyplot(rg)
xyplot(rg)
xyplot(re)
xyerror(expRg,'o')
xyerror(expRe,'s')
legend('theoretical Rg/Lp','theoretical Re/Lp','computed Rg/Lp','computed Re/Lp','location','southeast')
legend boxoff
logxy
xlabel('L/Lp')
ylabel('size of chain')
title('Validation of DNA bulk properties')
axis tight
zoomout(.1)
saveps(gcf,'validateModel')