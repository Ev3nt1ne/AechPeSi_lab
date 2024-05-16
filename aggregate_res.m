
%{
colorplot = [ 0.9843    0.7059    0.6824;
    0.7020    0.8039    0.8902;
    0.8000    0.9216    0.7725;
    0.8706    0.7961    0.8941;
    0.9961    0.8510    0.6510;
    1.0000    1.0000    0.8000;
    0.8980    0.8471    0.7412;
    0.9922    0.8549    0.9255;
    0.9490    0.9490    0.9490;
];
%}

%TODO not working (Add folders and subfolders
addpath = ('Libraries/DataViz-3.2.3.0/daviolinplot')

colorplot = [0.8941    0.1020    0.1098;
    0.2157    0.4941    0.7216;
    0.3020    0.6863    0.2902;
    0.5961    0.3059    0.6392;
    1.0000    0.4980         0;
    1.0000    1.0000    0.2000;
    0.6510    0.3373    0.1569;
    0.9686    0.5059    0.7490;
    0.6000    0.6000    0.6000;];

mycolors = colorplot([2 5 3], :);
mycolors = [];


% TODO: At the moment I assumed you could not separate per algo or per test
%       (since it does not make much sense)

test_pos = 1;
dom_pos = 2;
mod_pos = 3;
alg_pos = 4;
wl_pos = 5;

% Main Aggregator metric: *test*/domain/model/*alg*/wl
mam = wl_pos;

%secondary Aggregator metric or Comparison metric (i.e. the one in the bars)
scm = alg_pos;

% show values:
show_tip_values = 1;
show_violin = 0;
font_size_tips = 10;
font_size_tit = 18;
font_size_axes = 12;

% Vector dimension:
vdim = size(tres);

mdim = vdim(mam);
sdim = vdim(scm);

% indexes
for i=1:length(vdim)
    aidx{i} = 1:vdim(i); 
end

%%


% Initialize
T_M = zeros(3, mdim,sdim);
T_AvEv = zeros(3, mdim,sdim);
T_AvEt = zeros(3, mdim,sdim);
P_AvEp = zeros(3, mdim,sdim);
P_MEp = zeros(3, mdim,sdim);
P_AvEt = zeros(3, mdim,sdim);
W_AvFD = zeros(3, mdim,sdim);
W_MWL = zeros(3, mdim,sdim);
W_AvWL = zeros(3, mdim,sdim);
W_mWL = zeros(3, mdim,sdim);


for m=1:mdim
    
    aidx{mam} = m;

    for s=1:sdim
        aidx{scm} = s;

        % Max temp
        data = (extractCell(tres, aidx, "temp", "Max"));
        T_M(1,m,s) = max(data);
        T_M(2,m,s) = mean(data);
        T_M(3,m,s) = min(data);

        % Avearge Exceeding Temp
        data = (extractCell(tres, aidx, "temp", "exMn", "MeanAv"));
        T_AvEv(2,m,s) = mean(data);
        T_AvEv(1,m,s) = max(data);
        T_AvEv(3,m,s) = min(data);

        % Average Exceeding Temp time
        Et = (extractCell(tres, aidx, "temp", "exMn", "MeanTime"));
        T_AvEt(2,m,s) = mean(Et);
        T_AvEt(1,m,s) = max(Et);
        T_AvEt(3,m,s) = min(Et);

        %%% POWER

        % Average Exceeding Power 
        data = (extractCell(tres, aidx, "power", "exAvP"));
        P_AvEp(2,m,s) = mean(data);
        P_AvEp(1,m,s) = max(data);
        P_AvEp(3,m,s) = min(data);

        % 95P Exceeding Power 
        data = (extractCell(tres, aidx, "power", "ex95P"));
        P_MEp(1,m,s) = max(data);
        P_MEp(2,m,s) = mean(data);
        P_MEp(3,m,s) = min(data);


        % Average Exceeding Power time
        Et = (extractCell(tres, aidx, "power", "exTime"));
        P_AvEt(2,m,s) = mean(Et);
        P_AvEt(1,m,s) = max(Et);
        P_AvEt(3,m,s) = min(Et);

        %%% PERFORMANCE

        % fd norm
        l2nfd = (extractCell(tres, aidx, "perf", "fd", "l2norm"));
        % MAX
        W_AvFD(1,m,s) = max(l2nfd);
        % Mean
        W_AvFD(2,m,s) = mean(l2nfd);
        % min
        W_AvFD(3,m,s) = min(l2nfd);

        % wl Max
        data = (extractCell(tres, aidx, "perf", "wl", "Max"));
        W_MWL(1,m,s) = max(data);
        W_MWL(2,m,s) = mean(data);
        W_MWL(3,m,s) = min(data);

        % wl Mean
        data = (extractCell(tres, aidx, "perf", "wl", "Av"));
        W_AvWL(2,m,s) = mean(data);
        W_AvWL(1,m,s) = max(data);
        W_AvWL(3,m,s) = min(data);

        % wl min
        data = (extractCell(tres, aidx,  "perf", "wl", "min"));
        W_mWL(3,m,s) = min(data);
        W_mWL(1,m,s) = max(data);
        W_mWL(2,m,s) = mean(data);

    end
end

%% Plot

f = figure();
t = tiledlayout(3,12);

%TODO this!
switch mam
    case test_pos
        %TODO this
        xt = round(median(init_cond)-273.15);
        xmaml = split(int2str(xt), '  ');
        title_name = "Tests Iterations";
    case dom_pos
        xmaml = {'1D', '4D', '9D', 'AD'};
        title_name = "Number of Domains";
    case mod_pos
        xmaml = {'Water', 'Air', 'Rack'};
        title_name = "Thermal Model";
    case alg_pos
    case wl_pos
        xmaml = {'MAX-WL', 'MULTI-WL', 'CLOUD-WL'};
        title_name = "Workload Type";
end

switch scm
    case test_pos
        %TODO this
        xt = round(median(init_cond)-273.15);
        xsaml = split(int2str(xt), '  ');        
    case dom_pos
        xsaml = {'1D', '4D', '9D', 'AD'};        
    case mod_pos
        xsaml = {'Water', 'Air', 'Rack'};
        title_name = "Thermal Model";
    case alg_pos
        xsaml = {'FCA', 'EBA', 'VBA'};
    case wl_pos
        xsaml = {'MAX-WL', 'MULTI-WL', 'CLOUD-WL'};
        title_name = "Workload Type";
end

x = 1:mdim;

t.TileSpacing = 'loose'; %'compact';
t.Padding = 'tight';

Tit = [];

%%%%%%%%%%%%%%%%%%%%%%
%%% TEMP
%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%
% Max Exceeding Value
nexttile([1 6]);

%TODO: FIX THIS
model_series = T_M - (mean(ctrl.T_target)-273.15);
model_series(model_series<0) = 0;

%model_series = T_AvEv;
%TODO: FIX THIS
%err_high = T_M-mean(ctrl.T_target)+273.15;
%err_low = zeros(size(T_M));

err_high = squeeze(model_series(1,:,:));
err_low = squeeze(model_series(3,:,:));
model_series = squeeze(model_series(2,:,:));

err_high = err_high - model_series;
err_low = model_series - err_low;

b = bar(model_series, 'grouped');

if show_tip_values
for bi=1:sdim
    xtips = b(bi).XEndPoints;
    ytips = b(bi).YEndPoints;
    labels = string(round(b(bi).YData,2));
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
		    'VerticalAlignment','bottom','FontSize',font_size_tips);
end
end

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);

% Get the x coordinate of the bars
xc = nan(nbars, ngroups);
for i = 1:nbars
    xc(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(xc',model_series, err_low, err_high, 'k','linestyle','none');
hold off

set(gca, 'YGrid', 'on', 'XGrid', 'off')
Tit = [Tit title('Max Exceeding Temperature','FontSize',font_size_tit)];
xticklabels(xmaml);
ylabel('[°C]','FontSize',font_size_axes);
ax = gca;
ax.XAxis.FontSize = font_size_axes;
ax.YAxis.FontSize = font_size_axes;

%%%%%%%%%%%
% Exceeding time
nexttile([1 6]);

%TODO fix this
model_series = T_AvEt/2*100;

err_high = squeeze(model_series(1,:,:));
err_low = squeeze(model_series(3,:,:));
model_series = squeeze(model_series(2,:,:));

err_high = err_high - model_series;
err_low = model_series - err_low;

b = bar(model_series, 'grouped');

if show_tip_values
for bi=1:sdim
    xtips = b(bi).XEndPoints;
    ytips = b(bi).YEndPoints;
    labels = string(round(b(bi).YData,2));
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
		    'VerticalAlignment','bottom','FontSize',font_size_tips);
end
end

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);

% Get the x coordinate of the bars
xc = nan(nbars, ngroups);
for i = 1:nbars
    xc(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(xc',model_series, err_low, err_high, 'k','linestyle','none');
hold off

set(gca, 'YGrid', 'on', 'XGrid', 'off')
Tit = [Tit title('Total Exceeding Time - Temperature','FontSize',font_size_tit)];
xticklabels(xmaml);
ylabel('[%]','FontSize',font_size_axes);
ax = gca;
ax.XAxis.FontSize = font_size_axes;
ax.YAxis.FontSize = font_size_axes;


%%%%%%%%%%%%%%%%%%%%%%
%%% POWER
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
% Exceeding Value
nexttile([1 5]);

model_series = P_AvEp*100;

err_high = squeeze(model_series(1,:,:));
err_low = squeeze(model_series(3,:,:));
model_series = squeeze(model_series(2,:,:));

err_high = err_high - model_series;
err_low = model_series - err_low;

b = bar(model_series, 'grouped');

if show_tip_values
for bi=1:sdim
    xtips = b(bi).XEndPoints;
    ytips = b(bi).YEndPoints;
    labels = string(round(b(bi).YData,2));
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
		    'VerticalAlignment','bottom','FontSize',font_size_tips);
end
end

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);

% Get the x coordinate of the bars
xc = nan(nbars, ngroups);
for i = 1:nbars
    xc(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(xc',model_series, err_low, err_high, 'k','linestyle','none');
hold off

set(gca, 'YGrid', 'on', 'XGrid', 'off')
Tit = [Tit title('Average Exceeding Power','FontSize',font_size_tit)];
xticklabels(xmaml);
ylabel('[%]','FontSize',font_size_axes);
ax = gca;
ax.XAxis.FontSize = font_size_axes;
ax.YAxis.FontSize = font_size_axes;

%%%%%%%%%%%
% Exceeding time
nexttile([1 5]);

%TODO fix this
model_series = P_AvEt/1.5*100;

err_high = squeeze(model_series(1,:,:));
err_low = squeeze(model_series(3,:,:));
model_series = squeeze(model_series(2,:,:));

err_high = err_high - model_series;
err_low = model_series - err_low;

b = bar(model_series, 'grouped');

if show_tip_values
for bi=1:sdim
    xtips = b(bi).XEndPoints;
    ytips = b(bi).YEndPoints;
    labels = string(round(b(bi).YData,2));
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
		    'VerticalAlignment','bottom','FontSize',font_size_tips);
end
end

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);

% Get the x coordinate of the bars
xc = nan(nbars, ngroups);
for i = 1:nbars
    xc(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(xc',model_series, err_low, err_high, 'k','linestyle','none');
hold off

set(gca, 'YGrid', 'on', 'XGrid', 'off')
Tit = [Tit title('Total Exceeding Time - Power','FontSize',font_size_tit)];
xticklabels(xmaml);
ylabel('[%]','FontSize',font_size_axes);
ax = gca;
ax.XAxis.FontSize = font_size_axes;
ax.YAxis.FontSize = font_size_axes;

%%%%%%%%%%%%%%%%%%%%%%
%%% PERFORMANCE
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
% Norm
nexttile([1 4]);

model_series = W_AvFD;

err_high = squeeze(model_series(1,:,:));
err_low = squeeze(model_series(3,:,:));
model_series = squeeze(model_series(2,:,:));

err_high = err_high - model_series;
err_low = model_series - err_low;

b = bar(model_series, 'grouped');

if show_tip_values
for bi=1:sdim
    xtips = b(bi).XEndPoints;
    ytips = b(bi).YEndPoints;
    labels = string(round(b(bi).YData,2));
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
		    'VerticalAlignment','bottom','FontSize',font_size_tips);
end
end

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);

% Get the x coordinate of the bars
xc = nan(nbars, ngroups);
for i = 1:nbars
    xc(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(xc',model_series, err_low, err_high, 'k','linestyle','none');
hold off

set(gca, 'YGrid', 'on', 'XGrid', 'off')
Tit = [Tit title('L2-Norm of Frequency Difference','(less is better)','FontSize',font_size_tit)];
xticklabels(xmaml);
ylabel('[GHz/s]','FontSize',font_size_axes);
ax = gca;
ax.XAxis.FontSize = font_size_axes;
ax.YAxis.FontSize = font_size_axes;

%%%%%%%%%%%
% WL Av
ax1 = nexttile([1 4]);

model_series = W_AvWL;

err_high = squeeze(model_series(1,:,:));
err_low = squeeze(model_series(3,:,:));
model_series = squeeze(model_series(2,:,:));

err_high = err_high - model_series;
err_low = model_series - err_low;

b = bar(model_series, 'grouped');

if show_tip_values
for bi=1:sdim
    xtips = b(bi).XEndPoints;
    ytips = b(bi).YEndPoints;
    labels = string(round(b(bi).YData,2));
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
		    'VerticalAlignment','bottom','FontSize',font_size_tips);
end
end


hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);

% Get the x coordinate of the bars
xc = nan(nbars, ngroups);
for i = 1:nbars
    xc(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(xc',model_series, err_low, err_high, 'k','linestyle','none');
hold off

set(gca, 'YGrid', 'on', 'XGrid', 'off')
Tit = [Tit title('Average Executed Workload','FontSize',font_size_tit)];
xticklabels(xmaml);
ylabel('[%]','FontSize',font_size_axes);
ax = gca;
ax.XAxis.FontSize = font_size_axes;
ax.YAxis.FontSize = font_size_axes;

yl1 = ylim();

%%%%%%%%%%%
% WL min
ax2 = nexttile([1 4]);

model_series = W_mWL;

err_high = squeeze(model_series(1,:,:));
err_low = squeeze(model_series(3,:,:));
model_series = squeeze(model_series(2,:,:));

err_high = err_high - model_series;
err_low = model_series - err_low;

b = bar(model_series, 'grouped');

if show_tip_values
for bi=1:sdim
    xtips = b(bi).XEndPoints;
    ytips = b(bi).YEndPoints;
    labels = string(round(b(bi).YData,2));
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
		    'VerticalAlignment','bottom','FontSize',font_size_tips);
end
end

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);

% Get the x coordinate of the bars
xc = nan(nbars, ngroups);
for i = 1:nbars
    xc(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(xc',model_series, err_low, err_high, 'k','linestyle','none');
hold off

set(gca, 'YGrid', 'on', 'XGrid', 'off')
Tit = [Tit title('Minimum Executed Workload','FontSize',font_size_tit)];
xticklabels(xmaml);
ylabel('[%]','FontSize',font_size_axes);
ax = gca;
ax.XAxis.FontSize = font_size_axes;
ax.YAxis.FontSize = font_size_axes;

yl2 = ylim();

m = min(yl1, yl2);
M = max(yl1, yl2);

ylim(ax1, [m(1) M(2)]);
ylim(ax2, [m(1) M(2)]);

%
%

%fontsize(f, "increase");
%ax = findall(gcf,'-property','FontSize');
%set(ax,'FontSize',12)


%for i=1:length(Tit)
%    Tit(i).FontSize = 18;
%end

if ~(isempty(mycolors))
    ax = findall(gcf,'-property','ColorOrder');
    set(ax, 'ColorOrder', mycolors);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% LEGEND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TODO 
lgd = legend(xsaml);
lgd.Layout.Tile = 22;
lgd.FontSize = 18;
lgd.Title.String = 'Algorithms';
lgd.Title.FontSize = 22;

%title(t,title_name,' ', 'FontSize', 28);
title(t,title_name, 'FontSize', 24);


% SAVE
% This example assumes that 'fig' is the handle to your figure. 
% In R2020a and later, you can use the exportgraphics command and specify 
%   the 'ContentType' as 'vector' (to ensure tightly cropped, scalable output):
f.Position = [1, 1, 1920*1,1080*1];
exportgraphics(f, strcat(title_name,'.pdf'), 'ContentType', 'vector');


%% VIOLINS

% reset indexes
for i=1:length(vdim)
    aidx{i} = 1:vdim(i); 
end

% initialize
T_M_d = [];
T_AvEv_d = [];
T_AvEt_d = [];
P_AvEp_d = [];
P_MEp_d = [];
P_AvEt_d = [];
W_AvFD_d = [];
W_MWL_d = [];
W_AvWL_d = [];
W_mWL_d = [];

for s=1:sdim
    aidx{scm} = s;

    % Max temp
    data = (extractCell(tres, aidx, "temp", "Max"));
    T_M_d(:,s) = data;

    % Avearge Exceeding Temp
    data = (extractCell(tres, aidx, "temp", "exMn", "MeanAv"));
    T_AvEv_d(:,s) = (data);

    % Average Exceeding Temp time
    Et = (extractCell(tres, aidx, "temp", "exMn", "MeanTime"));
    T_AvEt_d(:,s) = Et;

    %%% POWER

    % Average Exceeding Power 
    data = (extractCell(tres, aidx, "power", "exAvP"));
    P_AvEp_d(:,s) = data;

    % 95P Exceeding Power 
    data = (extractCell(tres, aidx, "power", "ex95P"));
    P_MEp_d(:,s) = data;

    % Average Exceeding Power time
    Et = (extractCell(tres, aidx, "power", "exTime"));
    P_AvEt_d(:,s) = Et;

    %%% PERFORMANCE

    % fd norm
    l2nfd = (extractCell(tres, aidx, "perf", "fd", "l2norm"));
    % MAX
    W_AvFD_d(:,s) = l2nfd;

    % wl Max
    data = (extractCell(tres, aidx, "perf", "wl", "Max"));
    W_MWL_d(:,s) = data;

    % wl Mean
    data = (extractCell(tres, aidx, "perf", "wl", "Av"));
    W_AvWL_d(:,s) = data;

    % wl min
    data = (extractCell(tres, aidx,  "perf", "wl", "min"));
    W_mWL_d(:,s) = data;
end

if (show_violin)
group_inx = ones(size(T_M_d));
for s=1:sdim
    group_inx(:,s) = s;
end

smoothing_l = 0.4; %0.75;

f2 = figure();
t2 = tiledlayout(3,12);

%TODO
xsaml = {'FCA', 'EBA', 'VBA'};

Tit = [];

%%%%%%%%%%%%%%%%%%%%%%
%%% TEMP
%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%
% Max Exceeding Value
nexttile([1 6]);

%TODO: FIX THIS
model_series = T_M_d - (mean(ctrl.T_target)-273.15);
model_series(model_series<0) = 0;

daviolinplot(model_series(:),'groups',group_inx(:),'smoothing',smoothing_l,'xtlabels',xsaml);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
Tit = [Tit title('Max Exceeding Temperature')];
ylabel('[°C]');

ylf = ylim();
ylf = max(ylf, [0 ylf(2)]);
%ylf = min(ylf, [ylf(1) 100]);
ylim(ylf);

%%%%%%%%%%%
% Exceeding time
nexttile([1 6]);

%TODO fix this
model_series = T_AvEt_d/2*100;

daviolinplot(model_series(:),'groups',group_inx(:),'smoothing',smoothing_l,'xtlabels',xsaml);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
Tit = [Tit title('Total Exceeding Time - Temperature')];
ylabel('[%]');

ylf = ylim();
ylf = max(ylf, [0 ylf(2)]);
ylf = min(ylf, [ylf(1) 100]);
ylim(ylf);


%%%%%%%%%%%%%%%%%%%%%%
%%% POWER
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
% Exceeding Value
nexttile([1 6]);

model_series = P_AvEp_d*100;

daviolinplot(model_series(:),'groups',group_inx(:),'smoothing',smoothing_l,'xtlabels',xsaml);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
Tit = [Tit title('Average Exceeding Power')];
ylabel('[%]');

ylf = ylim();
ylf = max(ylf, [0 ylf(2)]);
ylf = min(ylf, [ylf(1) 100]);
ylim(ylf);

%%%%%%%%%%%
% Exceeding time
nexttile([1 6]);

%TODO fix this
model_series = P_AvEt_d/1.5*100;

daviolinplot(model_series(:),'groups',group_inx(:),'smoothing',smoothing_l,'xtlabels',xsaml);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
Tit = [Tit title('Total Exceeding Time - Power')];
ylabel('[%]');

ylf = ylim();
ylf = max(ylf, [0 ylf(2)]);
ylf = min(ylf, [ylf(1) 100]);
ylim(ylf);

%%%%%%%%%%%%%%%%%%%%%%
%%% PERFORMANCE
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
% Norm
nexttile([1 4]);

model_series = W_AvFD_d;

daviolinplot(model_series(:),'groups',group_inx(:),'smoothing',smoothing_l/20,'xtlabels',xsaml);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
Tit = [Tit title('L2-Norm of Frequency Difference','(less is better)')];
ylabel('[GHz/s]');

ylf = ylim();
ylf = max(ylf, [0 ylf(2)]);
%ylf = min(ylf, [ylf(1) 100]);
ylim(ylf);

%%%%%%%%%%%
% WL Av
ax1 = nexttile([1 4]);

model_series = W_AvWL_d;

daviolinplot(model_series(:),'groups',group_inx(:),'smoothing',smoothing_l,'xtlabels',xsaml);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
Tit = [Tit title('Average Executed Workload')];
ylabel('[%]');

yl1 = ylim();

%%%%%%%%%%%
% WL min
ax2 = nexttile([1 4]);

model_series = W_mWL_d;

h = daviolinplot(model_series(:),'groups',group_inx(:),'smoothing',smoothing_l,'xtlabels',xsaml);%,'legend',xsaml);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
Tit = [Tit title('Minimum Executed Workload')];
ylabel('[%]');

yl2 = ylim();

m = min(yl1, yl2);
m = max(m, [0 m(2)]);
M = max(yl1, yl2);
M = min(M, [M(1) 100]);

ylim(ax1, [m(1) M(2)]);
ylim(ax2, [m(1) M(2)]);

%
%

%fontsize(f, "increase");
ax = findall(gcf,'-property','FontSize');
set(ax,'FontSize',12)


for i=1:length(Tit)
    Tit(i).FontSize = 18;
end

if ~(isempty(mycolors))
    ax = findall(gcf,'-property','ColorOrder');
    set(ax, 'ColorOrder', mycolors);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% LEGEND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TODO 
%{
h.lg.Layout.Tile = 22;
h.lg.FontSize = 18;
h.lg.Title.String = 'Algorithms';
h.lg.Title.FontSize = 22;
%}

% SAVE
% This example assumes that 'fig' is the handle to your figure. 
% In R2020a and later, you can use the exportgraphics command and specify 
%   the 'ContentType' as 'vector' (to ensure tightly cropped, scalable output): 
f2.Position = [1, 1, 1920*1,1080*1];
exportgraphics(f2, 'violin.pdf', 'ContentType', 'vector');

end %show_violin

%%

mean(T_M_d- (mean(ctrl.T_target)-273.15))
mean(T_AvEv_d)
mean(T_AvEt_d/2*100)
mean(P_AvEp_d*100)
mean(P_MEp_d)
mean(P_AvEt_d/1.5*100)
mean(W_AvFD_d)
mean(W_MWL_d)
mean(W_AvWL_d)
mean(W_mWL_d)

%lower target compliance
data = mean(W_AvFD_d);
(data(1) - data(2))/data(2)*100
(data(1) - data(3))/data(3)*100
%(data(2) - data(1))/data(1)*100
%(data(3) - data(1))/data(1)*100

%higher average wl
data = mean(W_AvWL_d);
(data(1) - data(2))/data(2)*100
(data(1) - data(3))/data(3)*100

%higher average wl
data = mean(W_mWL_d);
(data(1) - data(2))/data(2)*100
(data(1) - data(3))/data(3)*100

%max temp
data = max(T_M_d- (mean(ctrl.T_target)-273.15));
data(2) / data(1)
data(3) / data(1)
(data(2) - data(1))/data(1)*100
(data(3) - data(1))/data(1)*100

%lower exceeded temp time
data =mean(T_AvEt_d/2*100);
(data(1) - data(2))/data(2)*100
(data(1) - data(3))/data(3)*100
(data(2) - data(1))/data(1)*100
(data(3) - data(1))/data(1)*100


%%

function res = extractCell(inc, mask, field1, field2, field3)

    dim = size(inc);

    qq = cell2mat(inc);

    maskdim = length(mask);

    if maskdim ~= length(dim)
        error("wrong mask length");
    end

    tot = 1;
    divmaskdim = ones(maskdim,1);
    for i=1:maskdim
        tot = tot * length(mask{i});
        divmaskdim(i) = length(mask{i});
    end

    res = zeros(tot,1);
    index = 1;
    mask_idx = ones(maskdim,1);

    % k(d) holds the accumulated No. of array elements
    % up to the d'th subscript (or dimension)
    k = [1 cumprod(dim(1:end-1))];                                     

    for nf=1:tot

        %idx = sub2ind(dim, mask_idx)
        %  Compute linear indices, e.g. in 3D, with edges L1, L2, L3:
        %  idx = x1 + (x2-1)*L1 + (x3-1)*L1*L2,
        idx = 1;
        for i=1:maskdim
            idx = bsxfun( @plus, idx, (mask{i}(mask_idx(i))-1)*k(i) );  % iteratively calculate the sum shown above
        end

        if (nargin == 5) && ~isempty(field3)
            res(index) = qq(idx).(field1).(field2).(field3);
        elseif (nargin == 4) && ~isempty(field2)
            res(index) = qq(idx).(field1).(field2);
        else
            res(index) = qq(idx).(field1);
        end

        %squeeze(res(res~=0))

        index = index + 1;

        if index > tot
            break;
        end

        mask_idx(end) = mask_idx(end) + 1;
        while sum(mask_idx > divmaskdim)
            tr = mask_idx > divmaskdim;
            [~,pp] = max(tr);
            mask_idx(pp) = 1;
            mask_idx(pp-1) = mask_idx(pp-1)+1;
        end
    end
    

end
