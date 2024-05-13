

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

% Vector dimension:
vdim = size(tres);

mdim = vdim(mam);
sdim = vdim(scm);

% indexes
for i=1:length(vdim)
    aidx{i} = 1:vdim(i); 
end

%%


% Max Temp
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

switch mam
    case test_pos
        %TODO this
        xt = round(median(init_cond)-273.15);
        xmaml = split(int2str(xt), '  ');
    case dom_pos
        xmaml = {'1D', '4D', '9D', 'AD'};
    case mod_pos
        xmaml = {'Water', 'Air', 'Rack'};
    case alg_pos
    case wl_pos
        xmaml = {'MAX-WL', 'MULTI-WL', 'CLOUD-WL'};
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
    case alg_pos
        xsaml = {'FCA', 'EBA', 'VBA'};
    case wl_pos
        xsaml = {'MAX-WL', 'MULTI-WL', 'CLOUD-WL'};
end

x = 1:mdim;

t.TileSpacing = 'compact';
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
		    'VerticalAlignment','bottom');
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
Tit = [Tit title('Max Exceeding Temperature Value')];
xticklabels(xmaml);
ylabel('[°C]');

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
		    'VerticalAlignment','bottom');
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
Tit = [Tit title('Exceeding Temperature Time Percentage')];
xticklabels(xmaml);
ylabel('[%]');


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
		    'VerticalAlignment','bottom');
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
Tit = [Tit title('Ratio between Exceeding Power and Power Budget')];
xticklabels(xmaml);
ylabel('[%]');

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
		    'VerticalAlignment','bottom');
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
Tit = [Tit title('Exceeding Power Time Percentage')];
xticklabels(xmaml);
ylabel('[%]');

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
		    'VerticalAlignment','bottom');
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
Tit = [Tit title('L2-Norm of Frequency Difference \n (less is better)')];
xticklabels(xmaml);
ylabel('[]');

%%%%%%%%%%%
% WL Av
nexttile([1 4]);

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
		    'VerticalAlignment','bottom');
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
Tit = [Tit title('Cores Average Executed Workload')];
xticklabels(xmaml);
ylabel('[%]');

%%%%%%%%%%%
% WL min
nexttile([1 4]);

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
		    'VerticalAlignment','bottom');
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
Tit = [Tit title('Cores minimum Executed Workload')];
xticklabels(xmaml);
ylabel('[%]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% LEGEND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lgd = legend();
lgd.Layout.Tile = [];

%
%

%fontsize(f, "increase");
set(findall(gcf,'-property','FontSize'),'FontSize',12)


for i=1:length(Tit)
    Tit(i).FontSize = 16;
end




%%
%hist(extractCell(tres, aidx, "temp", "exMn", "MeanAv"))
%extractCell(tres, aidx, "perf", "wl", "Av")


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
