

% TODO: At the moment I assumed you could not separate per algo or per test
%       (since it does not make much sense)

test_pos = 1;
dom_pos = 2;
mod_pos = 3;
alg_pos = 4;
wl_pos = 5;

% Main Aggregator metric: *test*/domain/model/*alg*/wl
mam = test_pos;

%secondary Aggregator metric or Comparison metric (i.e. the one in the bars)
scm = alg_pos;

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
T_M = zeros(mdim,sdim);
T_AvEv = zeros(mdim,sdim);
T_AvEt = zeros(mdim,sdim);
T_MEt = zeros(mdim,sdim);
T_mEt = zeros(mdim,sdim);
P_AvEp = zeros(mdim,sdim);
P_MEp = zeros(mdim,sdim);
P_AvEt = zeros(mdim,sdim);
P_MEt = zeros(mdim,sdim);
P_mEt = zeros(mdim,sdim);
W_MFD = zeros(mdim,sdim);
W_AvFD = zeros(mdim,sdim);
W_mFD = zeros(mdim,sdim);
W_MWL = zeros(mdim,sdim);
W_AvWL = zeros(mdim,sdim);
W_mWL = zeros(mdim,sdim);


for m=1:mdim
    
    aidx{mam} = m;

    for s=1:sdim
        aidx{scm} = s;

        % Max temp
        T_M(m,s) = max(extractCell(tres, aidx, "temp", "Max"));

        % Avearge Exceeding Temp
        T_AvEv(m,s) = mean(extractCell(tres, aidx, "temp", "exMn", "MeanAv"));

        % Average Exceeding Temp time
        Et = (extractCell(tres, aidx, "temp", "exMn", "MeanTime"));
        T_AvEt(m,s) = mean(Et);
        T_MEt(m,s) = max(Et);
        T_mEt(m,s) = min(Et);

        %%% POWER

        % Average Exceeding Power 
        P_AvEp(m,s) = mean(extractCell(tres, aidx, "power", "exAvP"));

        % 95P Exceeding Power 
        P_MEp(m,s) = max(extractCell(tres, aidx, "power", "ex95P"));

        % Average Exceeding Power time
        Et = (extractCell(tres, aidx, "power", "exTime"));
        P_AvEt(m,s) = mean(Et);
        P_MEt(m,s) = max(Et);
        P_mEt(m,s) = min(Et);

        %%% PERFORMANCE

        % fd norm
        l2nfd = (extractCell(tres, aidx, "perf", "fd", "l2norm"));
        % MAX
        W_MFD(m,s) = max(l2nfd);
        % Mean
        W_AvFD(m,s) = mean(l2nfd);
        % min
        W_mFD(m,s) = min(l2nfd);

        % wl Max
        W_MWL(m,s) = max(extractCell(tres, aidx, "perf", "wl", "Max"));

        % wl Mean
        W_AvWL(m,s) = mean(extractCell(tres, aidx, "perf", "wl", "Av"));

        % wl min
        W_mWL(m,s) = min(extractCell(tres, aidx,  "perf", "wl", "min"));

    end
end

%% Plot

f = figure();
t = tiledlayout(3,2);

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

x = 1:mdim;

t.TileSpacing = 'compact';
t.Padding = 'tight';

%%%%%%%%%%%%%%%%%%%%%%
%%% TEMP
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
% Exceeding Value
nexttile;

model_series = T_AvEv;
%TODO: FIX THIS
err_high = T_M-mean(ctrl.T_target)+273.15;
err_low = zeros(size(T_M));

err_high = err_high - model_series;
err_low = model_series - err_low;

b = bar(model_series, 'grouped');
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

title('Exceeding Temperature Value');
xticklabels(xmaml);
ylabel('[°C]');

%%%%%%%%%%%
% Exceeding time
nexttile;

model_series = T_AvEt/0.75*100;
err_high = T_MEt/0.75*100;
err_low = T_mEt/0.75*100;

err_high = err_high - model_series;
err_low = model_series - err_low;

b = bar(model_series, 'grouped');
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

title('Exceeding Temperature Time Percentage');
xticklabels(xmaml);
ylabel('[%]');


%%%%%%%%%%%%%%%%%%%%%%
%%% POWER
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
% Exceeding Value
nexttile;

model_series = P_AvEp*100;
err_high = P_MEp*100;
err_low = zeros(size(P_AvEp));

err_high = err_high - model_series;
err_low = model_series - err_low;

b = bar(model_series, 'grouped');
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

title('Ratio between Exceeding Power and Power Budget');
xticklabels(xmaml);
ylabel('[%]');

%%%%%%%%%%%
% Exceeding time
nexttile;

model_series = P_AvEt/1.5*100;
err_high = P_MEt/1.5*100;
err_low = P_mEt/1.5*100;

err_high = err_high - model_series;
err_low = model_series - err_low;

b = bar(model_series, 'grouped');
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

title('Exceeding Power Time Percentage');
xticklabels(xmaml);
ylabel('[%]');


%%%%%%%%%%%%%%%%%%%%%%
%%% PERFORMANCE
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
% Norm
nexttile;

model_series = W_AvFD;
err_high = W_MFD;
err_low = W_mFD;

err_high = err_high - model_series;
err_low = model_series - err_low;

b = bar(model_series, 'grouped');
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

title('L2-Norm of Frequency Difference (min is better)');
xticklabels(xmaml);
ylabel('[]');

%%%%%%%%%%%
% WL
nexttile;

model_series = W_AvWL;
err_high = W_MWL;
err_low = W_mWL;

err_high = err_high - model_series;
err_low = model_series - err_low;

b = bar(model_series, 'grouped');
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

title('Executed Workload');
xticklabels(xmaml);
ylabel('[%]');




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
