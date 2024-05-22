function [wlres, wlop] = base_ideal_unr(obj)
%BASE_IDEAL_UNR Summary of this function goes here
%   Detailed explanation goes here


ctrl = ideal_unr;
ctrl.C = obj.Cc;
ctrl.Ts_ctrl = 500e-6;


[~, ~, ~, ~, wlop] = obj.simulation(ctrl, 0);

wlres = ctrl.widx;

%Clear and update stuff
obj.perf_max_check = wlres;
obj.pmc_need_update = 0;

end

