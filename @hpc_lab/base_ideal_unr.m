function [wlres, wlop] = base_ideal_unr(obj, chip, cid)
%BASE_IDEAL_UNR Summary of this function goes here
%   Detailed explanation goes here


ctrl = ideal_unr(chip);
ctrl.C = chip.Cc;
ctrl.Ts_ctrl = 500e-6;


simres = obj.simulation(ctrl, chip,1,0);

wlres = ctrl.widx;
wlop = simres.wlop;

%Clear and update stuff
obj.perf_max_check{cid} = wlres;
obj.pmc_need_update = 0;

end

