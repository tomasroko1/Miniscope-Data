function s = caliali_settings()
% CALIALI_SETTINGS  Settings of the caliali_* pipeline that are not CaliAli parameters.
% One place for the whole pipeline: change a value here and every function uses it. Keep these
% the same for all days of an animal (the summary of each day records the smin it used).
%
%   smin           Threshold of the final deconvolution (C and S): events smaller than
%                  abs(smin) noise units are dropped. CaliAli uses -5. Lower abs(smin) detects
%                  smaller events (and more false ones). Try values with
%                  caliali_redeconv(day, [-5 -4 -3 -2]); to apply a new value to days that are
%                  already extracted, run caliali_run(animal, 'summary') after changing it here.
%   min_events     events a neuron needs in a session to count as active there (summary)
%   contour_level  neuron contours are drawn at this fraction of the footprint's peak (summary)
s.smin          = -5;
s.min_events    = 3;
s.contour_level = 0.5;

% what CaliAli itself used for the saved extraction result (runCNMFe, line 106): do not change
s.smin_caliali  = -5;
end
