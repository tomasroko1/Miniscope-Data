function caliali_redeconv(basedir, smin, min_events)
% CALIALI_REDECONV  Redo only the deconvolution (C and S) of an extracted day with another smin.
%   caliali_redeconv(day, [-5 -4 -3 -2])   % SCAN: events per session for each smin; saves nothing
%   caliali_redeconv(day, -3)              % apply smin = -3 and rebuild the day's summary
%   caliali_redeconv(day, -3, 3)           % ... with 'active' = at least 3 events (default: caliali_settings)
%
% What it is: CaliAli ends an extraction by deconvolving each neuron's raw trace C_raw into the
% denoised trace C and the events S, keeping only events larger than abs(smin) noise units
% (CaliAli uses smin = -5). A lower abs(smin) detects smaller events (and more false ones), a
% higher one fewer. Nothing before that step depends on it, so footprints, alignment and the
% extraction are not touched, and C_raw is not changed.
%
% Applying a value (single smin): the CNMF-E result file is NOT modified. The deconvolution is
% redone in memory with CaliAli's own postprocessDeconvolvedTraces on the saved neuron, and
% caliali_summary rewrites the day's summary (.mat/.csv/.png) and traces in CaliAli_results with
% the new C and S (summary.smin records the value). Neuron numbers stay the same as before. To go
% back: caliali_redeconv(day, -5). Scan mode (several values) only prints a table.
% To make a value the pipeline's for all days: set s.smin in caliali_settings.m and run
% caliali_run(animal, 'summary') (new extractions use it automatically). Calling caliali_summary
% for a day afterwards, or caliali_run(..., 'summary'), applies the value in caliali_settings.m.
% Run it on the server (it needs the work folder and the stripe records), in a normal MATLAB
% session: the deconvolution uses the parallel pool and takes a few minutes per value.
st = caliali_settings();
if nargin < 3 || isempty(min_events), min_events = st.min_events; end
caliali_setup();
p = caliali_day_paths(basedir);
lock = caliali_lock(p); %#ok<NASGU>
result = caliali_result_file(p);
fprintf('Result: %s\n', result);
R = load(result, 'neuron');
neuron = R.neuron;
F = double(neuron.CaliAli_options.inter_session_alignment.F(:)');      % frames per session
assert(sum(F) == size(neuron.C_raw, 2), 'Session frame counts do not add up to the trace length');
edges = [0 cumsum(F)];

if numel(smin) > 1
    fprintf('\nScan, nothing is saved. Events = onsets of S > 0; active = at least %d events in the session.\n', min_events);
    fprintf('%6s %12s %14s %s\n', 'smin', 'total events', 'neurons w/o events', 'active neurons per session');
    for v = smin(:)'
        neuron = postprocessDeconvolvedTraces(neuron, 'foopsi', 'ar2', v);   % redone from C_raw each time
        ev = events_per_session(neuron.S, edges);
        fprintf('%6g %12d %14d      %s\n', v, sum(ev(:)), sum(sum(ev, 2) == 0), mat2str(sum(ev >= min_events, 1)));
    end
    fprintf('\nPick a value and apply it: caliali_redeconv(day, <smin>)\n');
    return
end

before = events_per_session(neuron.S, edges);
fprintf('Before: %d events in total (smin as saved in the result).\n', sum(before(:)));
neuron = postprocessDeconvolvedTraces(neuron, 'foopsi', 'ar2', smin);
after = events_per_session(neuron.S, edges);
fprintf('After smin = %g: %d events in total.\n', smin, sum(after(:)));
caliali_summary(basedir, min_events, [], neuron, smin);
end

function ev = events_per_session(S, edges)
S = full(S);
K = size(S, 1);
ev = zeros(K, numel(edges) - 1);
for s = 1:numel(edges) - 1
    on = S(:, edges(s)+1:edges(s+1)) > 0;
    ev(:, s) = sum(diff([zeros(K, 1) on], 1, 2) == 1, 2);
end
end
