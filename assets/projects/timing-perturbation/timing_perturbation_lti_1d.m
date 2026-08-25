function results = timing_perturbation_lti_1d(options)
%TIMING_PERTURBATION_LTI_1D Reproduce the scalar LTI timing-perturbation demo.
%
% This file reproduces the one-dimensional example used by the website:
%
%       x_dot(t) = A*x(t) + B*u(t),    A = 1, B = 1, x(0) = 1.
%
% The nominal piecewise-constant input is recomputed from the same convex
% finite-horizon problem used in timing-perturbation-demo-v2.html. At an
% actual (delayed) update gamma(k), the runtime input is
%
%   u(k) = u_n(k) - K(Delta(k))*(x(gamma(k)) - x_n(gamma(k))),
%
%   K(Delta(k)) = (A/B)*exp(A*T_k)/(exp(A*T_k) - 1),
%   T_k = gamma_tilde(k+1) - gamma(k).
%
% All state propagation uses the exact zero-order-hold solution. The plot
% sampling interval therefore changes only the visual smoothness, not the
% simulated update-time values.
%
% Basic use (plays an animation suitable for screen recording):
%   timing_perturbation_lti_1d
%
% Directly create an MP4 and a final PNG:
%   timing_perturbation_lti_1d(MakeVideo=true, SaveFigure=true)
%
% Run without animation (useful for checking the numbers):
%   r = timing_perturbation_lti_1d(PlayAnimation=false)
%
% Reference:
% R. Wang and N. Yao, "Robust Model Predictive Control for Networked
% Control Systems with Timing Perturbations," ACC, 2024.

arguments
    options.PlayAnimation (1,1) logical = true
    options.MakeVideo (1,1) logical = false
    options.SaveFigure (1,1) logical = true
    options.ShowFigure (1,1) logical = true
    options.PlaybackDuration (1,1) double {mustBePositive} = 15
    options.FrameRate (1,1) double {mustBePositive} = 30
    options.SampleTime (1,1) double {mustBePositive} = 0.002
    options.VideoFile (1,1) string = "timing_perturbation_lti_1d.mp4"
    options.FigureFile (1,1) string = "timing_perturbation_lti_1d.png"
    options.Deltas (1,:) double = [0, 0.097, 0.075, 0.069, 0.029, 0.026, 0, 0]
end

%% Paper system and the schedule shown in the website video
A = 1;
B = 1;
x0 = 1;
nominalTimes = [0, 1, 1.3, 2.5, 3.5, 4.35, 5.3, 6];
horizon = nominalTimes(end);
nInputs = numel(nominalTimes) - 1;
deltas = options.Deltas;

if numel(deltas) ~= numel(nominalTimes)
    error("timing_perturbation_lti_1d:DeltaCount", ...
        "Deltas must contain %d entries, including t=0 and the horizon endpoint.", ...
        numel(nominalTimes));
end
if any(deltas < 0)
    error("timing_perturbation_lti_1d:NegativeDelta", ...
        "Timing perturbations must be nonnegative.");
end

% Match Condition 1 in the web demo: exclude initialization and the added
% terminal plotting endpoint when finding the minimum scheduled interval.
internalIntervals = diff(nominalTimes(2:end-1));
minimumInterval = min(internalIntervals);
deltaBound = minimumInterval - log((exp(A*minimumInterval) + 1)/2)/A;

if any(deltas(2:nInputs) >= deltaBound)
    error("timing_perturbation_lti_1d:BoundViolation", ...
        "Each runtime perturbation must satisfy Delta(k) < %.8f s.", deltaBound);
end

actualTimes = nominalTimes + deltas;
if any(actualTimes(1:nInputs) >= nominalTimes(2:end))
    error("timing_perturbation_lti_1d:CausalityViolation", ...
        "Every delayed update gamma(k) must occur before gamma_tilde(k+1).");
end
if any(diff(actualTimes) <= 0)
    error("timing_perturbation_lti_1d:EventOrder", ...
        "Actual update times must be strictly increasing.");
end

%% Nominal MPC plan and exact runtime correction
nominalInputs = solveNominalMpc(A, B, nominalTimes, x0);
nominalEventStates = eventStates(nominalTimes, nominalInputs, A, B, x0);

actualInputs = zeros(1, nInputs);
actualEventStates = zeros(1, nInputs + 1);
nominalAtActual = zeros(1, nInputs);
eventErrors = zeros(1, nInputs);
gains = zeros(1, nInputs);

actualInputs(1) = nominalInputs(1);
actualEventStates(1) = x0;

for k = 1:nInputs
    remainingTime = nominalTimes(k+1) - actualTimes(k);
    gains(k) = dynamicGain(A, B, remainingTime);
    nominalAtActual(k) = evaluateZoh(actualTimes(k), nominalTimes, ...
        nominalInputs, A, B, nominalEventStates);

    if k > 1
        eventErrors(k) = actualEventStates(k) - nominalAtActual(k);
        actualInputs(k) = nominalInputs(k) - gains(k)*eventErrors(k);
    end

    actualEventStates(k+1) = zohStep(actualEventStates(k), ...
        actualInputs(k), actualTimes(k+1) - actualTimes(k), A, B);
end

% A plotting grid that contains every update time exactly.
time = unique([0:options.SampleTime:horizon, nominalTimes, actualTimes]);
nominalState = evaluateZoh(time, nominalTimes, nominalInputs, A, B, nominalEventStates);
actualState = evaluateZoh(time, actualTimes, actualInputs, A, B, actualEventStates);
nominalInput = evaluateInput(time, nominalTimes, nominalInputs);
actualInput = evaluateInput(time, actualTimes, actualInputs);

nominalAtSchedule = evaluateZoh(nominalTimes, nominalTimes, ...
    nominalInputs, A, B, nominalEventStates);
actualAtSchedule = evaluateZoh(nominalTimes, actualTimes, ...
    actualInputs, A, B, actualEventStates);
alignmentError = actualAtSchedule - nominalAtSchedule;

eventTable = table((0:nInputs-1).', nominalTimes(1:nInputs).', ...
    deltas(1:nInputs).', actualTimes(1:nInputs).', gains.', ...
    eventErrors.', nominalInputs.', actualInputs.', ...
    'VariableNames', {'k','NominalTime','Delta','ActualTime','DynamicGain', ...
    'StateErrorAtUpdate','NominalInput','CorrectedInput'});

fprintf('\nOne-dimensional timing-perturbation reproduction\n');
fprintf('System: x_dot = %.3g x + %.3g u,  x(0) = %.3g\n', A, B, x0);
fprintf('Admissible bound: Delta(k) < %.8f s\n', deltaBound);
fprintf('Maximum |x(gamma_tilde(k))-x_n(gamma_tilde(k))|: %.3e\n\n', ...
    max(abs(alignmentError)));
disp(eventTable);

%% Publication-style animation figure
plotData = struct( ...
    'A', A, 'B', B, 'x0', x0, 'horizon', horizon, ...
    'time', time, 'nominalTimes', nominalTimes, 'actualTimes', actualTimes, ...
    'deltas', deltas, 'deltaBound', deltaBound, ...
    'nominalInputs', nominalInputs, 'actualInputs', actualInputs, ...
    'nominalState', nominalState, 'actualState', actualState, ...
    'nominalEventStates', nominalEventStates, ...
    'actualEventStates', actualEventStates, ...
    'gains', gains, 'eventErrors', eventErrors);

[fig, graphics] = createAnimationFigure(plotData, options.ShowFigure);

thisFolder = fileparts(mfilename('fullpath'));
videoPath = resolveOutputPath(options.VideoFile, thisFolder);
figurePath = resolveOutputPath(options.FigureFile, thisFolder);

makeFrames = options.PlayAnimation || options.MakeVideo;
if makeFrames
    frameCount = max(2, round(options.PlaybackDuration*options.FrameRate));
    animationTimes = linspace(0, horizon, frameCount);
else
    animationTimes = horizon;
end

writer = [];
if options.MakeVideo
    writer = VideoWriter(videoPath, 'MPEG-4');
    writer.FrameRate = options.FrameRate;
    writer.Quality = 100;
    open(writer);
end

wallClock = tic;
try
    for frameIndex = 1:numel(animationTimes)
        updateAnimation(graphics, plotData, animationTimes(frameIndex));
        drawnow;

        if options.MakeVideo
            writeVideo(writer, getframe(fig));
        elseif options.PlayAnimation && frameIndex < numel(animationTimes)
            targetElapsed = frameIndex/options.FrameRate;
            pause(max(0, targetElapsed - toc(wallClock)));
        end
    end
catch animationError
    if ~isempty(writer)
        close(writer);
    end
    rethrow(animationError);
end

if ~isempty(writer)
    close(writer);
    fprintf('Video written to: %s\n', videoPath);
end

if options.SaveFigure
    exportgraphics(fig, figurePath, 'Resolution', 220);
    fprintf('Final figure written to: %s\n', figurePath);
end

if ~options.ShowFigure
    close(fig);
end

results = struct( ...
    'A', A, 'B', B, 'x0', x0, ...
    'nominalTimes', nominalTimes, 'actualTimes', actualTimes, ...
    'deltas', deltas, 'deltaBound', deltaBound, ...
    'time', time, 'nominalState', nominalState, 'actualState', actualState, ...
    'nominalInput', nominalInput, 'actualInput', actualInput, ...
    'nominalInputs', nominalInputs, 'correctedInputs', actualInputs, ...
    'gains', gains, 'eventErrors', eventErrors, ...
    'alignmentError', alignmentError, 'eventTable', eventTable, ...
    'figureFile', figurePath, 'videoFile', videoPath);
end


function input = solveNominalMpc(A, B, times, x0)
% Coordinate descent for the same strictly convex, box-constrained problem
% used in the HTML demo. u(0) is held at zero and -3 <= u(k) <= 3.
nInputs = numel(times) - 1;
constantState = zeros(nInputs + 1, 1);
sensitivity = zeros(nInputs + 1, nInputs);
constantState(1) = x0;

for k = 1:nInputs
    interval = times(k+1) - times(k);
    transition = exp(A*interval);
    if abs(A) < 1e-12
        inputMap = B*interval;
    else
        inputMap = ((transition - 1)/A)*B;
    end
    constantState(k+1) = transition*constantState(k);
    sensitivity(k+1,:) = transition*sensitivity(k,:);
    sensitivity(k+1,k) = sensitivity(k+1,k) + inputMap;
end

input = zeros(nInputs, 1);
state = constantState;
stateWeights = 0.5*ones(nInputs + 1, 1);
stateWeights(end) = 1;
inputWeightGradient = 1e-4;

for iteration = 1:2500
    maximumChange = 0;
    for control = 2:nInputs
        coefficient = sensitivity(:,control);
        gradient = inputWeightGradient*input(control) + ...
            2*coefficient.'*(stateWeights.*state);
        curvature = inputWeightGradient + ...
            2*sum(stateWeights.*coefficient.^2);
        nextInput = min(3, max(-3, input(control) - gradient/curvature));
        change = nextInput - input(control);
        if change ~= 0
            input(control) = nextInput;
            state = state + coefficient*change;
            maximumChange = max(maximumChange, abs(change));
        end
    end
    if maximumChange < 1e-10
        break;
    end
end

input = input.';
end


function gain = dynamicGain(A, B, remainingTime)
if remainingTime <= 0
    error("timing_perturbation_lti_1d:InvalidRemainingTime", ...
        "The dynamic-gain interval must be positive.");
end
if abs(A) < 1e-12
    gain = 1/(B*remainingTime);
else
    exponential = exp(A*remainingTime);
    gain = (A/B)*exponential/(exponential - 1);
end
end


function nextState = zohStep(state, input, step, A, B)
if abs(A) < 1e-12
    nextState = state + B*input*step;
else
    exponential = exp(A*step);
    nextState = exponential*state + ((exponential - 1)/A)*B*input;
end
end


function states = eventStates(times, inputs, A, B, x0)
states = zeros(1, numel(times));
states(1) = x0;
for k = 1:numel(inputs)
    states(k+1) = zohStep(states(k), inputs(k), times(k+1)-times(k), A, B);
end
end


function state = evaluateZoh(queryTimes, eventTimes, inputs, A, B, statesAtEvents)
state = zeros(size(queryTimes));
for q = 1:numel(queryTimes)
    segment = find(eventTimes(1:end-1) <= queryTimes(q) + 1e-12, 1, 'last');
    segment = min(segment, numel(inputs));
    state(q) = zohStep(statesAtEvents(segment), inputs(segment), ...
        queryTimes(q)-eventTimes(segment), A, B);
end
end


function input = evaluateInput(queryTimes, eventTimes, inputs)
input = zeros(size(queryTimes));
for q = 1:numel(queryTimes)
    segment = find(eventTimes(1:end-1) <= queryTimes(q) + 1e-12, 1, 'last');
    segment = min(segment, numel(inputs));
    input(q) = inputs(segment);
end
end


function [fig, graphics] = createAnimationFigure(data, showFigure)
black = [0.08, 0.10, 0.13];
red = [0.76, 0.15, 0.15];
paper = [0.985, 0.98, 0.965];
gridColor = [0.84, 0.86, 0.88];

visibility = "off";
if showFigure
    visibility = "on";
end
fig = figure('Color', paper, 'Position', [80, 80, 1280, 720], ...
    'Name', '1-D LTI Timing-Perturbation Compensation', ...
    'NumberTitle', 'off', 'Visible', visibility);
layout = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

axSchedule = nexttile(layout, 1);
hold(axSchedule, 'on');
plot(axSchedule, [0, data.horizon], [0.68, 0.68], '-', 'Color', black, 'LineWidth', 1.2);
plot(axSchedule, [0, data.horizon], [0.28, 0.28], '-', 'Color', [0.45, 0.49, 0.53], 'LineWidth', 1.0);

nEvents = numel(data.nominalTimes);
actualTicks = gobjects(1, nEvents);
actualLabels = gobjects(1, nEvents);
delayPatches = gobjects(1, nEvents);
for k = 1:nEvents
    nominalTime = data.nominalTimes(k);
    plot(axSchedule, [nominalTime, nominalTime], [0.61, 0.75], '-', ...
        'Color', black, 'LineWidth', 1.5);
    text(axSchedule, nominalTime, 0.81, sprintf('%.3g', nominalTime), ...
        'HorizontalAlignment', 'center', 'Color', black, 'FontSize', 9);

    if data.actualTimes(k) > nominalTime
        delayPatches(k) = patch(axSchedule, ...
            [nominalTime, data.actualTimes(k), data.actualTimes(k), nominalTime], ...
            [0.38, 0.38, 0.57, 0.57], red, ...
            'FaceAlpha', 0.12, 'EdgeColor', 'none', 'Visible', 'off');
    else
        delayPatches(k) = patch(axSchedule, nan, nan, red, ...
            'FaceAlpha', 0, 'EdgeColor', 'none', 'Visible', 'off');
    end
    actualTicks(k) = plot(axSchedule, [data.actualTimes(k), data.actualTimes(k)], ...
        [0.21, 0.35], '-', 'Color', red, 'LineWidth', 1.5, 'Visible', 'off');
    actualLabels(k) = text(axSchedule, data.actualTimes(k), 0.12, ...
        sprintf('%.3f', data.actualTimes(k)), 'HorizontalAlignment', 'center', ...
        'Color', red, 'FontSize', 9, 'Visible', 'off');
end

scheduleDot = plot(axSchedule, data.actualTimes(1), 0.28, 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', red, 'MarkerEdgeColor', red);
xlim(axSchedule, [0, data.horizon]);
ylim(axSchedule, [0, 1]);
yticks(axSchedule, [0.28, 0.68]);
yticklabels(axSchedule, {'actual \gamma[k]', 'nominal \tilde{\gamma}[k]'});
axSchedule.TickLabelInterpreter = 'latex';
axSchedule.XTick = [];
axSchedule.Box = 'off';
axSchedule.Color = paper;
title(axSchedule, { ...
    'Timing-Perturbation Compensation for a 1-D LTI System', ...
    sprintf('x-dot = x + u,   x(0) = 1,   admissible bound: Delta(k) < %.4f s', data.deltaBound)}, ...
    'FontWeight', 'bold');

axState = nexttile(layout, 2, [2, 1]);
hold(axState, 'on');
yline(axState, 0, '-', 'Target  x = 0', 'Color', [0.40, 0.43, 0.46], ...
    'LabelHorizontalAlignment', 'right');
nominalStateLine = plot(axState, data.time, data.nominalState, '-', ...
    'Color', black, 'LineWidth', 2.2, 'DisplayName', 'Nominal / precomputed');
actualStateLine = plot(axState, data.time(1), data.actualState(1), '--', ...
    'Color', red, 'LineWidth', 2.2, 'DisplayName', 'Perturbed / runtime');
nominalDot = plot(axState, 0, data.x0, 'o', 'MarkerSize', 7, ...
    'MarkerFaceColor', black, 'MarkerEdgeColor', paper, 'HandleVisibility', 'off');
actualDot = plot(axState, 0, data.x0, 'o', 'MarkerSize', 7, ...
    'MarkerFaceColor', red, 'MarkerEdgeColor', paper, 'HandleVisibility', 'off');
stateCursor = plot(axState, [0, 0], [-1, 1], ':', ...
    'Color', [0.40, 0.43, 0.46], 'LineWidth', 1.0, 'HandleVisibility', 'off');
xlim(axState, [0, data.horizon]);
stateLimits = paddedLimits([data.nominalState, data.actualState], 0.12);
ylim(axState, stateLimits);
stateCursor.YData = stateLimits;
ylabel(axState, 'State  x(t)');
title(axState, 'Exact zero-order-hold state trajectory', 'FontWeight', 'normal');
legend(axState, [nominalStateLine, actualStateLine], 'Location', 'northeast', 'Box', 'off');
styleAxes(axState, paper, gridColor);

axControl = nexttile(layout, 4);
hold(axControl, 'on');
[nominalStepTime, nominalStepInput] = stepCoordinates(data.nominalTimes, data.nominalInputs);
[actualStepTime, actualStepInput] = stepCoordinates(data.actualTimes, data.actualInputs);
plot(axControl, nominalStepTime, nominalStepInput, '-', ...
    'Color', black, 'LineWidth', 2.0, 'DisplayName', 'Nominal');
actualControlLine = plot(axControl, actualStepTime(1), actualStepInput(1), '--', ...
    'Color', red, 'LineWidth', 2.0, 'DisplayName', 'Corrected');
controlCursor = plot(axControl, [0, 0], [-1, 1], ':', ...
    'Color', [0.40, 0.43, 0.46], 'LineWidth', 1.0, 'HandleVisibility', 'off');
xlim(axControl, [0, data.horizon]);
controlLimits = paddedLimits([data.nominalInputs, data.actualInputs], 0.18);
ylim(axControl, controlLimits);
controlCursor.YData = controlLimits;
xlabel(axControl, 'Time (s)');
ylabel(axControl, 'Input  u(t)');
title(axControl, 'Piecewise-constant control input', 'FontWeight', 'normal');
styleAxes(axControl, paper, gridColor);

graphics = struct( ...
    'figure', fig, 'scheduleAxes', axSchedule, 'stateAxes', axState, ...
    'controlAxes', axControl, 'actualTicks', actualTicks, ...
    'actualLabels', actualLabels, 'delayPatches', delayPatches, ...
    'scheduleDot', scheduleDot, 'actualStateLine', actualStateLine, ...
    'nominalDot', nominalDot, 'actualDot', actualDot, ...
    'stateCursor', stateCursor, 'actualControlLine', actualControlLine, ...
    'controlCursor', controlCursor, 'actualStepTime', actualStepTime, ...
    'actualStepInput', actualStepInput);
end


function updateAnimation(graphics, data, currentTime)
stateMask = data.time <= currentTime + 1e-12;
stateTime = data.time(stateMask);
actualState = data.actualState(stateMask);

if stateTime(end) < currentTime
    stateTime(end+1) = currentTime;
    actualState(end+1) = evaluateZoh(currentTime, data.actualTimes, ...
        data.actualInputs, data.A, data.B, data.actualEventStates);
end

graphics.actualStateLine.XData = stateTime;
graphics.actualStateLine.YData = actualState;

nominalNow = evaluateZoh(currentTime, data.nominalTimes, ...
    data.nominalInputs, data.A, data.B, data.nominalEventStates);
actualNow = evaluateZoh(currentTime, data.actualTimes, ...
    data.actualInputs, data.A, data.B, data.actualEventStates);
graphics.nominalDot.XData = currentTime;
graphics.nominalDot.YData = nominalNow;
graphics.actualDot.XData = currentTime;
graphics.actualDot.YData = actualNow;
graphics.stateCursor.XData = [currentTime, currentTime];
graphics.controlCursor.XData = [currentTime, currentTime];

[stepTime, stepInput] = clipStep(graphics.actualStepTime, ...
    graphics.actualStepInput, currentTime);
graphics.actualControlLine.XData = stepTime;
graphics.actualControlLine.YData = stepInput;

revealed = find(data.actualTimes <= currentTime + 1e-12);
for k = 1:numel(data.actualTimes)
    isVisible = any(revealed == k);
    graphics.actualTicks(k).Visible = onOff(isVisible);
    graphics.actualLabels(k).Visible = onOff(isVisible);
    if isgraphics(graphics.delayPatches(k)) && data.deltas(k) > 0
        graphics.delayPatches(k).Visible = onOff(isVisible);
    end
end

latest = revealed(end);
graphics.scheduleDot.XData = data.actualTimes(latest);
if latest == 1
    eventText = sprintf('t = %.2f s  |  nominal MPC plan available; waiting for the first delayed update', currentTime);
elseif latest > numel(data.gains)
    eventText = sprintf(['t = %.2f s  |  horizon reached:  ' ...
        'x = %.3e,  x_n = %.3e,  alignment error = %.3e'], ...
        currentTime, actualNow, nominalNow, actualNow-nominalNow);
else
    eventText = sprintf(['t = %.2f s  |  latest update k = %d:  Delta = %.3f s,  ' ...
        'K = %.3f,  e = %+.3f,  u = %+.3f'], currentTime, latest-1, ...
        data.deltas(latest), data.gains(latest), data.eventErrors(latest), ...
        data.actualInputs(latest));
end
graphics.scheduleAxes.Title.String{3} = eventText;
end


function [stepTime, stepValue] = stepCoordinates(eventTimes, inputs)
nInputs = numel(inputs);
stepTime = zeros(1, 2*nInputs);
stepValue = zeros(1, 2*nInputs);
for k = 1:nInputs
    indices = 2*k-1:2*k;
    stepTime(indices) = eventTimes(k:k+1);
    stepValue(indices) = inputs(k);
end
end


function [clippedTime, clippedValue] = clipStep(stepTime, stepValue, limit)
mask = stepTime <= limit + 1e-12;
clippedTime = stepTime(mask);
clippedValue = stepValue(mask);
if isempty(clippedTime)
    clippedTime = 0;
    clippedValue = stepValue(1);
elseif clippedTime(end) < limit
    segment = find(stepTime <= limit, 1, 'last');
    clippedTime(end+1) = limit;
    clippedValue(end+1) = stepValue(segment);
end
end


function limits = paddedLimits(values, fraction)
lower = min(values);
upper = max(values);
span = upper - lower;
if span < 1e-9
    span = max(1, abs(upper));
end
limits = [lower - fraction*span, upper + fraction*span];
end


function styleAxes(ax, background, gridColor)
ax.Color = background;
ax.Box = 'off';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridColor = gridColor;
ax.GridAlpha = 0.65;
ax.FontName = 'Times New Roman';
ax.FontSize = 11;
end


function value = onOff(condition)
if condition
    value = 'on';
else
    value = 'off';
end
end


function outputPath = resolveOutputPath(requestedPath, baseFolder)
[folder, name, extension] = fileparts(char(requestedPath));
if isempty(folder)
    outputPath = fullfile(baseFolder, [name, extension]);
else
    outputPath = char(requestedPath);
end
end
