function results = timing_perturbation_hydrotank_2d(options)
%TIMING_PERTURBATION_HYDROTANK_2D Animate the two-tank timing example.
%
% The model, schedule, and nominal MPC inputs are reproduced from the
% Hydro_tank/MPC_simu MATLAB files. The source script's signed timing values
% are not used: this implementation enforces 0 <= Delta(k) < Delta_bar.
% It also replaces forward Euler integration with the exact ZOH transition.
%
% Screen-record the MATLAB animation:
%   timing_perturbation_hydrotank_2d
%
% Create an MP4 directly:
%   timing_perturbation_hydrotank_2d(MakeVideo=true)
%
% Run the numerical check without playing the animation:
%   r = timing_perturbation_hydrotank_2d(PlayAnimation=false)

arguments
    options.PlayAnimation (1,1) logical = true
    options.MakeVideo (1,1) logical = false
    options.SaveFigure (1,1) logical = true
    options.ShowFigure (1,1) logical = true
    options.PlaybackDuration (1,1) double {mustBePositive} = 18
    options.FrameRate (1,1) double {mustBePositive} = 30
    options.SampleTime (1,1) double {mustBePositive} = 0.05
    options.VideoFile (1,1) string = "timing_perturbation_hydrotank_2d.mp4"
    options.FigureFile (1,1) string = "timing_perturbation_hydrotank_2d.png"
end

%% Source model and experiment data
A = [-0.0146, 0.0146; 0.0146, -0.0272];
B = [64.9351, 0; 0, 64.9351];
x0 = [0.1; 0.1];
setpoint = [0.4; 0.23];
horizon = 240;

% The first input is active at t=0. The remaining inputs are scheduled at
% Time(k); the final value 245 s is the end of the last ZOH interval.
scheduledTimes = [45, 60, 110, 130, 150, 175, 200, 245];
deltas = [3, 3.5, 3, 3, 3, 4, 3, 0];
nominalEvents = [0, scheduledTimes];
actualEvents = [0, scheduledTimes + deltas];

u1 = [0.000120009988644764, 9.00888475692411e-08, ...
      2.68462944561538e-05, 5.00088595363246e-05, ...
      3.49203125633098e-05, 3.89838601619980e-05, ...
      3.79749413781706e-05, 3.82777294815020e-05];
u2 = [0.000118898849526939, 9.00153373236349e-08, ...
      9.00492773566249e-08, 3.73114908157035e-06, ...
      7.06310764084674e-06, 6.28908224803212e-06, ...
      6.42811457975071e-06, 6.40711365534216e-06];
nominalInputs = [u1; u2];

% Both eigenvalues of A are negative. For this case, Condition 1 gives
% Delta_bar = Delta_t_min/2. The shortest scheduled interval is 15 s.
minimumInterval = min(diff(nominalEvents));
deltaBound = minimumInterval/2;
if any(deltas < 0)
    error("timing_perturbation_hydrotank_2d:NegativeDelta", ...
        "Timing perturbations must be nonnegative.");
end
if any(deltas(1:end-1) >= deltaBound)
    error("timing_perturbation_hydrotank_2d:BoundViolation", ...
        "Each timing perturbation must satisfy Delta(k) < %.6g s.", deltaBound);
end
if any(actualEvents(2:end-1) >= nominalEvents(3:end))
    error("timing_perturbation_hydrotank_2d:CausalityViolation", ...
        "Each actual update must occur before the next scheduled update.");
end
if any(diff(actualEvents) <= 0)
    error("timing_perturbation_hydrotank_2d:EventOrder", ...
        "The perturbed update times must be strictly increasing.");
end

%% Reproduce the source controller without inv(...)
% The source code recursively computes B*u(k). Backslash solves are used
% here instead of explicit matrix inverses, preserving the same formula
% while improving numerical conditioning.
nInputs = size(nominalInputs, 2);
correctedGeneralizedInput = zeros(2, nInputs);
correctedInputs = zeros(2, nInputs);
correctedGeneralizedInput(:,1) = B*nominalInputs(:,1);
correctedInputs(:,1) = nominalInputs(:,1);

for k = 2:nInputs
    remainingActual = scheduledTimes(k) - actualEvents(k);
    nominalInterval = scheduledTimes(k) - scheduledTimes(k-1);
    actualTransition = expm(A*remainingActual);
    nominalTransition = expm(A*nominalInterval);
    previousMismatch = correctedGeneralizedInput(:,k-1) - B*nominalInputs(:,k);
    correction = A*((actualTransition-eye(2)) \ ...
        ((actualTransition-nominalTransition)*(A \ previousMismatch)));
    correctedGeneralizedInput(:,k) = B*nominalInputs(:,k) + correction;
    correctedInputs(:,k) = B \ correctedGeneralizedInput(:,k);
end

%% Exact ZOH trajectories
nominalEventStates = eventStates(nominalEvents, nominalInputs, A, B, x0);
actualEventStates = eventStates(actualEvents, correctedInputs, A, B, x0);
eventGrid = [nominalEvents(nominalEvents <= horizon), ...
             actualEvents(actualEvents <= horizon)];
time = unique([0:options.SampleTime:horizon, eventGrid]);
nominalState = evaluateZoh(time, nominalEvents, nominalInputs, ...
    A, B, nominalEventStates);
actualState = evaluateZoh(time, actualEvents, correctedInputs, ...
    A, B, actualEventStates);
nominalInput = evaluateInput(time, nominalEvents, nominalInputs);
actualInput = evaluateInput(time, actualEvents, correctedInputs);

% Report how the perturbed and nominal states compare at scheduled updates.
checkTimes = scheduledTimes(scheduledTimes <= horizon);
nominalAtSchedule = evaluateZoh(checkTimes, nominalEvents, nominalInputs, ...
    A, B, nominalEventStates);
actualAtSchedule = evaluateZoh(checkTimes, actualEvents, correctedInputs, ...
    A, B, actualEventStates);
scheduleError = actualAtSchedule - nominalAtSchedule;

updateTable = table((1:numel(checkTimes)).', checkTimes.', ...
    deltas(1:numel(checkTimes)).', actualEvents(2:numel(checkTimes)+1).', ...
    'VariableNames', {'k','ScheduledTime','TimingPerturbation','ActualTime'});

fprintf('\nTwo-tank timing-perturbation reproduction\n');
fprintf('Initial water levels: h1 = %.3f m, h2 = %.3f m\n', x0(1), x0(2));
fprintf('Target water levels:  h1 = %.3f m, h2 = %.3f m\n', setpoint(1), setpoint(2));
fprintf('Admissible timing bound: 0 <= Delta(k) < %.3f s\n', deltaBound);
fprintf('Maximum scheduled-time state mismatch: %.3e m\n\n', ...
    max(abs(scheduleError), [], 'all'));
disp(updateTable);

data = struct( ...
    'A', A, 'B', B, 'x0', x0, 'setpoint', setpoint, ...
    'horizon', horizon, 'time', time, ...
    'nominalEvents', nominalEvents, 'actualEvents', actualEvents, ...
    'scheduledTimes', scheduledTimes, 'deltas', deltas, 'deltaBound', deltaBound, ...
    'nominalInputs', nominalInputs, 'correctedInputs', correctedInputs, ...
    'nominalInput', nominalInput, 'actualInput', actualInput, ...
    'nominalState', nominalState, 'actualState', actualState, ...
    'nominalEventStates', nominalEventStates, ...
    'actualEventStates', actualEventStates);

[fig, graphics] = createFigure(data, options.ShowFigure);

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
videoFrameSize = [];
if options.MakeVideo
    writer = VideoWriter(videoPath, 'MPEG-4');
    writer.FrameRate = options.FrameRate;
    writer.Quality = 100;
    open(writer);
end

wallClock = tic;
try
    for frameIndex = 1:numel(animationTimes)
        updateFigure(graphics, data, animationTimes(frameIndex));
        drawnow;
        if options.MakeVideo
            videoFrame = getframe(fig);
            if isempty(videoFrameSize)
                videoFrameSize = size(videoFrame.cdata, [1,2]);
            else
                videoFrame = normalizeVideoFrame(videoFrame, videoFrameSize);
            end
            writeVideo(writer, videoFrame);
        elseif options.PlayAnimation && frameIndex < numel(animationTimes)
            targetElapsed = frameIndex/options.FrameRate;
            pause(max(0, targetElapsed-toc(wallClock)));
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
    'A', A, 'B', B, 'x0', x0, 'setpoint', setpoint, ...
    'time', time, 'nominalState', nominalState, 'actualState', actualState, ...
    'nominalInput', nominalInput, 'actualInput', actualInput, ...
    'nominalEvents', nominalEvents, 'actualEvents', actualEvents, ...
    'deltaBound', deltaBound, ...
    'nominalInputs', nominalInputs, 'correctedInputs', correctedInputs, ...
    'scheduleError', scheduleError, 'updateTable', updateTable, ...
    'figureFile', figurePath, 'videoFile', videoPath);
end


function [fig, g] = createFigure(data, showFigure)
paper = [0.975, 0.98, 0.985];
ink = [0.08, 0.11, 0.15];
blue = [0.05, 0.42, 0.78];
teal = [0.02, 0.58, 0.54];
muted = [0.46, 0.51, 0.57];
gridColor = [0.82, 0.85, 0.88];

visibility = "off";
if showFigure
    visibility = "on";
end
fig = figure('Color', paper, 'Position', [60, 60, 1400, 780], ...
    'Name', 'Two-Tank Timing-Perturbation Visualization', ...
    'NumberTitle', 'off', 'Visible', visibility, 'Scrollable', 'off');

annotation(fig, 'textbox', [0.035, 0.925, 0.93, 0.06], ...
    'String', 'Two-Tank Water-Level Control under Timing Perturbations', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold', 'FontSize', 18, 'Color', ink);
statusText = annotation(fig, 'textbox', [0.06, 0.875, 0.88, 0.045], ...
    'String', 'Initialization', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', muted);

%% Animated physical schematic
axTank = axes(fig, 'Position', [0.035, 0.10, 0.43, 0.74]);
hold(axTank, 'on');
axis(axTank, [0, 10, 0, 7]);
axis(axTank, 'equal');
axis(axTank, 'off');
title(axTank, 'Physical water levels', 'FontWeight', 'normal', 'Color', ink);

tank1 = [1.0, 4.1, 1.0, 5.5];
tank2 = [5.8, 8.9, 1.0, 5.5];
drawTank(axTank, tank1, ink);
drawTank(axTank, tank2, ink);

maxLevel = 0.5;
fill1 = waterPatch(axTank, tank1, data.x0(1), maxLevel, blue);
fill2 = waterPatch(axTank, tank2, data.x0(2), maxLevel, teal);
surface1 = plot(axTank, tank1(1:2), levelY(tank1, data.x0(1), maxLevel)*[1,1], ...
    '-', 'Color', blue, 'LineWidth', 2.2);
surface2 = plot(axTank, tank2(1:2), levelY(tank2, data.x0(2), maxLevel)*[1,1], ...
    '-', 'Color', teal, 'LineWidth', 2.2);
nominalMarker1 = plot(axTank, tank1(1:2), ...
    levelY(tank1, data.x0(1), maxLevel)*[1,1], '--', 'Color', ink, 'LineWidth', 1.4);
nominalMarker2 = plot(axTank, tank2(1:2), ...
    levelY(tank2, data.x0(2), maxLevel)*[1,1], '--', 'Color', ink, 'LineWidth', 1.4);

% Setpoint ticks.
plot(axTank, [tank1(1)-0.15, tank1(1)+0.35], ...
    levelY(tank1, data.setpoint(1), maxLevel)*[1,1], ':', 'Color', muted, 'LineWidth', 1.5);
plot(axTank, [tank2(1)-0.15, tank2(1)+0.35], ...
    levelY(tank2, data.setpoint(2), maxLevel)*[1,1], ':', 'Color', muted, 'LineWidth', 1.5);

text(axTank, mean(tank1(1:2)), 0.55, 'Tank 1', 'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold', 'Color', ink);
text(axTank, mean(tank2(1:2)), 0.55, 'Tank 2', 'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold', 'Color', ink);
levelText1 = text(axTank, mean(tank1(1:2)), 3.2, 'h_1 = 0.100 m', ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', ink);
levelText2 = text(axTank, mean(tank2(1:2)), 3.2, 'h_2 = 0.100 m', ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', ink);
text(axTank, 0.25, levelY(tank1, data.setpoint(1), maxLevel), '0.40 m target', ...
    'Color', muted, 'FontSize', 9, 'VerticalAlignment', 'middle');
text(axTank, 9.0, levelY(tank2, data.setpoint(2), maxLevel), '0.23 m target', ...
    'Color', muted, 'FontSize', 9, 'VerticalAlignment', 'middle');

% Inlet, connection, and outlet piping.
plot(axTank, [2.55, 2.55], [6.55, 5.65], '-', 'Color', muted, 'LineWidth', 2.0);
plot(axTank, [7.35, 7.35], [6.55, 5.65], '-', 'Color', muted, 'LineWidth', 2.0);
inletHead1 = plot(axTank, 2.55, 5.65, 'v', 'MarkerSize', 9, ...
    'MarkerFaceColor', blue, 'MarkerEdgeColor', blue);
inletHead2 = plot(axTank, 7.35, 5.65, 'v', 'MarkerSize', 9, ...
    'MarkerFaceColor', teal, 'MarkerEdgeColor', teal);
inletText1 = text(axTank, 2.75, 6.2, 'q_1', 'Color', blue, 'FontSize', 11);
inletText2 = text(axTank, 7.55, 6.2, 'q_2', 'Color', teal, 'FontSize', 11);

plot(axTank, [tank1(2), tank2(1)], [1.45, 1.45], '-', 'Color', ink, 'LineWidth', 2.5);
connectionArrow = plot(axTank, [4.5, 5.35], [1.45, 1.45], '-', ...
    'Color', blue, 'LineWidth', 2.2, 'Marker', '>', 'MarkerIndices', 2, ...
    'MarkerFaceColor', blue);
connectionText = text(axTank, 4.95, 1.75, 'q_{12}', ...
    'HorizontalAlignment', 'center', 'Color', ink);

plot(axTank, [tank2(2), 9.65], [1.25, 1.25], '-', 'Color', ink, 'LineWidth', 2.5);
outletArrow = plot(axTank, [9.15, 9.65], [1.25, 1.25], '-', ...
    'Color', teal, 'LineWidth', 2.2, 'Marker', '>', 'MarkerIndices', 2, ...
    'MarkerFaceColor', teal);
text(axTank, 9.35, 1.55, 'q_{out}', 'HorizontalAlignment', 'center', 'Color', ink);

legendLines = [plot(axTank, nan, nan, '-', 'Color', blue, 'LineWidth', 5), ...
               plot(axTank, nan, nan, '--', 'Color', ink, 'LineWidth', 1.4)];
legend(axTank, legendLines, {'Perturbed/actual level', 'Nominal level marker'}, ...
    'Location', 'south', 'Box', 'off');

%% Water-level trajectories
axState = axes(fig, 'Position', [0.54, 0.53, 0.42, 0.31]);
hold(axState, 'on');
yline(axState, data.setpoint(1), ':', 'h_1 target', 'Color', blue, ...
    'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
yline(axState, data.setpoint(2), ':', 'h_2 target', 'Color', teal, ...
    'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
plot(axState, data.time, data.nominalState(1,:), '-', 'Color', ...
    0.58*blue + 0.42*[1,1,1], 'LineWidth', 1.6, 'DisplayName', 'h_1^n');
plot(axState, data.time, data.nominalState(2,:), '-', 'Color', ...
    0.58*teal + 0.42*[1,1,1], 'LineWidth', 1.6, 'DisplayName', 'h_2^n');
actualState1 = plot(axState, 0, data.x0(1), '--', 'Color', blue, ...
    'LineWidth', 2.2, 'DisplayName', 'h_1');
actualState2 = plot(axState, 0, data.x0(2), '--', 'Color', teal, ...
    'LineWidth', 2.2, 'DisplayName', 'h_2');
stateDot1 = plot(axState, 0, data.x0(1), 'o', 'MarkerFaceColor', blue, ...
    'MarkerEdgeColor', paper, 'MarkerSize', 7, 'HandleVisibility', 'off');
stateDot2 = plot(axState, 0, data.x0(2), 'o', 'MarkerFaceColor', teal, ...
    'MarkerEdgeColor', paper, 'MarkerSize', 7, 'HandleVisibility', 'off');
stateCursor = plot(axState, [0,0], [0,0.5], ':', 'Color', muted, ...
    'LineWidth', 1.0, 'HandleVisibility', 'off');
addScheduleLines(axState, data.scheduledTimes, data.horizon, gridColor);
xlim(axState, [0, data.horizon]);
ylim(axState, [0.075, 0.47]);
ylabel(axState, 'Water level (m)');
title(axState, 'Nominal and timing-perturbed water levels', 'FontWeight', 'normal');
legend(axState, 'Location', 'northeast', 'NumColumns', 2, 'Box', 'off');
styleAxes(axState, paper, gridColor);

%% Pump-flow trajectories
axInput = axes(fig, 'Position', [0.54, 0.105, 0.42, 0.31]);
hold(axInput, 'on');
[tn, un] = stepCoordinates(data.nominalEvents, 1e5*data.nominalInputs, data.horizon);
[ta, ua] = stepCoordinates(data.actualEvents, 1e5*data.correctedInputs, data.horizon);
plot(axInput, tn, un(1,:), '-', 'Color', 0.58*blue + 0.42*[1,1,1], ...
    'LineWidth', 1.5, 'DisplayName', 'q_1^n');
plot(axInput, tn, un(2,:), '-', 'Color', 0.58*teal + 0.42*[1,1,1], ...
    'LineWidth', 1.5, 'DisplayName', 'q_2^n');
actualInput1 = plot(axInput, ta(1), ua(1,1), '--', 'Color', blue, ...
    'LineWidth', 2.0, 'DisplayName', 'q_1');
actualInput2 = plot(axInput, ta(1), ua(2,1), '--', 'Color', teal, ...
    'LineWidth', 2.0, 'DisplayName', 'q_2');
inputCursor = plot(axInput, [0,0], [-4,14], ':', 'Color', muted, ...
    'LineWidth', 1.0, 'HandleVisibility', 'off');
addScheduleLines(axInput, data.scheduledTimes, data.horizon, gridColor);
xlim(axInput, [0, data.horizon]);
inputLimits = paddedLimits(1e5*[data.nominalInputs(:); data.correctedInputs(:)], 0.12);
ylim(axInput, inputLimits);
inputCursor.YData = inputLimits;
xlabel(axInput, 'Time (s)');
ylabel(axInput, 'Inlet flow (10^{-5} m^3/s)');
title(axInput, 'Piecewise-constant pump commands', 'FontWeight', 'normal');
legend(axInput, 'Location', 'northeast', 'NumColumns', 2, 'Box', 'off');
styleAxes(axInput, paper, gridColor);

g = struct( ...
    'statusText', statusText, 'tankAxes', axTank, ...
    'tank1', tank1, 'tank2', tank2, 'maxLevel', maxLevel, ...
    'fill1', fill1, 'fill2', fill2, 'surface1', surface1, 'surface2', surface2, ...
    'nominalMarker1', nominalMarker1, 'nominalMarker2', nominalMarker2, ...
    'levelText1', levelText1, 'levelText2', levelText2, ...
    'inletHead1', inletHead1, 'inletHead2', inletHead2, ...
    'inletText1', inletText1, 'inletText2', inletText2, ...
    'connectionArrow', connectionArrow, 'connectionText', connectionText, ...
    'outletArrow', outletArrow, ...
    'actualState1', actualState1, 'actualState2', actualState2, ...
    'stateDot1', stateDot1, 'stateDot2', stateDot2, 'stateCursor', stateCursor, ...
    'actualInput1', actualInput1, 'actualInput2', actualInput2, ...
    'inputCursor', inputCursor, 'actualStepTime', ta, 'actualStepInput', ua);
end


function updateFigure(g, data, currentTime)
actualNow = evaluateZoh(currentTime, data.actualEvents, data.correctedInputs, ...
    data.A, data.B, data.actualEventStates);
nominalNow = evaluateZoh(currentTime, data.nominalEvents, data.nominalInputs, ...
    data.A, data.B, data.nominalEventStates);
inputNow = evaluateInput(currentTime, data.actualEvents, data.correctedInputs);

% Update physical water levels.
setWaterPatch(g.fill1, g.tank1, actualNow(1), g.maxLevel);
setWaterPatch(g.fill2, g.tank2, actualNow(2), g.maxLevel);
y1 = levelY(g.tank1, actualNow(1), g.maxLevel);
y2 = levelY(g.tank2, actualNow(2), g.maxLevel);
yn1 = levelY(g.tank1, nominalNow(1), g.maxLevel);
yn2 = levelY(g.tank2, nominalNow(2), g.maxLevel);
g.surface1.YData = [y1,y1];
g.surface2.YData = [y2,y2];
g.nominalMarker1.YData = [yn1,yn1];
g.nominalMarker2.YData = [yn2,yn2];
g.levelText1.String = sprintf('h_1 = %.3f m', actualNow(1));
g.levelText2.String = sprintf('h_2 = %.3f m', actualNow(2));
g.inletText1.String = sprintf('q_1 = %+.2f x10^{-5}', 1e5*inputNow(1));
g.inletText2.String = sprintf('q_2 = %+.2f x10^{-5}', 1e5*inputNow(2));
updateInletMarker(g.inletHead1, inputNow(1), 2.55);
updateInletMarker(g.inletHead2, inputNow(2), 7.35);

flow12 = data.A(2,1)*(actualNow(1)-actualNow(2));
if flow12 >= 0
    g.connectionArrow.XData = [4.5, 5.35];
    g.connectionArrow.Marker = '>';
else
    g.connectionArrow.XData = [5.35, 4.5];
    g.connectionArrow.Marker = '<';
end
g.connectionText.String = sprintf('q_{12}  %+.4f', flow12);
g.outletArrow.LineWidth = 1.4 + 3.0*min(1, max(0, actualNow(2)/g.maxLevel));

% Update water-level traces.
mask = data.time <= currentTime + 1e-12;
traceTime = data.time(mask);
traceState = data.actualState(:,mask);
if traceTime(end) < currentTime
    traceTime(end+1) = currentTime;
    traceState(:,end+1) = actualNow;
end
g.actualState1.XData = traceTime;
g.actualState1.YData = traceState(1,:);
g.actualState2.XData = traceTime;
g.actualState2.YData = traceState(2,:);
g.stateDot1.XData = currentTime;
g.stateDot1.YData = actualNow(1);
g.stateDot2.XData = currentTime;
g.stateDot2.YData = actualNow(2);
g.stateCursor.XData = [currentTime,currentTime];
g.inputCursor.XData = [currentTime,currentTime];

% Update piecewise-constant input traces.
[stepTime, stepInput] = clipStep(g.actualStepTime, g.actualStepInput, currentTime);
g.actualInput1.XData = stepTime;
g.actualInput1.YData = stepInput(1,:);
g.actualInput2.XData = stepTime;
g.actualInput2.YData = stepInput(2,:);

revealed = find(data.actualEvents(2:end) <= currentTime + 1e-12, 1, 'last');
if isempty(revealed)
    detail = 'nominal plan active; no perturbed update has occurred';
else
    detail = sprintf('latest update k=%d: scheduled %.1f s, actual %.1f s, Delta=%+.1f s', ...
        revealed, data.scheduledTimes(revealed), data.actualEvents(revealed+1), ...
        data.deltas(revealed));
end
g.statusText.String = sprintf(['t = %.1f s  |  %s  |  ' ...
    'h = [%.3f, %.3f] m'], currentTime, detail, actualNow(1), actualNow(2));
end


function updateInletMarker(marker, flow, xPosition)
marker.XData = xPosition;
if flow >= 0
    marker.YData = 5.65;
    marker.Marker = 'v';
else
    marker.YData = 6.55;
    marker.Marker = '^';
end
marker.MarkerSize = 7 + 5*min(1, abs(flow)/1.2e-4);
end


function drawTank(ax, tank, color)
plot(ax, [tank(1),tank(1)], [tank(3),tank(4)], '-', 'Color', color, 'LineWidth', 3);
plot(ax, [tank(2),tank(2)], [tank(3),tank(4)], '-', 'Color', color, 'LineWidth', 3);
plot(ax, tank(1:2), [tank(3),tank(3)], '-', 'Color', color, 'LineWidth', 3);
end


function patchHandle = waterPatch(ax, tank, level, maxLevel, color)
surface = levelY(tank, level, maxLevel);
patchHandle = patch(ax, [tank(1),tank(2),tank(2),tank(1)], ...
    [tank(3),tank(3),surface,surface], color, ...
    'FaceAlpha', 0.48, 'EdgeColor', 'none');
end


function setWaterPatch(patchHandle, tank, level, maxLevel)
surface = levelY(tank, level, maxLevel);
patchHandle.XData = [tank(1),tank(2),tank(2),tank(1)];
patchHandle.YData = [tank(3),tank(3),surface,surface];
end


function y = levelY(tank, level, maxLevel)
usableHeight = tank(4)-tank(3)-0.18;
y = tank(3) + usableHeight*min(1, max(0, level/maxLevel));
end


function addScheduleLines(ax, times, horizon, color)
for time = times(times <= horizon)
    xline(ax, time, '-', 'Color', color, 'LineWidth', 0.7, ...
        'HandleVisibility', 'off');
end
end


function states = eventStates(events, inputs, A, B, x0)
states = zeros(numel(x0), numel(events));
states(:,1) = x0;
for k = 1:size(inputs,2)
    states(:,k+1) = zohStep(states(:,k), inputs(:,k), ...
        events(k+1)-events(k), A, B);
end
end


function nextState = zohStep(state, input, step, A, B)
n = size(A,1);
augmentedTransition = expm([A, B; zeros(size(B,2), n+size(B,2))]*step);
nextState = augmentedTransition(1:n,1:n)*state + ...
    augmentedTransition(1:n,n+1:end)*input;
end


function state = evaluateZoh(queryTimes, events, inputs, A, B, statesAtEvents)
state = zeros(size(A,1), numel(queryTimes));
for q = 1:numel(queryTimes)
    segment = find(events(1:end-1) <= queryTimes(q)+1e-12, 1, 'last');
    segment = min(segment, size(inputs,2));
    state(:,q) = zohStep(statesAtEvents(:,segment), inputs(:,segment), ...
        queryTimes(q)-events(segment), A, B);
end
end


function input = evaluateInput(queryTimes, events, inputs)
input = zeros(size(inputs,1), numel(queryTimes));
for q = 1:numel(queryTimes)
    segment = find(events(1:end-1) <= queryTimes(q)+1e-12, 1, 'last');
    segment = min(segment, size(inputs,2));
    input(:,q) = inputs(:,segment);
end
end


function [stepTime, stepInput] = stepCoordinates(events, inputs, horizon)
nInputs = size(inputs,2);
stepTime = zeros(1,2*nInputs);
stepInput = zeros(size(inputs,1),2*nInputs);
for k = 1:nInputs
    indices = 2*k-1:2*k;
    stepTime(indices) = [events(k), min(events(k+1),horizon)];
    stepInput(:,indices) = inputs(:,[k,k]);
end
keep = stepTime <= horizon + 1e-12;
stepTime = stepTime(keep);
stepInput = stepInput(:,keep);
end


function [clippedTime, clippedInput] = clipStep(stepTime, stepInput, limit)
mask = stepTime <= limit + 1e-12;
clippedTime = stepTime(mask);
clippedInput = stepInput(:,mask);
if clippedTime(end) < limit
    segment = find(stepTime <= limit, 1, 'last');
    clippedTime(end+1) = limit;
    clippedInput(:,end+1) = stepInput(:,segment);
end
end


function limits = paddedLimits(values, fraction)
lower = min(values);
upper = max(values);
span = upper-lower;
if span < 1e-9
    span = max(1,abs(upper));
end
limits = [lower-fraction*span, upper+fraction*span];
end


function styleAxes(ax, background, gridColor)
ax.Color = background;
ax.Box = 'off';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridColor = gridColor;
ax.GridAlpha = 0.55;
ax.FontName = 'Times New Roman';
ax.FontSize = 10.5;
end


function frame = normalizeVideoFrame(frame, targetSize)
% Windows display scaling can change getframe dimensions when figure
% scrollability changes. Resample only when needed so VideoWriter receives
% a constant frame size for the entire animation.
currentSize = size(frame.cdata, [1,2]);
if isequal(currentSize, targetSize)
    return;
end
rowIndex = round(linspace(1, currentSize(1), targetSize(1)));
columnIndex = round(linspace(1, currentSize(2), targetSize(2)));
frame.cdata = frame.cdata(rowIndex, columnIndex, :);
frame.colormap = [];
end


function outputPath = resolveOutputPath(requestedPath, baseFolder)
[folder, name, extension] = fileparts(char(requestedPath));
if isempty(folder)
    outputPath = fullfile(baseFolder,[name,extension]);
else
    outputPath = char(requestedPath);
end
end
