(function () {
  'use strict';

  var root = document.querySelector('[data-blimp-demo]');
  if (!root) return;

  var canvas = root.querySelector('[data-canvas]');
  var context = canvas.getContext('2d');
  var toggleButton = root.querySelector('[data-action="toggle"]');
  var restartButton = root.querySelector('[data-action="restart"]');
  var resampleButton = root.querySelector('[data-action="resample"]');
  var blimpSelect = root.querySelector('[data-control="blimp"]');
  var speedSelect = root.querySelector('[data-control="speed"]');
  var sampleSelect = root.querySelector('[data-control="samples"]');
  var layerInputs = root.querySelectorAll('[data-layer]');
  var statusElement = root.querySelector('[data-status]');
  var scheduleSvg = root.querySelector('[data-schedule-chart]');
  var scheduleCurrentElement = root.querySelector('[data-schedule-current]');
  var scheduleMarker = null;
  var scheduleMarkerLabel = null;
  var scheduleGeometry = null;
  var scheduleBlocks = [];

  var WORLD_WIDTH = 9;
  var WORLD_HEIGHT = 5;
  var GRID_LENGTH_CM = 100;
  var BLIMP_SPEED = 10;
  var BASE_EXECUTION = 2;
  var RESOURCE_EFFICIENCY = [1.2, 1.0];
  var WIND_MEAN = 5;
  var WIND_SIGMA = 0.5;
  var ANGLE_SIGMA = 0.1;
  var BASE_PLAYBACK_RATE = 5;
  var TWO_PI = Math.PI * 2;
  var SVG_NS = 'http://www.w3.org/2000/svg';

  var routes = [
    {
      id: 1,
      colorToken: '--blimp-1',
      points: [[5, 4.5], [6, 4.5], [7, 4.5], [8, 4.5], [8.5, 4], [8, 3.5], [7, 3.5], [6, 3.5], [5, 3.5], [4.5, 4], [5, 4.5]],
      betas: [Math.PI / 4, Math.PI / 4, Math.PI / 4, 0, -Math.PI / 2, -3 * Math.PI / 4, -3 * Math.PI / 4, -3 * Math.PI / 4, -Math.PI, Math.PI / 2]
    },
    {
      id: 2,
      colorToken: '--blimp-2',
      points: [[1, 0.5], [2, 0.5], [3, 0.5], [4, 0.5], [4.5, 1], [4, 1.5], [3, 1.5], [2, 1.5], [1, 1.5], [0.5, 1], [1, 0.5]],
      betas: [Math.PI / 4, Math.PI / 4, Math.PI / 4, Math.PI / 2, -Math.PI, -3 * Math.PI / 4, -3 * Math.PI / 4, -3 * Math.PI / 4, -Math.PI / 2, 0]
    },
    {
      id: 3,
      colorToken: '--blimp-3',
      points: [[1, 2.5], [2, 2.5], [3, 2.5], [3.5, 3], [3, 3.5], [2, 3.5], [1, 3.5], [0.5, 4], [1, 4.5], [2, 4.5], [3, 4.5]],
      betas: [Math.PI / 4, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4, -Math.PI, -3 * Math.PI / 4, -3 * Math.PI / 4, -Math.PI / 2, -Math.PI / 4, 0]
    },
    {
      id: 4,
      colorToken: '--blimp-4',
      points: [[6, 0.5], [7, 0.5], [8, 0.5], [8.5, 1], [8.5, 2], [8, 2.5], [7, 2.5], [6, 2.5], [5.5, 2], [5.5, 1], [6, 0.5]],
      betas: [Math.PI / 4, Math.PI / 4, Math.PI / 2, -Math.PI, -3 * Math.PI / 4, -3 * Math.PI / 4, -Math.PI, Math.PI / 2, Math.PI / 4, Math.PI / 4]
    }
  ];

  var state = {
    time: 0,
    playing: true,
    speed: 1,
    selectedBlimp: 0,
    sampleCount: 24,
    seed: 20260824,
    showPredicted: true,
    showMonteCarlo: true,
    showWind: true,
    showProcessorLinks: true,
    lastFrame: performance.now(),
    schedule: null,
    scenario: null,
    ensembles: null,
    colors: null
  };

  function mulberry32(seed) {
    var value = seed >>> 0;
    return function () {
      value += 0x6D2B79F5;
      var result = value;
      result = Math.imul(result ^ (result >>> 15), result | 1);
      result ^= result + Math.imul(result ^ (result >>> 7), result | 61);
      return ((result ^ (result >>> 14)) >>> 0) / 4294967296;
    };
  }

  function createNormal(random) {
    var spare = null;
    return function () {
      if (spare !== null) {
        var saved = spare;
        spare = null;
        return saved;
      }
      var u = 0;
      var v = 0;
      while (u === 0) u = random();
      while (v === 0) v = random();
      var magnitude = Math.sqrt(-2 * Math.log(u));
      var angle = TWO_PI * v;
      spare = magnitude * Math.sin(angle);
      return magnitude * Math.cos(angle);
    };
  }

  function clamp(value, minimum, maximum) {
    return Math.max(minimum, Math.min(maximum, value));
  }

  function normalizeAngle(angle) {
    var normalized = angle;
    while (normalized > Math.PI) normalized -= TWO_PI;
    while (normalized <= -Math.PI) normalized += TWO_PI;
    return normalized;
  }

  function distance(pointA, pointB) {
    return Math.hypot(pointB[0] - pointA[0], pointB[1] - pointA[1]);
  }

  function segmentAngle(pointA, pointB) {
    return Math.atan2(-(pointB[1] - pointA[1]), pointB[0] - pointA[0]);
  }

  function safeRatio(numerator, denominator) {
    if (Math.abs(denominator) < 1e-7) return Number.POSITIVE_INFINITY;
    return numerator / denominator;
  }

  function taskPeriod(delay, execution, windMagnitude, beta, segmentLength) {
    var distanceCm = segmentLength * GRID_LENGTH_CM;
    var delta = delay + execution;
    var kx = windMagnitude * Math.cos(beta) / BLIMP_SPEED;
    var ky = windMagnitude * Math.sin(beta) / BLIMP_SPEED;
    var k = 1 - kx * kx - ky * ky;
    var radicand =
      Math.pow(1 + kx, 2) * Math.pow(BLIMP_SPEED, 2) * Math.pow(delta, 2) -
      2 * distanceCm * (1 + kx - ky * ky) * BLIMP_SPEED * delta +
      (1 - ky * ky) * Math.pow(distanceCm, 2);
    var numerator =
      -(kx * distanceCm - (1 + kx) * BLIMP_SPEED * delta) +
      Math.sqrt(Math.max(0, radicand));
    var result = numerator / (Math.max(0.05, k) * BLIMP_SPEED);
    var lowerBound = delta + 0.35;
    var nominalUpper = distanceCm / BLIMP_SPEED * 1.9;
    if (!Number.isFinite(result)) result = lowerBound;
    return clamp(result, lowerBound, Math.max(lowerBound, nominalUpper));
  }

  function robustDelayBound(route, taskIndex, execution) {
    var start = route.points[taskIndex];
    var end = route.points[taskIndex + 1];
    var dx = end[0] - start[0];
    var dyUp = -(end[1] - start[1]);
    var diagonal = Math.abs(Math.abs(dx) - Math.abs(dyUp)) < 1e-6 && Math.abs(dx) > 0.1;
    var pathAngle = Math.atan2(dyUp, dx);
    var betaH = normalizeAngle(pathAngle + route.betas[taskIndex]);
    var sinAngle = Math.sin(betaH);
    var cosAngle = Math.cos(betaH);
    var rootTwo = Math.sqrt(2);
    var d1 = safeRatio(rootTwo * GRID_LENGTH_CM, BLIMP_SPEED * (rootTwo + cosAngle));
    var d2 = safeRatio(GRID_LENGTH_CM, rootTwo * BLIMP_SPEED * sinAngle);
    var d3 = safeRatio(rootTwo * GRID_LENGTH_CM, BLIMP_SPEED * (1 + cosAngle));
    var d4 = safeRatio(GRID_LENGTH_CM, rootTwo * BLIMP_SPEED * (1 + sinAngle));
    var d5 = safeRatio(GRID_LENGTH_CM, rootTwo * BLIMP_SPEED * (1 - sinAngle));
    var deltaMaximum;

    if (!diagonal) {
      if (betaH >= 0 && betaH <= Math.PI / 2) deltaMaximum = Math.min(d1, d2);
      else if (betaH > Math.PI / 2) deltaMaximum = Math.min(GRID_LENGTH_CM / BLIMP_SPEED, d2);
      else if (betaH >= -Math.PI / 2) deltaMaximum = Math.min(d1, -d2);
      else deltaMaximum = Math.min(GRID_LENGTH_CM / BLIMP_SPEED, -d2);
    } else {
      var typeTwo = Math.sin(2 * pathAngle) >= 0;
      if (typeTwo) {
        if (betaH >= 0 && betaH <= Math.PI / 2) deltaMaximum = Math.min(d3, d4);
        else if (betaH > Math.PI / 2) deltaMaximum = Math.min(rootTwo * GRID_LENGTH_CM / BLIMP_SPEED, d4);
        else deltaMaximum = GRID_LENGTH_CM / (rootTwo * BLIMP_SPEED);
      } else {
        if (betaH >= 0) deltaMaximum = GRID_LENGTH_CM / (rootTwo * BLIMP_SPEED);
        else if (betaH >= -Math.PI / 2) deltaMaximum = Math.min(d3, d5);
        else deltaMaximum = Math.min(rootTwo * GRID_LENGTH_CM / BLIMP_SPEED, d5);
      }
    }

    if (!Number.isFinite(deltaMaximum)) deltaMaximum = GRID_LENGTH_CM / BLIMP_SPEED;
    return Math.max(0.25, deltaMaximum - execution);
  }

  function buildSchedule() {
    var tasks = routes.map(function () { return new Array(10); });
    var resourceRecords = [[], []];
    var resources = [{ current: null }, { current: null }];
    var waiting = [];
    var events = routes.map(function (_, blimpIndex) {
      return { time: 0, type: 'arrival', blimpIndex: blimpIndex, taskIndex: 0 };
    });
    var scheduledCount = 0;
    var guard = 0;

    function priority(job, now) {
      var route = routes[job.blimpIndex];
      var beta = route.betas[job.taskIndex];
      var windSensitivity = Math.abs(Math.sin(beta)) * WIND_MEAN;
      var wait = Math.max(0, now - job.arrival);
      return wait * 10 + windSensitivity - route.id * 0.001;
    }

    while ((events.length || waiting.length) && scheduledCount < 40 && guard < 1000) {
      guard += 1;
      events.sort(function (a, b) { return a.time - b.time; });
      if (!events.length) break;
      var currentTime = events[0].time;
      var currentEvents = [];
      while (events.length && Math.abs(events[0].time - currentTime) < 1e-7) {
        currentEvents.push(events.shift());
      }

      currentEvents.forEach(function (event) {
        if (event.type === 'complete') {
          resources[event.resourceIndex].current = null;
        }
      });
      currentEvents.forEach(function (event) {
        if (event.type === 'arrival') {
          waiting.push({
            blimpIndex: event.blimpIndex,
            taskIndex: event.taskIndex,
            arrival: event.time
          });
        }
      });

      var idleResources = resources
        .map(function (resource, index) { return resource.current === null ? index : -1; })
        .filter(function (index) { return index >= 0; })
        .sort(function (a, b) { return RESOURCE_EFFICIENCY[b] - RESOURCE_EFFICIENCY[a]; });

      waiting.sort(function (a, b) {
        var scoreDifference = priority(b, currentTime) - priority(a, currentTime);
        if (Math.abs(scoreDifference) > 1e-7) return scoreDifference;
        return routes[a.blimpIndex].id - routes[b.blimpIndex].id;
      });

      idleResources.forEach(function (resourceIndex) {
        if (!waiting.length) return;
        var job = waiting.shift();
        var route = routes[job.blimpIndex];
        var execution = BASE_EXECUTION / RESOURCE_EFFICIENCY[resourceIndex];
        var delay = Math.max(0, currentTime - job.arrival);
        var segmentLength = distance(route.points[job.taskIndex], route.points[job.taskIndex + 1]);
        var period = taskPeriod(delay, execution, WIND_MEAN, route.betas[job.taskIndex], segmentLength);
        var exit = job.arrival + period;
        var completion = currentTime + execution;
        var bound = robustDelayBound(route, job.taskIndex, execution);
        var record = {
          blimpIndex: job.blimpIndex,
          taskIndex: job.taskIndex,
          arrival: job.arrival,
          start: currentTime,
          gamma: completion,
          exit: exit,
          delay: delay,
          execution: execution,
          resourceIndex: resourceIndex,
          period: period,
          bound: bound
        };

        tasks[job.blimpIndex][job.taskIndex] = record;
        resources[resourceIndex].current = record;
        resourceRecords[resourceIndex].push(record);
        scheduledCount += 1;

        events.push({
          time: completion,
          type: 'complete',
          resourceIndex: resourceIndex
        });
        if (job.taskIndex < 9) {
          events.push({
            time: exit,
            type: 'arrival',
            blimpIndex: job.blimpIndex,
            taskIndex: job.taskIndex + 1
          });
        }
      });
    }

    var totalTime = Math.max.apply(null, tasks.map(function (blimpTasks) {
      return blimpTasks[blimpTasks.length - 1].exit;
    }));

    return {
      tasks: tasks,
      resourceRecords: resourceRecords,
      totalTime: totalTime
    };
  }

  function sampleWind(normal, meanAngle) {
    var maximumWind = BLIMP_SPEED / Math.sqrt(2) * 0.995;
    return {
      magnitude: clamp(WIND_MEAN + WIND_SIGMA * normal(), 0.2, maximumWind),
      beta: normalizeAngle(meanAngle + ANGLE_SIGMA * normal())
    };
  }

  function driftPoint(start, end, delta, wind) {
    var segmentLength = distance(start, end);
    var ux = (end[0] - start[0]) / segmentLength;
    var uy = (end[1] - start[1]) / segmentLength;
    var perpendicularX = -uy;
    var perpendicularY = ux;
    var forward = (BLIMP_SPEED + wind.magnitude * Math.cos(wind.beta)) * delta / GRID_LENGTH_CM;
    var lateral = wind.magnitude * Math.sin(wind.beta) * delta / GRID_LENGTH_CM;
    return [
      start[0] + ux * forward + perpendicularX * lateral,
      start[1] + uy * forward + perpendicularY * lateral
    ];
  }

  function createScenario(seed) {
    var random = mulberry32(seed);
    var normal = createNormal(random);
    return routes.map(function (route, blimpIndex) {
      return state.schedule.tasks[blimpIndex].map(function (task, taskIndex) {
        var start = route.points[taskIndex];
        var end = route.points[taskIndex + 1];
        var predictedWind = { magnitude: WIND_MEAN, beta: route.betas[taskIndex] };
        var actualWind = sampleWind(normal, route.betas[taskIndex]);
        var delta = task.delay + task.execution;
        return {
          start: start,
          end: end,
          predictedWind: predictedWind,
          actualWind: actualWind,
          predictedDrift: driftPoint(start, end, delta, predictedWind),
          actualDrift: driftPoint(start, end, delta, actualWind)
        };
      });
    });
  }

  function createEnsembles(seed, sampleCount) {
    var ensembles = [];
    for (var sampleIndex = 0; sampleIndex < sampleCount; sampleIndex += 1) {
      var random = mulberry32((seed + (sampleIndex + 1) * 104729) >>> 0);
      var normal = createNormal(random);
      var sample = routes.map(function (route, blimpIndex) {
        return state.schedule.tasks[blimpIndex].map(function (task, taskIndex) {
          var wind = sampleWind(normal, route.betas[taskIndex]);
          return driftPoint(
            route.points[taskIndex],
            route.points[taskIndex + 1],
            task.delay + task.execution,
            wind
          );
        });
      });
      ensembles.push(sample);
    }
    return ensembles;
  }

  function readColors() {
    var styles = getComputedStyle(document.documentElement);
    function color(name, fallback) {
      return styles.getPropertyValue(name).trim() || fallback;
    }
    return {
      background: color('--color-bg', '#ffffff'),
      text: color('--color-text', '#202522'),
      soft: color('--color-text-soft', '#515a54'),
      muted: color('--color-text-muted', '#737b76'),
      border: color('--color-border', '#e2e5e1'),
      borderStrong: color('--color-border-strong', '#cbd1cc'),
      accentSoft: color('--color-accent-soft', '#edf5f0'),
      wind: color('--blimp-wind', '#2878a9'),
      routeColors: routes.map(function (route) {
        return color(route.colorToken, '#0b6241');
      })
    };
  }

  function createSvgNode(name, attributes, text) {
    var node = document.createElementNS(SVG_NS, name);
    Object.keys(attributes || {}).forEach(function (key) {
      node.setAttribute(key, attributes[key]);
    });
    if (typeof text === 'string') node.textContent = text;
    return node;
  }

  function scheduleX(time) {
    if (!scheduleGeometry) return 0;
    return scheduleGeometry.left + clamp(time / scheduleGeometry.totalTime, 0, 1) * scheduleGeometry.width;
  }

  function renderScheduleChart() {
    if (!scheduleSvg) return;

    while (scheduleSvg.firstChild) scheduleSvg.removeChild(scheduleSvg.firstChild);
    scheduleBlocks = [];

    var width = 1240;
    var height = 220;
    var plotLeft = 150;
    var plotRight = 28;
    var plotWidth = width - plotLeft - plotRight;
    var rowHeight = 38;
    var rowYs = [60, 122];
    var gridTop = 48;
    var gridBottom = 174;
    var totalTime = state.schedule.totalTime;
    var tickCount = 8;

    scheduleGeometry = {
      left: plotLeft,
      width: plotWidth,
      totalTime: totalTime,
      top: gridTop,
      bottom: gridBottom
    };

    scheduleSvg.setAttribute('viewBox', '0 0 ' + width + ' ' + height);
    scheduleSvg.appendChild(createSvgNode('title', {}, 'Two-processor wind-estimation task schedule'));
    scheduleSvg.appendChild(createSvgNode('desc', {}, 'Colored blocks show the blimp tasks processed by Processor A and Processor B. Dashed segments show waiting time before processing, and a moving vertical line shows the current animation time.'));

    for (var tickIndex = 0; tickIndex <= tickCount; tickIndex += 1) {
      var tickTime = totalTime * tickIndex / tickCount;
      var tickX = plotLeft + plotWidth * tickIndex / tickCount;
      scheduleSvg.appendChild(createSvgNode('line', {
        x1: tickX,
        y1: gridTop,
        x2: tickX,
        y2: gridBottom,
        'class': 'blimp-schedule-grid'
      }));
      scheduleSvg.appendChild(createSvgNode('text', {
        x: tickX,
        y: 198,
        'text-anchor': tickIndex === 0 ? 'start' : (tickIndex === tickCount ? 'end' : 'middle'),
        'class': 'blimp-schedule-axis-label'
      }, tickTime.toFixed(0) + ' s'));
    }

    RESOURCE_EFFICIENCY.forEach(function (efficiency, resourceIndex) {
      var rowY = rowYs[resourceIndex];
      var processorName = resourceIndex === 0 ? 'Processor A' : 'Processor B';
      var processorMeta = resourceIndex === 0 ? 'fast · 1.2×' : 'standard · 1.0×';

      scheduleSvg.appendChild(createSvgNode('rect', {
        x: plotLeft,
        y: rowY,
        width: plotWidth,
        height: rowHeight,
        rx: 2,
        'class': 'blimp-schedule-row'
      }));
      scheduleSvg.appendChild(createSvgNode('text', {
        x: 18,
        y: rowY + 15,
        'class': 'blimp-schedule-row-label'
      }, processorName));
      scheduleSvg.appendChild(createSvgNode('text', {
        x: 18,
        y: rowY + 31,
        'class': 'blimp-schedule-row-meta'
      }, processorMeta));

      state.schedule.resourceRecords[resourceIndex].forEach(function (record) {
        var waitStartX = scheduleX(record.arrival);
        var blockStartX = scheduleX(record.start);
        var blockEndX = scheduleX(record.gamma);
        var blockWidth = Math.max(3, blockEndX - blockStartX);
        var route = routes[record.blimpIndex];
        var label = 'B' + route.id;
        var fullLabel = 'Blimp ' + route.id + ' · Cell ' + (record.taskIndex + 1) + ' · ' + processorName + ' · released ' + record.arrival.toFixed(2) + ' s · processed ' + record.start.toFixed(2) + '–' + record.gamma.toFixed(2) + ' s';

        if (record.start - record.arrival > 0.02) {
          scheduleSvg.appendChild(createSvgNode('line', {
            x1: waitStartX,
            y1: rowY + rowHeight / 2,
            x2: blockStartX,
            y2: rowY + rowHeight / 2,
            'class': 'blimp-schedule-wait'
          }));
        }

        var block = createSvgNode('rect', {
          x: blockStartX,
          y: rowY + 5,
          width: blockWidth,
          height: rowHeight - 10,
          rx: 2,
          fill: state.colors.routeColors[record.blimpIndex],
          'class': 'blimp-schedule-block',
          tabindex: '0',
          'aria-label': fullLabel
        });
        block.appendChild(createSvgNode('title', {}, fullLabel));
        scheduleSvg.appendChild(block);
        scheduleBlocks.push({ element: block, record: record });

        if (blockWidth >= 14) {
          scheduleSvg.appendChild(createSvgNode('text', {
            x: blockStartX + blockWidth / 2,
            y: rowY + 24,
            'text-anchor': 'middle',
            'class': 'blimp-schedule-block-label'
          }, label));
        }
      });
    });

    scheduleMarker = createSvgNode('line', {
      x1: plotLeft,
      y1: gridTop - 5,
      x2: plotLeft,
      y2: gridBottom,
      'class': 'blimp-schedule-marker'
    });
    scheduleMarkerLabel = createSvgNode('text', {
      x: plotLeft + 5,
      y: 38,
      'text-anchor': 'start',
      'class': 'blimp-schedule-time-label'
    }, 't = 0.0 s');
    scheduleSvg.appendChild(scheduleMarker);
    scheduleSvg.appendChild(scheduleMarkerLabel);
    updateScheduleMarker();
  }

  function updateScheduleMarker() {
    if (!scheduleMarker || !scheduleMarkerLabel || !scheduleGeometry) return;
    var markerX = scheduleX(state.time);
    var nearRightEdge = markerX > scheduleGeometry.left + scheduleGeometry.width * 0.9;
    scheduleMarker.setAttribute('x1', markerX);
    scheduleMarker.setAttribute('x2', markerX);
    scheduleMarkerLabel.setAttribute('x', markerX + (nearRightEdge ? -5 : 5));
    scheduleMarkerLabel.setAttribute('text-anchor', nearRightEdge ? 'end' : 'start');
    scheduleMarkerLabel.textContent = 't = ' + state.time.toFixed(1) + ' s';

    scheduleBlocks.forEach(function (item) {
      var active = state.time >= item.record.start && state.time < item.record.gamma;
      item.element.classList.toggle('is-active', active);
    });

    if (scheduleCurrentElement) {
      var assignments = RESOURCE_EFFICIENCY.map(function (_, resourceIndex) {
        var record = activeResourceRecord(resourceIndex, state.time);
        var processor = resourceIndex === 0 ? 'A' : 'B';
        return record
          ? processor + ': Blimp ' + routes[record.blimpIndex].id + ' · Cell ' + (record.taskIndex + 1)
          : processor + ': idle';
      });
      scheduleCurrentElement.textContent = 't = ' + state.time.toFixed(1) + ' s · ' + assignments.join(' · ');
    }
  }

  function layout() {
    var width = canvas.clientWidth || 960;
    var height = canvas.clientHeight || Math.round(width * 0.56);
    var dpr = Math.min(window.devicePixelRatio || 1, 2);
    if (canvas.width !== Math.round(width * dpr) || canvas.height !== Math.round(height * dpr)) {
      canvas.width = Math.round(width * dpr);
      canvas.height = Math.round(height * dpr);
    }
    context.setTransform(dpr, 0, 0, dpr, 0, 0);
    var padding = {
      left: width < 560 ? 24 : 42,
      right: width < 560 ? 18 : 30,
      top: width < 560 ? 34 : 42,
      bottom: width < 560 ? 28 : 38
    };
    var scale = Math.min(
      (width - padding.left - padding.right) / WORLD_WIDTH,
      (height - padding.top - padding.bottom) / WORLD_HEIGHT
    );
    var plotWidth = scale * WORLD_WIDTH;
    var plotHeight = scale * WORLD_HEIGHT;
    var originX = padding.left + (width - padding.left - padding.right - plotWidth) / 2;
    var originY = padding.top + (height - padding.top - padding.bottom - plotHeight) / 2;
    return {
      width: width,
      height: height,
      scale: scale,
      originX: originX,
      originY: originY,
      plotWidth: plotWidth,
      plotHeight: plotHeight
    };
  }

  function mapPoint(point, dimensions) {
    return [
      dimensions.originX + point[0] * dimensions.scale,
      dimensions.originY + point[1] * dimensions.scale
    ];
  }

  function currentTaskIndex(blimpIndex, time) {
    var tasks = state.schedule.tasks[blimpIndex];
    for (var index = 0; index < tasks.length; index += 1) {
      if (time < tasks[index].exit) return index;
    }
    return tasks.length - 1;
  }

  function pointOnTask(blimpIndex, taskIndex, time, usePredicted) {
    var task = state.schedule.tasks[blimpIndex][taskIndex];
    var geometry = state.scenario[blimpIndex][taskIndex];
    var drift = usePredicted ? geometry.predictedDrift : geometry.actualDrift;
    if (time <= task.arrival) {
      return { point: geometry.start, angle: Math.atan2(geometry.end[1] - geometry.start[1], geometry.end[0] - geometry.start[0]) };
    }
    if (time < task.gamma) {
      var firstFraction = clamp((time - task.arrival) / Math.max(0.001, task.gamma - task.arrival), 0, 1);
      return {
        point: [
          geometry.start[0] + (drift[0] - geometry.start[0]) * firstFraction,
          geometry.start[1] + (drift[1] - geometry.start[1]) * firstFraction
        ],
        angle: Math.atan2(drift[1] - geometry.start[1], drift[0] - geometry.start[0])
      };
    }
    var secondFraction = clamp((time - task.gamma) / Math.max(0.001, task.exit - task.gamma), 0, 1);
    return {
      point: [
        drift[0] + (geometry.end[0] - drift[0]) * secondFraction,
        drift[1] + (geometry.end[1] - drift[1]) * secondFraction
      ],
      angle: Math.atan2(geometry.end[1] - drift[1], geometry.end[0] - drift[0])
    };
  }

  function currentPosition(blimpIndex, time, usePredicted) {
    var taskIndex = currentTaskIndex(blimpIndex, time);
    var task = state.schedule.tasks[blimpIndex][taskIndex];
    if (time >= task.exit && taskIndex === 9) {
      var finalGeometry = state.scenario[blimpIndex][taskIndex];
      return {
        point: finalGeometry.end,
        angle: Math.atan2(finalGeometry.end[1] - finalGeometry.actualDrift[1], finalGeometry.end[0] - finalGeometry.actualDrift[0]),
        taskIndex: taskIndex
      };
    }
    var result = pointOnTask(blimpIndex, taskIndex, time, usePredicted);
    result.taskIndex = taskIndex;
    return result;
  }

  function drawGrid(dimensions) {
    context.fillStyle = state.colors.background;
    context.fillRect(0, 0, dimensions.width, dimensions.height);
    context.fillStyle = state.colors.accentSoft;
    context.globalAlpha = 0.28;
    context.fillRect(dimensions.originX, dimensions.originY, dimensions.plotWidth, dimensions.plotHeight);
    context.globalAlpha = 1;
    context.strokeStyle = state.colors.border;
    context.lineWidth = 1;
    for (var column = 0; column <= WORLD_WIDTH; column += 1) {
      var x = dimensions.originX + column * dimensions.scale;
      context.beginPath();
      context.moveTo(x, dimensions.originY);
      context.lineTo(x, dimensions.originY + dimensions.plotHeight);
      context.stroke();
    }
    for (var row = 0; row <= WORLD_HEIGHT; row += 1) {
      var y = dimensions.originY + row * dimensions.scale;
      context.beginPath();
      context.moveTo(dimensions.originX, y);
      context.lineTo(dimensions.originX + dimensions.plotWidth, y);
      context.stroke();
    }
    context.strokeStyle = state.colors.borderStrong;
    context.lineWidth = 1.3;
    context.strokeRect(dimensions.originX, dimensions.originY, dimensions.plotWidth, dimensions.plotHeight);
    context.fillStyle = state.colors.muted;
    context.font = '600 12px Lato, sans-serif';
    context.textAlign = 'left';
    context.fillText('9 m × 5 m inspection grid', dimensions.originX, Math.max(15, dimensions.originY - 14));
  }

  function drawNominalRoutes(dimensions) {
    routes.forEach(function (route, blimpIndex) {
      context.beginPath();
      route.points.forEach(function (point, index) {
        var mapped = mapPoint(point, dimensions);
        if (index === 0) context.moveTo(mapped[0], mapped[1]);
        else context.lineTo(mapped[0], mapped[1]);
      });
      context.strokeStyle = blimpIndex === state.selectedBlimp ? state.colors.routeColors[blimpIndex] : state.colors.borderStrong;
      context.globalAlpha = blimpIndex === state.selectedBlimp ? 0.5 : 0.7;
      context.lineWidth = blimpIndex === state.selectedBlimp ? 2 : 1.15;
      context.setLineDash([2, 5]);
      context.stroke();
      context.setLineDash([]);
      context.globalAlpha = 1;

      var start = mapPoint(route.points[0], dimensions);
      context.fillStyle = state.colors.routeColors[blimpIndex];
      context.font = '700 11px Lato, sans-serif';
      context.textAlign = 'center';
      context.fillText('B' + route.id, start[0], start[1] - 10);
    });
  }

  function drawMonteCarlo(dimensions) {
    if (!state.showMonteCarlo) return;
    var blimpIndex = state.selectedBlimp;
    var route = routes[blimpIndex];
    context.save();
    context.strokeStyle = state.colors.routeColors[blimpIndex];
    context.lineWidth = 0.8;
    context.globalAlpha = state.sampleCount > 24 ? 0.075 : 0.1;
    state.ensembles.forEach(function (sample) {
      context.beginPath();
      route.points.slice(0, -1).forEach(function (start, taskIndex) {
        var mappedStart = mapPoint(start, dimensions);
        var mappedDrift = mapPoint(sample[blimpIndex][taskIndex], dimensions);
        var mappedEnd = mapPoint(route.points[taskIndex + 1], dimensions);
        if (taskIndex === 0) context.moveTo(mappedStart[0], mappedStart[1]);
        else context.lineTo(mappedStart[0], mappedStart[1]);
        context.lineTo(mappedDrift[0], mappedDrift[1]);
        context.lineTo(mappedEnd[0], mappedEnd[1]);
      });
      context.stroke();
    });
    context.restore();
  }

  function drawPredictedRoute(dimensions) {
    if (!state.showPredicted) return;
    var blimpIndex = state.selectedBlimp;
    var route = routes[blimpIndex];
    context.beginPath();
    route.points.slice(0, -1).forEach(function (start, taskIndex) {
      var mappedStart = mapPoint(start, dimensions);
      var mappedDrift = mapPoint(state.scenario[blimpIndex][taskIndex].predictedDrift, dimensions);
      var mappedEnd = mapPoint(route.points[taskIndex + 1], dimensions);
      if (taskIndex === 0) context.moveTo(mappedStart[0], mappedStart[1]);
      else context.lineTo(mappedStart[0], mappedStart[1]);
      context.lineTo(mappedDrift[0], mappedDrift[1]);
      context.lineTo(mappedEnd[0], mappedEnd[1]);
    });
    context.strokeStyle = state.colors.routeColors[blimpIndex];
    context.globalAlpha = 0.72;
    context.lineWidth = 2;
    context.setLineDash([8, 5]);
    context.stroke();
    context.setLineDash([]);
    context.globalAlpha = 1;
  }

  function drawActualTrail(blimpIndex, dimensions) {
    var tasks = state.schedule.tasks[blimpIndex];
    var scenario = state.scenario[blimpIndex];
    var route = routes[blimpIndex];
    context.beginPath();
    var mappedStart = mapPoint(route.points[0], dimensions);
    context.moveTo(mappedStart[0], mappedStart[1]);

    for (var taskIndex = 0; taskIndex < tasks.length; taskIndex += 1) {
      var task = tasks[taskIndex];
      var geometry = scenario[taskIndex];
      if (state.time <= task.arrival) break;
      if (state.time < task.gamma) {
        var firstPosition = pointOnTask(blimpIndex, taskIndex, state.time, false).point;
        var mappedFirst = mapPoint(firstPosition, dimensions);
        context.lineTo(mappedFirst[0], mappedFirst[1]);
        break;
      }
      var mappedDrift = mapPoint(geometry.actualDrift, dimensions);
      context.lineTo(mappedDrift[0], mappedDrift[1]);
      if (state.time < task.exit) {
        var secondPosition = pointOnTask(blimpIndex, taskIndex, state.time, false).point;
        var mappedSecond = mapPoint(secondPosition, dimensions);
        context.lineTo(mappedSecond[0], mappedSecond[1]);
        break;
      }
      var mappedEnd = mapPoint(geometry.end, dimensions);
      context.lineTo(mappedEnd[0], mappedEnd[1]);
    }

    context.strokeStyle = state.colors.routeColors[blimpIndex];
    context.globalAlpha = blimpIndex === state.selectedBlimp ? 1 : 0.58;
    context.lineWidth = blimpIndex === state.selectedBlimp ? 3 : 1.7;
    context.stroke();
    context.globalAlpha = 1;
  }

  function drawArrow(startX, startY, endX, endY, color) {
    var angle = Math.atan2(endY - startY, endX - startX);
    var headLength = 6;
    context.strokeStyle = color;
    context.fillStyle = color;
    context.lineWidth = 1.5;
    context.beginPath();
    context.moveTo(startX, startY);
    context.lineTo(endX, endY);
    context.stroke();
    context.beginPath();
    context.moveTo(endX, endY);
    context.lineTo(endX - headLength * Math.cos(angle - Math.PI / 6), endY - headLength * Math.sin(angle - Math.PI / 6));
    context.lineTo(endX - headLength * Math.cos(angle + Math.PI / 6), endY - headLength * Math.sin(angle + Math.PI / 6));
    context.closePath();
    context.fill();
  }

  function drawBlimp(blimpIndex, position, dimensions) {
    var mapped = mapPoint(position.point, dimensions);
    var selected = blimpIndex === state.selectedBlimp;
    var length = selected ? 19 : 15;
    var height = selected ? 9 : 7;
    context.save();
    context.translate(mapped[0], mapped[1]);
    context.rotate(position.angle);
    context.fillStyle = state.colors.routeColors[blimpIndex];
    context.strokeStyle = state.colors.background;
    context.lineWidth = 1.5;
    context.beginPath();
    context.ellipse(0, 0, length / 2, height / 2, 0, 0, TWO_PI);
    context.fill();
    context.stroke();
    context.beginPath();
    context.moveTo(-length / 2 + 2, 0);
    context.lineTo(-length / 2 - 4, -4);
    context.lineTo(-length / 2 - 4, 4);
    context.closePath();
    context.fill();
    context.restore();

    context.fillStyle = state.colors.text;
    context.font = selected ? '700 11px Lato, sans-serif' : '600 10px Lato, sans-serif';
    context.textAlign = 'left';
    context.fillText('B' + routes[blimpIndex].id, mapped[0] + 10, mapped[1] - 8);

    if (state.showWind && state.time < state.schedule.totalTime) {
      var taskIndex = position.taskIndex;
      var wind = state.scenario[blimpIndex][taskIndex].actualWind;
      var start = routes[blimpIndex].points[taskIndex];
      var end = routes[blimpIndex].points[taskIndex + 1];
      var pathCanvasAngle = Math.atan2(end[1] - start[1], end[0] - start[0]);
      var windCanvasAngle = pathCanvasAngle + wind.beta;
      var arrowLength = 9 + wind.magnitude * 2.1;
      drawArrow(
        mapped[0] + 4,
        mapped[1] + 10,
        mapped[0] + 4 + Math.cos(windCanvasAngle) * arrowLength,
        mapped[1] + 10 + Math.sin(windCanvasAngle) * arrowLength,
        state.colors.wind
      );
    }
  }

  function processorBadge(resourceIndex, dimensions) {
    var compact = dimensions.width < 560;
    var width = compact ? 58 : 78;
    var height = compact ? 18 : 21;
    var gap = compact ? 5 : 8;
    var right = dimensions.originX + dimensions.plotWidth;
    return {
      x: right - (2 - resourceIndex) * width - (1 - resourceIndex) * gap,
      y: Math.max(5, dimensions.originY - height - 9),
      width: width,
      height: height,
      centerX: right - (2 - resourceIndex) * width - (1 - resourceIndex) * gap + width / 2,
      centerY: Math.max(5, dimensions.originY - height - 9) + height / 2
    };
  }

  function drawProcessorCommunication(dimensions) {
    if (!state.showProcessorLinks) return;

    var assignments = RESOURCE_EFFICIENCY.map(function (_, resourceIndex) {
      return {
        resourceIndex: resourceIndex,
        record: activeResourceRecord(resourceIndex, state.time),
        badge: processorBadge(resourceIndex, dimensions)
      };
    });

    context.save();
    assignments.forEach(function (assignment) {
      var record = assignment.record;
      if (!record) return;

      var position = currentPosition(record.blimpIndex, state.time, false);
      var mapped = mapPoint(position.point, dimensions);
      var badge = assignment.badge;
      var duration = Math.max(0.001, record.gamma - record.start);
      var fadeWindow = Math.min(0.45, duration / 3);
      var fadeIn = clamp((state.time - record.start) / Math.max(0.001, fadeWindow), 0, 1);
      var fadeOut = clamp((record.gamma - state.time) / Math.max(0.001, fadeWindow), 0, 1);
      var pulse = 0.86 + 0.14 * Math.sin((state.time - record.start) * 8);

      context.beginPath();
      context.moveTo(mapped[0], mapped[1]);
      context.lineTo(badge.centerX, badge.centerY);
      context.strokeStyle = state.colors.routeColors[record.blimpIndex];
      context.lineWidth = 1.05;
      context.globalAlpha = (0.12 + 0.22 * Math.min(fadeIn, fadeOut)) * pulse;
      context.setLineDash([4, 5]);
      context.stroke();
      context.setLineDash([]);
    });

    assignments.forEach(function (assignment) {
      var record = assignment.record;
      var badge = assignment.badge;
      var activeColor = record ? state.colors.routeColors[record.blimpIndex] : state.colors.borderStrong;

      context.globalAlpha = 1;
      context.fillStyle = state.colors.background;
      context.fillRect(badge.x, badge.y, badge.width, badge.height);
      if (record) {
        context.globalAlpha = 0.1;
        context.fillStyle = activeColor;
        context.fillRect(badge.x, badge.y, badge.width, badge.height);
      }
      context.globalAlpha = record ? 0.75 : 0.9;
      context.strokeStyle = activeColor;
      context.lineWidth = 1;
      context.strokeRect(badge.x, badge.y, badge.width, badge.height);
      context.globalAlpha = 1;
      context.fillStyle = record ? activeColor : state.colors.muted;
      context.font = (dimensions.width < 560 ? '600 9px ' : '700 10px ') + 'Lato, sans-serif';
      context.textAlign = 'center';
      context.textBaseline = 'middle';
      context.fillText(
        'P-' + (assignment.resourceIndex === 0 ? 'A' : 'B') + (record ? ' · B' + routes[record.blimpIndex].id : ' · idle'),
        badge.centerX,
        badge.centerY
      );
    });
    context.restore();
  }

  function draw() {
    var dimensions = layout();
    context.clearRect(0, 0, dimensions.width, dimensions.height);
    drawGrid(dimensions);
    drawNominalRoutes(dimensions);
    drawMonteCarlo(dimensions);
    drawPredictedRoute(dimensions);
    routes.forEach(function (_, blimpIndex) { drawActualTrail(blimpIndex, dimensions); });
    drawProcessorCommunication(dimensions);
    routes.forEach(function (_, blimpIndex) {
      drawBlimp(blimpIndex, currentPosition(blimpIndex, state.time, false), dimensions);
    });
    updateReadout();
  }

  function setText(selector, value) {
    var element = root.querySelector(selector);
    if (element) element.textContent = value;
  }

  function activeResourceRecord(resourceIndex, time) {
    var records = state.schedule.resourceRecords[resourceIndex];
    for (var index = 0; index < records.length; index += 1) {
      if (time >= records[index].start && time < records[index].gamma) return records[index];
    }
    return null;
  }

  function updateReadout() {
    var blimpIndex = state.selectedBlimp;
    var taskIndex = currentTaskIndex(blimpIndex, state.time);
    var task = state.schedule.tasks[blimpIndex][taskIndex];
    var wind = state.scenario[blimpIndex][taskIndex].actualWind;
    var taskState;
    if (state.time >= state.schedule.totalTime) taskState = 'Inspection complete';
    else if (state.time < task.start) taskState = 'Waiting for resource';
    else if (state.time < task.gamma) taskState = 'Estimating wind';
    else taskState = 'Corrected tracking';

    setText('[data-value="time"]', state.time.toFixed(1) + ' s / ' + state.schedule.totalTime.toFixed(1) + ' s');
    setText('[data-value="task"]', 'Blimp ' + routes[blimpIndex].id + ' · Cell ' + (taskIndex + 1));
    setText('[data-value="state"]', taskState);
    setText('[data-value="wind"]', wind.magnitude.toFixed(2) + ' cm/s · ' + (wind.beta * 180 / Math.PI).toFixed(1) + '°');
    setText('[data-value="delay"]', task.delay.toFixed(2) + ' s');
    setText('[data-value="bound"]', task.bound.toFixed(2) + ' s');
    setText('[data-value="seed"]', String(state.seed));

    RESOURCE_EFFICIENCY.forEach(function (_, resourceIndex) {
      var record = activeResourceRecord(resourceIndex, state.time);
      var resourceElement = root.querySelector('[data-resource="' + resourceIndex + '"]');
      if (!resourceElement) return;
      resourceElement.textContent = record
        ? 'Blimp ' + routes[record.blimpIndex].id + ' · Cell ' + (record.taskIndex + 1)
        : 'Idle';
      resourceElement.classList.toggle('is-idle', !record);
    });
    updateScheduleMarker();
  }

  function rebuildScenario() {
    state.scenario = createScenario(state.seed);
    state.ensembles = createEnsembles(state.seed, state.sampleCount);
    state.time = 0;
    state.lastFrame = performance.now();
    draw();
  }

  function setPlaying(playing) {
    state.playing = playing;
    toggleButton.textContent = playing ? 'Pause' : 'Play';
    state.lastFrame = performance.now();
  }

  toggleButton.addEventListener('click', function () {
    setPlaying(!state.playing);
    statusElement.textContent = state.playing ? 'Simulation playing.' : 'Simulation paused.';
  });

  restartButton.addEventListener('click', function () {
    state.time = 0;
    setPlaying(true);
    statusElement.textContent = 'Inspection run restarted from t = 0.';
    draw();
  });

  resampleButton.addEventListener('click', function () {
    state.seed = (Math.imul(state.seed, 1664525) + 1013904223) >>> 0;
    rebuildScenario();
    setPlaying(true);
    statusElement.textContent = 'New independent wind realization generated with seed ' + state.seed + '.';
  });

  blimpSelect.addEventListener('change', function () {
    state.selectedBlimp = Number(blimpSelect.value);
    statusElement.textContent = 'Focused on Blimp ' + routes[state.selectedBlimp].id + '.';
    draw();
  });

  speedSelect.addEventListener('change', function () {
    state.speed = Number(speedSelect.value);
    state.lastFrame = performance.now();
    statusElement.textContent = 'Playback speed set to ' + state.speed + '×.';
  });

  sampleSelect.addEventListener('change', function () {
    state.sampleCount = Number(sampleSelect.value);
    state.ensembles = createEnsembles(state.seed, state.sampleCount);
    statusElement.textContent = state.sampleCount + ' independent wind-sample trajectories displayed.';
    draw();
  });

  layerInputs.forEach(function (input) {
    input.addEventListener('change', function () {
      var layer = input.getAttribute('data-layer');
      if (layer === 'predicted') state.showPredicted = input.checked;
      if (layer === 'montecarlo') state.showMonteCarlo = input.checked;
      if (layer === 'wind') state.showWind = input.checked;
      if (layer === 'processor-links') state.showProcessorLinks = input.checked;
      draw();
    });
  });

  function frame(now) {
    var elapsed = Math.min(0.1, (now - state.lastFrame) / 1000);
    state.lastFrame = now;
    if (state.playing) {
      state.time += elapsed * BASE_PLAYBACK_RATE * state.speed;
      if (state.time >= state.schedule.totalTime) {
        state.time = state.schedule.totalTime;
        setPlaying(false);
        statusElement.textContent = 'Inspection run complete.';
      }
      draw();
    }
    requestAnimationFrame(frame);
  }

  state.colors = readColors();
  state.schedule = buildSchedule();
  state.scenario = createScenario(state.seed);
  state.ensembles = createEnsembles(state.seed, state.sampleCount);
  renderScheduleChart();

  if (window.matchMedia('(prefers-reduced-motion: reduce)').matches) {
    setPlaying(false);
    statusElement.textContent = 'Animation paused because reduced motion is enabled.';
  }

  if ('ResizeObserver' in window) {
    new ResizeObserver(draw).observe(canvas.parentElement);
  } else {
    window.addEventListener('resize', draw);
  }

  draw();
  requestAnimationFrame(frame);
})();
