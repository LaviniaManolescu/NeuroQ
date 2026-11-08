// Etapa 4: circuitul QAOA (imaginea din export) și curba costului, desenată progresiv.
(function () {
  var NeuroQ = window.NeuroQ = window.NeuroQ || {};
  NeuroQ.views = NeuroQ.views || {};
  var S = NeuroQ.svg, C = NeuroQ.colors;
  var W = 640, H = 380, PAD = { l: 92, r: 20, t: 20, b: 56 };
  var START_COLORS = ['var(--accent)', 'var(--accent2)', 'var(--muted)'];

  function costChart(container, stage, startAtEnd) {
    var histories = stage.cost_history;
    var all = [].concat.apply([], histories).concat([stage.result.energy]);
    var lo = Math.min.apply(null, all), hi = Math.max.apply(null, all);
    var span = hi - lo || 1;
    lo -= span * 0.05; hi += span * 0.05;
    var longest = Math.max.apply(null, histories.map(function (h) { return h.length; }));
    var x = function (i) { return PAD.l + (W - PAD.l - PAD.r) * i / Math.max(1, longest - 1); };
    var y = function (v) { return PAD.t + (H - PAD.t - PAD.b) * (hi - v) / (hi - lo); };

    var svg = S.root(W, H, 'Costul QAOA pe iterații');
    for (var g = 0; g <= 4; g++) {
      var v = lo + (hi - lo) * g / 4;
      S.el('line', { x1: PAD.l, x2: W - PAD.r, y1: y(v), y2: y(v), stroke: C.line,
                     'stroke-opacity': 0.4 }, svg);
      S.text(svg, PAD.l - 8, y(v) + 5, S.fmt(v), { 'font-size': 14, 'text-anchor': 'end' });
    }
    S.text(svg, (PAD.l + W - PAD.r) / 2, H - 14, 'evaluarea COBYLA', { 'font-size': 15,
           'text-anchor': 'middle', class: 'muted' });
    var midY = (PAD.t + H - PAD.b) / 2;
    S.text(svg, 20, midY, 'energia medie ⟨H⟩', { 'font-size': 15, 'text-anchor': 'middle',
           class: 'muted', transform: 'rotate(-90 20 ' + midY + ')' });
    var target = y(stage.result.energy);
    S.el('line', { x1: PAD.l, x2: W - PAD.r, y1: target, y2: target, stroke: C.ink,
                   'stroke-dasharray': '6 5' }, svg);
    S.text(svg, W - PAD.r, target - 8, 'energia soluției alese = ' + S.fmt(stage.result.energy),
           { 'font-size': 14, 'text-anchor': 'end' });
    var lines = histories.map(function (_, s) {
      return S.el('polyline', { fill: 'none', stroke: START_COLORS[s % START_COLORS.length],
        'stroke-width': s === stage.best_start ? 3.5 : 2, 'stroke-linejoin': 'round' }, svg);
    });
    var legend = S.html('p', { class: 'nq-note' }, null);
    legend.textContent = histories.map(function (h, s) {
      return 'startul ' + (s + 1) + (s === stage.best_start ? ' (ales)' : '') + ': ' + h.length + ' evaluări';
    }).join('; ') + '.';
    container.appendChild(svg);
    container.appendChild(legend);

    NeuroQ.player(container, longest, function (k) {
      lines.forEach(function (line, s) {
        line.setAttribute('points', histories[s].slice(0, k).map(function (v, i) {
          return x(i) + ',' + y(v);
        }).join(' '));
      });
    }, { delay: 120, unit: 'evaluarea', step: true, startAtEnd: startAtEnd });
  }

  // Warm-start: înclinarea de start a fiecărui qubit (probabilitatea ca bitul să fie 1), înainte de
  // regularizare (relaxarea continuă) și după ea (tăiată în [ε, 1 − ε]); linia 0,5 = fără preferință.
  var TW = 640, TH = 250, TPAD = { l: 52, r: 16, t: 24, b: 66 };

  function tiltBars(container, stage, v) {
    var labels = S.variableLabels(stage), n = v.c.length;
    var slot = (TW - TPAD.l - TPAD.r) / n, bw = Math.min(22, slot * 0.38);
    var y = function (p) { return TPAD.t + (TH - TPAD.t - TPAD.b) * (1 - p); };
    var svg = S.root(TW, TH, 'Înclinarea de start a qubiților');
    [0, 0.5, 1].forEach(function (p) {
      S.el('line', { x1: TPAD.l, x2: TW - TPAD.r, y1: y(p), y2: y(p), stroke: C.line,
                     'stroke-opacity': p === 0.5 ? 0.95 : 0.4,
                     'stroke-dasharray': p === 0.5 ? '6 5' : 'none' }, svg);
      S.text(svg, TPAD.l - 8, y(p) + 5, S.fmt(p), { 'font-size': 14, 'text-anchor': 'end' });
    });
    S.text(svg, TW - TPAD.r, y(0.5) - 6, 'linia 0,5: fără preferință',
           { 'font-size': 13, 'text-anchor': 'end', class: 'muted' });
    labels.forEach(function (label, i) {
      var cx = TPAD.l + slot * (i + 0.5);
      [[v.c[i], -bw - 1, C.bit1, 0.65, 'relaxarea'], [v.c_regularized[i], 1, C.accent, 0.95, 'după ε']]
        .forEach(function (bar) {
          var top = y(bar[0]);
          var rect = S.el('rect', { x: cx + bar[1], y: top, width: bw, height: Math.max(1.5, y(0) - top),
                                    rx: 2, fill: bar[2], 'fill-opacity': bar[3] }, svg);
          S.el('title', {}, rect).textContent = 'qubitul ' + label + ', ' + bar[4] + ': ' + S.fmt(bar[0]);
        });
      S.text(svg, cx, y(0) + 20, label, { 'font-size': n > 12 ? 12 : 14, 'text-anchor': 'middle' });
    });
    [[C.bit1, 0.65, 'relaxarea continuă, înainte de ε'],
     [C.accent, 0.95, 'după regularizare, ε = ' + S.fmt(v.epsilon)]].forEach(function (item, k) {
      var x0 = TPAD.l + k * 270;
      S.el('rect', { x: x0, y: TH - 20, width: 14, height: 14, rx: 2, fill: item[0],
                     'fill-opacity': item[1] }, svg);
      S.text(svg, x0 + 20, TH - 8, item[2], { 'font-size': 13 });
    });
    container.appendChild(svg);
  }

  NeuroQ.views.qaoa = {
    title: 'QAOA',
    lead: function (stage, ctx) {
      var key = ctx && ctx.variantKey, v = ctx && ctx.variant;
      var n = v ? v.n_qubits : stage.metadata.n_qubits;
      var train = 'COBYLA ajustează cele două unghiuri ca energia medie ⟨H⟩ să scadă, din două puncte ' +
        'de pornire.';
      if (key === 'qaoa_fixed') {
        return 'Primul qubit este fixat la 0, deci circuitul are ' + n + ' qubiți, unul mai puțin. Cuplajele ' +
          'lui cu ceilalți devin câmpuri locale, adică porți RZ pe un singur qubit, iar oglinda A↔B dispare: ' +
          'soluția nu se mai împarte în două șiruri echivalente. ' + train;
      }
      if (key === 'ws_qaoa') {
        return 'Warm-start: nu pornim din suprapunerea uniformă. Rezolvăm întâi o variantă „moale” a ' +
          'problemei, în care biții pot fi între 0 și 1, iar fiecare qubit pornește înclinat spre valoarea ' +
          'găsită, printr-o rotație RY. Cele ' + v.n_rzz + ' porți RZZ rămân, iar mixerul RY-RZ-RY păstrează ' +
          'starea de start. ' + train;
      }
      return 'Circuitul are un strat: H pe toți cei ' + n + ' qubiți, o poartă RZZ pentru fiecare cuplaj ' +
        'nenul și RX pe fiecare qubit. ' + train;
    },
    render: function (container, stage, ctx) {
      var v = ctx.variant || stage;  // fără variante în stages.js, rămân câmpurile instanței
      var nQubits = ctx.variant ? v.n_qubits : stage.metadata.n_qubits;
      var nRzz = ctx.variant ? v.n_rzz : stage.ising.edges.length;
      var seconds = ctx.variant ? v.qaoa_seconds : stage.metadata.qaoa_seconds;
      var row = S.html('div', { class: 'nq-row' }, container);
      var left = S.html('div', { class: 'nq-col' }, row);
      S.html('div', { class: 'nq-caption' }, left, 'Circuitul cu unghiurile antrenate');
      if (v.circuit_png) {
        var img = S.html('img', { class: 'nq-circuit', src: ctx.resultsBase + v.circuit_png,
          alt: 'Circuitul QAOA cu ' + nQubits + ' qubiți' }, left);
        img.addEventListener('error', function () {
          img.replaceWith(S.html('p', { class: 'nq-note' }, null,
            'Imaginea circuitului lipsește (' + v.circuit_png + '). Rulează export_stages.'));
        });
      } else {
        S.html('p', {}, left, nQubits + ' qubiți, ' + nRzz +
          ' porți RZZ (cu doi qubiți) pe strat; circuitul e prea mare pentru desen.');
      }
      S.html('p', { class: 'nq-note' }, left, 'Timpul QAOA (antrenare + 2048 de măsurători): ' +
        S.fmt(seconds) + ' s pe simulator.');
      var right = S.html('div', { class: 'nq-col' }, row);
      S.html('div', { class: 'nq-caption' }, right, 'Costul pe iterații');
      costChart(right, v, ctx.frameEnd);
      if (ctx.variantKey === 'ws_qaoa') {
        S.html('div', { class: 'nq-caption' }, right, 'Înclinarea de start a fiecărui qubit (P bit = 1)');
        tiltBars(right, stage, v);
      }
    }
  };
})();
