// Etapa 4: circuitul QAOA (imaginea din export) și curba costului, desenată progresiv.
(function () {
  var NeuroQ = window.NeuroQ = window.NeuroQ || {};
  NeuroQ.views = NeuroQ.views || {};
  var S = NeuroQ.svg, C = NeuroQ.colors;
  var W = 640, H = 380, PAD = { l: 70, r: 20, t: 20, b: 56 };
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
    S.text(svg, 16, (PAD.t + H - PAD.b) / 2, '⟨H⟩', { 'font-size': 16, 'text-anchor': 'middle' });
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

  NeuroQ.views.qaoa = {
    title: 'QAOA',
    lead: function (stage) {
      return 'Circuitul are un strat: H pe toți cei ' + stage.metadata.n_qubits + ' qubiți, o poartă ' +
        'RZZ pentru fiecare cuplaj nenul și RX pe fiecare qubit. COBYLA ajustează cele două unghiuri ' +
        'ca energia medie ⟨H⟩ să scadă, din două puncte de pornire.';
    },
    render: function (container, stage, ctx) {
      var row = S.html('div', { class: 'nq-row' }, container);
      var left = S.html('div', { class: 'nq-col' }, row);
      S.html('div', { class: 'nq-caption' }, left, 'Circuitul cu unghiurile antrenate');
      if (stage.circuit_png) {
        var img = S.html('img', { class: 'nq-circuit', src: ctx.resultsBase + stage.circuit_png,
          alt: 'Circuitul QAOA cu ' + stage.metadata.n_qubits + ' qubiți' }, left);
        img.addEventListener('error', function () {
          img.replaceWith(S.html('p', { class: 'nq-note' }, null,
            'Imaginea circuitului lipsește (' + stage.circuit_png + '). Rulează export_stages.'));
        });
      } else {
        S.html('p', {}, left, stage.metadata.n_qubits + ' qubiți, ' + stage.ising.edges.length +
          ' porți RZZ (cu doi qubiți) pe strat; circuitul e prea mare pentru desen.');
      }
      S.html('p', { class: 'nq-note' }, left, 'Timpul QAOA (antrenare + 2048 de măsurători): ' +
        S.fmt(stage.metadata.qaoa_seconds) + ' s pe simulator.');
      var right = S.html('div', { class: 'nq-col' }, row);
      S.html('div', { class: 'nq-caption' }, right, 'Costul pe iterații');
      costChart(right, stage, ctx.frameEnd);
    }
  };
})();
