// Etapa 5 (top-10 și soluția pe graf) și etapa 6 (copiile, metricile, solverii).
(function () {
  var NeuroQ = window.NeuroQ = window.NeuroQ || {};
  NeuroQ.views = NeuroQ.views || {};
  var S = NeuroQ.svg, C = NeuroQ.colors;

  function topBars(container, stage) {
    var top = stage.top_samples, W = 640, ROW = 36, LEFT = 30 + 11 * top[0].bitstring.length;
    var H = top.length * ROW + 40, max = Math.max.apply(null, top.map(function (t) { return t.count; }));
    var best = Math.min.apply(null, top.map(function (t) { return t.energy; }));
    var svg = S.root(W, H, 'Cele mai probabile 10 șiruri');
    top.forEach(function (t, i) {
      var y = 10 + i * ROW, w = (W - LEFT - 150) * t.count / max;
      var chosen = t.bitstring === stage.result.bitstring;
      var optimal = Math.abs(t.energy - best) < 1e-9;
      S.text(svg, LEFT - 10, y + ROW / 2 + 2, t.bitstring, { 'font-size': 17, 'text-anchor': 'end',
             'font-weight': chosen ? 800 : 400, 'letter-spacing': 1 });
      S.el('rect', { x: LEFT, y: y + 4, width: Math.max(2, w), height: ROW - 10, rx: 4,
                     fill: optimal ? C.accent : C.bit1, 'fill-opacity': optimal ? 0.85 : 0.6,
                     stroke: chosen ? C.ink : 'none', 'stroke-width': 2.5 }, svg);
      S.text(svg, LEFT + w + 10, y + ROW / 2 + 2, t.count + ' ori, E = ' + S.fmt(t.energy) +
             (chosen ? '  (ales)' : ''), { 'font-size': 15, 'font-weight': chosen ? 700 : 400 });
    });
    container.appendChild(svg);
    S.html('p', { class: 'nq-note' }, container,
      'Lungimea barei = de câte ori a apărut șirul în 2048 de măsurători. Turcoaz = energia cea mai mică ' +
      'dintre ele; conturul negru = șirul ales.');
  }

  // Alege orientarea (A/B) care se potrivește mai bine cu adevărul.
  function aligned(decoded, truth) {
    function matches(a, b) {
      var n = 0;
      for (var i = 0; i < a.length; i++) { if (a[i] === b[i]) n++; }
      return n;
    }
    return matches(decoded[1], truth[0]) > matches(decoded[0], truth[0])
      ? [decoded[1], decoded[0]] : decoded;
  }

  function copiesFigure(container, stage) {
    var decoded = stage.haplotypes.decoded, truth = stage.haplotypes.truth;
    var rows = truth
      ? (function () {
        var found = aligned(decoded, truth);
        return [['A adevărată', truth[0], null], ['A găsită', found[0], truth[0]],
                ['B adevărată', truth[1], null], ['B găsită', found[1], truth[1]]];
      })()
      : [['A găsită', decoded[0], null], ['B găsită', decoded[1], null]];
    var cols = decoded[0].length, cell = Math.max(30, Math.min(46, Math.floor(860 / cols)));
    var LEFT = 140, TOP = 64, labels = S.snpLabels(stage);
    var rotate = labels.some(function (l) { return l.length > 3; });
    var top = TOP + (rotate ? 24 : 0);
    var svg = S.root(LEFT + cols * cell + 10, top + rows.length * cell + 16, 'Copiile');
    S.text(svg, LEFT, 22, S.snpAxisTitle(stage), { 'font-size': 15, class: 'muted' });
    labels.forEach(function (l, c) {
      var x = LEFT + c * cell + cell / 2, y = top - 10;
      S.text(svg, x, y, l, rotate ? { 'font-size': 14, transform: 'rotate(-45 ' + x + ' ' + y + ')' }
        : { 'font-size': 16, 'text-anchor': 'middle' });
    });
    var wrong = 0;
    rows.forEach(function (row, r) {
      var y = top + r * cell + (r >= 2 ? 10 : 0);
      S.text(svg, LEFT - 12, y + cell / 2 + 6, row[0], { 'font-size': 16, 'text-anchor': 'end' });
      row[1].split('').forEach(function (ch, c) {
        var x = LEFT + c * cell, v = ch === '-' ? null : Number(ch);
        var diff = row[2] !== null && ch !== row[2][c];
        wrong += diff ? 1 : 0;
        S.el('rect', { x: x + 1, y: y + 1, width: cell - 2, height: cell - 2, rx: 4,
                       fill: v === null ? 'none' : C.bit(v), stroke: diff ? C.violated : C.line,
                       'stroke-width': diff ? 3.5 : 1 }, svg);
        S.text(svg, x + cell / 2, y + cell / 2 + 6, ch, { 'font-size': 16, 'text-anchor': 'middle',
               'font-weight': 600, fill: v === null ? C.muted : C.bitInk(v) });
      });
    });
    container.appendChild(svg);
    var verdict = !truth ? 'Adevărul nu este cunoscut pentru această instanță.'
      : wrong === 0 ? 'Copiile găsite sunt identice cu cele adevărate (oglinda contează ca potrivire).'
      : wrong + (wrong === 1 ? ' poziție diferă' : ' poziții diferă') + ' de adevăr (chenar roșu).';
    S.html('p', { class: 'nq-note' }, container, verdict + ' „-” = SNP nefazat.');
  }

  function cards(container, stage) {
    var m = stage.metrics, wrap = S.html('div', { class: 'nq-cards' }, container);
    [[Math.round(m.accuracy * 100) + '%', 'acuratețe față de adevăr'],
     [S.fmt(m.switch_errors), 'switch errors'],
     [S.fmt(m.mec_cost), 'MEC: valori de corectat'],
     [S.fmt(stage.metadata.qaoa_seconds) + ' s', 'timp QAOA']].forEach(function (c) {
      var card = S.html('div', { class: 'nq-card' }, wrap);
      S.html('div', { class: 'v' }, card, c[0]);
      S.html('div', { class: 'k' }, card, c[1]);
    });
  }

  var METHODS = { exact: 'exact', classical: 'clasic', greedy: 'lipire lacomă', qaoa: 'QAOA',
                  rqaoa: 'RQAOA', resources: 'doar resurse' };
  var ENCODINGS = { position: 'poziții', read: 'citiri' };

  function solverTable(container, stage) {
    var rows = ((window.RESULTS && window.RESULTS.rows) || []).filter(function (r) {
      return r.instance && r.instance.id === stage.metadata.instance;
    });
    if (!rows.length) {
      S.html('p', { class: 'nq-note' }, container,
        'results.js nu are rânduri pentru această instanță.');
      return;
    }
    var table = S.html('table', { class: 'nq-table' }, container);
    var head = S.html('tr', {}, S.html('thead', {}, table));
    ['Metoda', 'Codificarea', 'Transformări', 'Cost', 'Acuratețe', 'Timp QAOA'].forEach(function (h, i) {
      S.html('th', i >= 3 ? { class: 'num' } : {}, head, h);
    });
    var body = S.html('tbody', {}, table);
    rows.forEach(function (r) {
      var m = r.metrics, tr = S.html('tr', {}, body);
      S.html('td', {}, tr, METHODS[r.method] || r.method);
      S.html('td', {}, tr, r.encoding ? ENCODINGS[r.encoding] : '–');
      S.html('td', {}, tr, r.transforms.length ? 'qubit fixat' : '–');
      S.html('td', { class: 'num' }, tr, r.method === 'resources'
        ? m.n_variables + ' qubiți, ' + m.rzz_per_layer + ' RZZ'
        : (m.pair_cost !== undefined ? S.fmt(m.pair_cost) : '–'));
      S.html('td', { class: 'num' }, tr, m.accuracy !== undefined ? Math.round(m.accuracy * 100) + '%' : '–');
      S.html('td', { class: 'num' }, tr, m.qaoa_seconds !== undefined ? S.fmt(m.qaoa_seconds) + ' s' : '–');
    });
  }

  NeuroQ.views.top = {
    title: 'Rezultatul',
    lead: function () {
      return 'După antrenare măsurăm circuitul de 2048 de ori. Alegem șirul cu energia cea mai mică ' +
        'dintre cele observate; pe graf, muchiile pe care acest șir le încalcă apar cu roșu.';
    },
    render: function (container, stage) {
      var row = S.html('div', { class: 'nq-row' }, container);
      var left = S.html('div', { class: 'nq-col' }, row);
      S.html('div', { class: 'nq-caption' }, left, 'Cele mai probabile 10 șiruri');
      topBars(left, stage);
      var right = S.html('div', { class: 'nq-col' }, row);
      S.html('div', { class: 'nq-caption' }, right, 'Șirul ales ' + stage.result.bitstring +
        ', energia ' + S.fmt(stage.result.energy));
      NeuroQ.views.graph.render(right, stage, { bitstring: stage.result.bitstring });
    }
  };

  NeuroQ.views.copies = {
    title: 'Copiile',
    lead: function (stage) {
      return stage.metadata.encoding === 'position'
        ? 'Șirul ales este chiar copia A; copia B este complementul ei.'
        : 'Fiecare citire intră în echipa dată de bitul ei; copiile se obțin prin vot majoritar pe fiecare SNP.';
    },
    render: function (container, stage) {
      copiesFigure(container, stage);
      cards(container, stage);
      S.html('div', { class: 'nq-caption' }, container, 'Solverii pe această instanță (results.js)');
      solverTable(container, stage);
    }
  };
})();
