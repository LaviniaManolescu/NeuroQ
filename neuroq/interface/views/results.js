// Etapa 5 (histograma măsurătorilor și soluția pe graf) și etapa 6 (copiile, metricile, solverii).
(function () {
  var NeuroQ = window.NeuroQ = window.NeuroQ || {};
  NeuroQ.views = NeuroQ.views || {};
  var S = NeuroQ.svg, C = NeuroQ.colors;

  // ---- Etapa 5: histograma tuturor șirurilor măsurate, cu oglinzile legate ------------------
  var MAX_BARS = 40;     // peste atâtea șiruri distincte, arătăm doar pe cele frecvente
  var RARE_SHARE = 0.01; // „frecvent” = cel puțin 1% din măsurători
  var HW = 760, HH = 420, PAD = { l: 68, r: 16, t: 62, b: 124 };
  var MONO = 'font-family:ui-monospace,Consolas,monospace';

  function mirrorOf(bits) {
    return bits.replace(/[01]/g, function (c) { return c === '0' ? '1' : '0'; });
  }

  function hamming(a, b) {
    var n = 0;
    for (var i = 0; i < a.length; i++) { if (a[i] !== b[i]) n++; }
    return n;
  }

  // Perechea (șir, oglindă) se orientează după șirul ales de QAOA: varianta cea mai apropiată de
  // el devine „șirul perechii”. Așa șirul ales dă exact read_assignment, iar cele două bare ale
  // unei perechi dau aceleași numere (A și B nu se inversează la oglindă).
  function orient(bits, chosen) {
    var other = mirrorOf(bits), d1 = hamming(bits, chosen), d2 = hamming(other, chosen);
    if (d1 !== d2) return d1 < d2 ? bits : other;
    return bits < other ? bits : other;
  }

  // Poziții: șirul e copia A, complementul e copia B; fiecare citire intră în copia cu mai
  // puține nepotriviri, iar la egalitate e suspectă (-1). Ca în encodings/position.py.
  function teamsByPositions(reads, bits) {
    var copyA = bits.split('').map(Number);
    return reads.map(function (r) {
      var missA = 0, missB = 0;
      r.snps.forEach(function (snp, j) {
        if (r.values[j] !== copyA[snp]) { missA += 1; } else { missB += 1; }
      });
      return missA < missB ? 0 : missB < missA ? 1 : -1;
    });
  }

  // Citiri: bitul i e echipa citirii i. Copiile vin din votul majoritar pe coloană al fiecărei
  // echipe (completat cu complementul celeilalte), iar o citire e suspectă dacă se potrivește
  // strict mai bine cu cealaltă copie. Ca în encodings/read.py.
  function teamsByReads(reads, bits, nSnps) {
    var teams = bits.split('').map(Number), votes = [[], []];
    for (var c = 0; c < nSnps; c++) { votes[0].push([0, 0]); votes[1].push([0, 0]); }
    reads.forEach(function (r, i) {
      r.snps.forEach(function (snp, j) { votes[teams[i]][snp][r.values[j]] += 1; });
    });
    function majority(v) { return v[1] > v[0] ? 1 : v[0] > v[1] ? 0 : null; }
    var ma = votes[0].map(majority), mb = votes[1].map(majority);
    var copies = [
      ma.map(function (a, k) { return a !== null ? a : (mb[k] !== null ? 1 - mb[k] : null); }),
      mb.map(function (b, k) { return b !== null ? b : (ma[k] !== null ? 1 - ma[k] : null); })
    ];
    function misses(r, copy) {
      var n = 0;
      r.snps.forEach(function (snp, j) { if (copy[snp] !== null && r.values[j] !== copy[snp]) n++; });
      return n;
    }
    return reads.map(function (r, i) {
      return misses(r, copies[1 - teams[i]]) < misses(r, copies[teams[i]]) ? -1 : teams[i];
    });
  }

  function pairTeams(stage, bits) {
    var oriented = orient(bits, stage.result.bitstring), reads = stage.read_assignment;
    return stage.metadata.encoding === 'position'
      ? teamsByPositions(reads, oriented)
      : teamsByReads(reads, oriented, stage.metadata.n_snps);
  }

  function tally(teams) {
    var t = { a: 0, b: 0, suspect: 0 };
    teams.forEach(function (x) { if (x === 0) { t.a++; } else if (x === 1) { t.b++; } else { t.suspect++; } });
    return t;
  }

  function niceStep(raw) {
    var pow = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / pow;
    return (f <= 1 ? 1 : f <= 2 ? 2 : f <= 5 ? 5 : 10) * pow;
  }

  function percent(count, total) {
    return (100 * count / total).toFixed(1).replace('.', ',') + '%';
  }

  function histogram(container, stage, onSelect) {
    var all = stage.all_samples, chosen = stage.result.bitstring, chosenMirror = mirrorOf(chosen);
    var total = all.reduce(function (s, x) { return s + x.count; }, 0), byBits = {};
    all.forEach(function (s) { byBits[s.bitstring] = s; });
    var shown = all;
    if (all.length > MAX_BARS) { // șirul ales și oglinda lui rămân mereu, chiar dacă sunt rare
      shown = all.filter(function (s) {
        return s.count >= RARE_SHARE * total || s.bitstring === chosen || s.bitstring === chosenMirror;
      });
    }
    var hiddenStrings = all.length - shown.length;
    var hiddenShots = total - shown.reduce(function (s, x) { return s + x.count; }, 0);

    // o culoare pe pereche, în ordinea frecvenței
    var best = {};
    shown.forEach(function (s) {
      var key = s.bitstring < mirrorOf(s.bitstring) ? s.bitstring : mirrorOf(s.bitstring);
      best[key] = Math.max(best[key] || 0, s.count);
    });
    var colorOfPair = {};
    Object.keys(best).sort(function (a, b) { return best[b] - best[a] || (a < b ? -1 : 1); })
      .forEach(function (key, rank) { colorOfPair[key] = C.pair(rank); });

    var maxCount = Math.max.apply(null, shown.map(function (s) { return s.count; }));
    var step = niceStep(maxCount / 4), yMax = step * Math.ceil(maxCount / step);
    var pw = HW - PAD.l - PAD.r, ph = HH - PAD.t - PAD.b, slot = pw / shown.length;
    var barW = Math.min(56, slot * 0.74), base = PAD.t + ph;
    var y = function (v) { return PAD.t + ph * (1 - v / yMax); };

    var wrap = S.html('div', { class: 'nq-hist' }, container);
    var svg = S.root(HW, HH, 'Histograma măsurătorilor');
    wrap.appendChild(svg);
    var tip = S.html('div', { class: 'nq-tip', role: 'tooltip', hidden: '' }, wrap);

    for (var v = 0; v <= yMax + 1e-9; v += step) {
      S.el('line', { x1: PAD.l, x2: HW - PAD.r, y1: y(v), y2: y(v), stroke: C.line,
                     'stroke-opacity': v === 0 ? 0.9 : 0.35 }, svg);
      S.text(svg, PAD.l - 8, y(v) + 5, S.fmt(v), { 'font-size': 14, 'text-anchor': 'end' });
    }
    var midY = PAD.t + ph / 2;
    S.text(svg, 18, midY, 'măsurători', { 'font-size': 15, 'text-anchor': 'middle', class: 'muted',
           transform: 'rotate(-90 18 ' + midY + ')' });

    var bars = {}, list = [];
    shown.forEach(function (s, i) {
      var x0 = PAD.l + slot * i, cx = x0 + slot / 2, h = Math.max(1.5, ph * s.count / yMax);
      var pairKey = s.bitstring < mirrorOf(s.bitstring) ? s.bitstring : mirrorOf(s.bitstring);
      var g = S.el('g', { class: 'nq-bar', tabindex: 0, role: 'button', 'data-bits': s.bitstring,
        'data-count': s.count, 'aria-label': 'șirul ' + s.bitstring + ', ' + s.count + ' măsurători' }, svg);
      S.el('rect', { class: 'nq-hit', x: x0, y: PAD.t - 6, width: slot, height: ph + 6, fill: 'transparent' }, g);
      var rect = S.el('rect', { x: cx - barW / 2, y: base - h, width: barW, height: h, rx: 2,
                                fill: colorOfPair[pairKey] }, g);
      var marker = S.el('path', { d: 'M' + cx + ' ' + (base + 4) + ' L' + (cx - 6) + ' ' + (base + 14) +
        ' L' + (cx + 6) + ' ' + (base + 14) + ' Z', fill: C.accent, visibility: 'hidden' }, g);
      var label = S.text(g, cx + 5, base + 22, s.bitstring, { 'font-size': shown.length > 20 ? 12 : 14,
        'text-anchor': 'end', transform: 'rotate(-90 ' + (cx + 5) + ' ' + (base + 22) + ')', style: MONO });
      var bar = { bits: s.bitstring, count: s.count, energy: s.energy, g: g, rect: rect,
                  marker: marker, label: label, cx: cx, top: base - h, color: colorOfPair[pairKey] };
      bars[s.bitstring] = bar;
      list.push(bar);
    });
    var link = S.el('path', { fill: 'none', 'stroke-width': 2, 'stroke-linecap': 'round',
                              'stroke-linejoin': 'round', visibility: 'hidden' }, svg);

    var selected = chosen, hot = null;
    function paint() {
      var pair = [selected, mirrorOf(selected)];
      list.forEach(function (b) {
        var inPair = pair.indexOf(b.bits) !== -1;
        var isHot = hot !== null && (b.bits === hot || b.bits === mirrorOf(hot));
        b.rect.setAttribute('fill-opacity', inPair || isHot ? 1 : 0.5);
        b.rect.setAttribute('stroke', b.bits === chosen || isHot ? C.ink : 'none');
        b.rect.setAttribute('stroke-width', b.bits === chosen ? 3.5 : 2);
        b.marker.setAttribute('visibility', inPair ? 'visible' : 'hidden');
        b.label.setAttribute('font-weight', inPair ? 800 : 400);
        b.g.setAttribute('aria-pressed', b.bits === selected ? 'true' : 'false');
      });
      var a = hot !== null ? bars[hot] : null, m = hot !== null ? bars[mirrorOf(hot)] : null;
      if (a && m) { // linie subțire între cele două bare, ca o paranteză deasupra lor
        var level = Math.min(a.top, m.top) - 16;
        link.setAttribute('d', 'M' + a.cx + ' ' + (a.top - 4) + ' V' + level + ' H' + m.cx + ' V' + (m.top - 4));
        link.setAttribute('stroke', a.color);
        link.setAttribute('visibility', 'visible');
      } else {
        link.setAttribute('visibility', 'hidden');
      }
    }

    function showTip(bar, clientX, clientY) {
      var mirror = mirrorOf(bar.bits), mirrorCount = byBits[mirror] ? byBits[mirror].count : 0;
      tip.textContent = '';
      S.html('strong', { style: MONO }, tip, bar.bits);
      S.html('div', {}, tip, bar.count + ' măsurători (' + percent(bar.count, total) + ')');
      S.html('div', {}, tip, 'energia E = ' + S.fmt(bar.energy));
      S.html('div', { class: 'nq-tip-mirror' }, tip, 'oglinda ' + mirror + ': ' + (mirrorCount
        ? mirrorCount + ' măsurători (' + percent(mirrorCount, total) + ')' : 'nemăsurată'));
      tip.hidden = false;
      var box = wrap.getBoundingClientRect(), w = tip.offsetWidth, h = tip.offsetHeight;
      var x = clientX - box.left + 14, yy = clientY - box.top + 14;
      if (x + w > box.width) { x = clientX - box.left - w - 14; }
      if (yy + h > box.height) { yy = box.height - h - 4; }
      tip.style.left = Math.max(0, x) + 'px';
      tip.style.top = Math.max(0, yy) + 'px';
    }

    function select(bits) {
      selected = bits;
      paint();
      onSelect(bits);
    }

    list.forEach(function (bar) {
      var g = bar.g;
      g.addEventListener('mouseenter', function (e) { hot = bar.bits; paint(); showTip(bar, e.clientX, e.clientY); });
      g.addEventListener('mousemove', function (e) { showTip(bar, e.clientX, e.clientY); });
      g.addEventListener('mouseleave', function () { hot = null; paint(); tip.hidden = true; });
      g.addEventListener('focus', function () {
        var box = bar.rect.getBoundingClientRect();
        hot = bar.bits; paint(); showTip(bar, box.left + box.width / 2, box.top);
      });
      g.addEventListener('blur', function () { hot = null; paint(); tip.hidden = true; });
      g.addEventListener('click', function () { select(bar.bits); });
      g.addEventListener('keydown', function (e) {
        if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); select(bar.bits); }
      });
    });

    S.html('p', { class: 'nq-note' }, wrap,
      'Înălțimea barei = de câte ori a apărut șirul în ' + total + ' de măsurători. ' +
      'Aceeași culoare = un șir și oglinda lui (aceeași împărțire a citirilor, cu A și B schimbate). ' +
      'Conturul gros = șirul ales de QAOA. Click pe o bară selectează perechea.');
    if (hiddenStrings > 0) {
      S.html('p', { class: 'nq-note' }, wrap, '+' + hiddenStrings + ' șiruri rare (sub 1% fiecare), ' +
        hiddenShots + ' din ' + total + ' măsurători, nu apar pe axă. Șirul ales și oglinda lui apar mereu.');
    }
    paint();
    return { select: select, sample: function (bits) { return byBits[bits]; }, bars: bars,
             total: total, shown: shown.length, hidden: hiddenStrings };
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
                  ws_qaoa: 'WS-QAOA', rqaoa: 'RQAOA', resources: 'doar resurse' };
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
      S.html('td', {}, tr, r.transforms.length ? 'qubit fixat'
        : r.method === 'ws_qaoa' ? 'start din relaxare' : '–');
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
      return 'După antrenare măsurăm circuitul de 2048 de ori. Fiecare bară e un șir măsurat, iar ' +
        'barele de aceeași culoare sunt oglinzi una alteia. Alege o pereche ca să vezi cum împarte ' +
        'citirile; graful arată șirul selectat, cu muchiile încălcate în roșu.';
    },
    pairTeams: pairTeams,
    orient: orient,
    render: function (container, stage) {
      var row = S.html('div', { class: 'nq-row' }, container);
      var left = S.html('div', { class: 'nq-col' }, row);
      S.html('div', { class: 'nq-caption' }, left, 'Toate șirurile măsurate');
      var right = S.html('div', { class: 'nq-col' }, row);
      var graphCaption = S.html('div', { class: 'nq-caption' }, right);
      var graphBox = S.html('div', {}, right);

      var numbers = {};
      var wrap = S.html('div', { class: 'nq-pairs' });
      S.html('div', { class: 'nq-caption' }, wrap, 'Citirile pe copii, pentru perechea selectată');
      var cardsBox = S.html('div', { class: 'nq-cards nq-pairnums' }, wrap);
      [['a', 'în copia A', C.teamA], ['b', 'în copia B', C.teamB], ['suspect', 'suspecte', C.suspect]]
        .forEach(function (def) {
          var card = S.html('div', { class: 'nq-card' }, cardsBox);
          card.style.setProperty('--nq-card-color', def[2]);
          numbers[def[0]] = S.html('div', { class: 'v' }, card, '0');
          S.html('div', { class: 'k' }, card, 'citiri ' + def[1]);
        });
      var sum = S.html('p', { class: 'nq-note' }, wrap);

      var hist = histogram(left, stage, function (bits) {
        var t = tally(pairTeams(stage, bits));
        numbers.a.textContent = t.a;
        numbers.b.textContent = t.b;
        numbers.suspect.textContent = t.suspect;
        sum.textContent = t.a + ' + ' + t.b + ' + ' + t.suspect + ' = ' + stage.metadata.n_reads +
          ' citiri. Copia A este varianta perechii cea mai apropiată de șirul ales de QAOA.';
        graphCaption.textContent = 'Șirul selectat ' + bits + ', energia ' +
          S.fmt(hist.sample(bits).energy);
        graphBox.textContent = '';
        NeuroQ.views.graph.render(graphBox, stage, { bitstring: bits });
      });
      left.appendChild(wrap);
      hist.select(stage.result.bitstring);
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
