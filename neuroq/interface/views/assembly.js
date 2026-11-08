// Etapa 7: asamblarea. Citirile pornesc amestecate și gri; „Rezolvă” le mută, ca un slide, în două
// panouri alăturate: copia A și copia B. Citirile suspecte rămân separat, sub panouri.
// Aceeași animație pentru orice instanță: doar datele diferă.
(function () {
  var NeuroQ = window.NeuroQ = window.NeuroQ || {};
  NeuroQ.views = NeuroQ.views || {};
  var S = NeuroQ.svg, C = NeuroQ.colors;
  var LANE = 30;                 // înălțimea unei benzi de citire
  var TITLE_Y = 26, LABEL_Y = 56, HAP_Y = 68;
  var SLIDE = 'transform 1.5s cubic-bezier(.45,.05,.25,1)';
  var FADE = 'opacity .6s';

  // Citirile ordonate după primul SNP, apoi așezate pe benzi, fără suprapunere.
  function layout(reads) {
    var order = reads.map(function (r, i) {
      return { read: r, index: i, first: Math.min.apply(null, r.snps), last: Math.max.apply(null, r.snps) };
    }).sort(function (a, b) { return a.first - b.first || a.last - b.last || a.index - b.index; });
    var laneEnds = [];
    order.forEach(function (item) {
      var lane = laneEnds.findIndex(function (end) { return end < item.first; });
      if (lane === -1) { lane = laneEnds.length; laneEnds.push(-1); }
      laneEnds[lane] = item.last;
      item.lane = lane;
    });
    return { order: order, lanes: laneEnds.length };
  }

  function plural(n) { return n === 1 ? '1 citire' : n + ' citiri'; }

  function teamName(team) { return team === 0 ? 'A' : team === 1 ? 'B' : 'suspectă'; }

  function draw(container, stage, startSolved) {
    var reads = stage.read_assignment, nSnps = stage.metadata.n_snps, labels = S.snpLabels(stage);
    var decoded = stage.haplotypes.decoded;

    // O singură figură SVG (se scalează la lățimea paginii, fără scroll orizontal): două panouri
    // și, între ele, un gol; blocul amestecat de la început stă centrat.
    var cw = Math.max(20, Math.min(40, Math.floor(1000 / (2 * nSnps + 3))));
    var panelW = nSnps * cw, gap = Math.max(48, 3 * cw), pad = 24;
    var W = 2 * pad + 2 * panelW + gap;
    var xStart = (W - panelW) / 2, xPanel = [pad, pad + panelW + gap];
    var hapH = Math.min(cw, 32), lanesTop = HAP_Y + hapH + 14;
    var digit = cw >= 36 ? 16 : 14;

    var groups = [0, 1, -1].map(function (team) {
      return reads.filter(function (r) { return team === -1 ? (r.team !== 0 && r.team !== 1) : r.team === team; });
    });
    var all = layout(reads), parts = groups.map(layout);
    var band = Math.max(all.lanes, parts[0].lanes, parts[1].lanes) * LANE;
    var suspectTop = lanesTop + band + 70;   // lasă un gol între panouri și chenarul suspectelor
    var H = (parts[2].lanes ? suspectTop + parts[2].lanes * LANE : lanesTop + band) + 22;

    var pos = {};
    all.order.forEach(function (it) { pos[it.read.read_id] = { first: it.first, last: it.last, startLane: it.lane }; });
    parts.forEach(function (part, k) {
      part.order.forEach(function (it) { pos[it.read.read_id].lane = it.lane; pos[it.read.read_id].group = k; });
    });
    function startXY(p) { return [xStart, lanesTop + p.startLane * LANE]; }
    function solvedXY(p) {
      return p.group === 2 ? [xStart, suspectTop + p.lane * LANE] : [xPanel[p.group], lanesTop + p.lane * LANE];
    }

    var bar = S.html('div', { class: 'nq-controls' }, container);
    var solve = S.html('button', { type: 'button', class: 'primary' }, bar, 'Rezolvă');
    var back = S.html('button', { type: 'button', 'aria-label': 'Înapoi: readu citirile în starea amestecată' },
      bar, 'Înapoi');
    var status = S.html('span', { class: 'nq-progress', 'aria-live': 'polite' }, bar);

    var svg = S.root(W, H, 'Asamblarea citirilor în copia A și copia B');
    var solvedEls = [], startEls = [];   // apar după „Rezolvă” / dispar după „Rezolvă”

    function axisLabels(x0, list) {
      labels.forEach(function (l, c) {
        list.push(S.text(svg, x0 + c * cw + cw / 2, LABEL_Y, l, { 'font-size': 15, 'text-anchor': 'middle',
          class: 'muted' }));
      });
    }
    var startAxis = [];
    axisLabels(xStart, startAxis);
    startAxis.forEach(function (el) { startEls.push(el); });

    [0, 1].forEach(function (k) {
      var color = C.team(k), n = parts[k].order.length;
      solvedEls.push(S.el('rect', { x: xPanel[k] - 12, y: 4, width: panelW + 24,
        height: lanesTop + band + 10 - 4, rx: 12, fill: 'none', stroke: color, 'stroke-width': 2,
        'stroke-dasharray': '7 5' }, svg));
      solvedEls.push(S.text(svg, xPanel[k], TITLE_Y, 'Copia ' + (k === 0 ? 'A' : 'B') + ': ' + plural(n),
        { 'font-size': 19, 'font-weight': 800, style: 'fill:' + color }));
      axisLabels(xPanel[k], solvedEls);
      decoded[k].split('').forEach(function (ch, c) {
        var v = ch === '-' ? null : Number(ch);
        solvedEls.push(S.el('rect', { x: xPanel[k] + c * cw + 1, y: HAP_Y, width: cw - 2, height: hapH,
          rx: 4, fill: v === null ? 'none' : C.bit(v), stroke: color, 'stroke-width': 1.5 }, svg));
        solvedEls.push(S.text(svg, xPanel[k] + c * cw + cw / 2, HAP_Y + hapH / 2 + 6, ch,
          { 'font-size': digit + 1, 'text-anchor': 'middle', 'font-weight': 700 }));
      });
    });
    if (parts[2].lanes) {
      solvedEls.push(S.el('rect', { x: xStart - 12, y: suspectTop - 36, width: panelW + 24,
        height: parts[2].lanes * LANE + 46, rx: 12, fill: 'none', stroke: C.violated, 'stroke-width': 2,
        'stroke-dasharray': '7 5' }, svg));
      solvedEls.push(S.text(svg, xStart, suspectTop - 12, 'Citiri suspecte: ' + parts[2].order.length,
        { 'font-size': 18, 'font-weight': 800, style: 'fill:' + C.violated }));
    }
    solvedEls.concat(startEls).forEach(function (el) { el.style.transition = FADE; });

    var drawn = reads.map(function (r) {
      var p = pos[r.read_id], g = S.el('g', {}, svg);
      g.style.transition = SLIDE;
      var outline = S.el('rect', { x: p.first * cw + 1, y: 2, width: (p.last - p.first + 1) * cw - 2,
        height: LANE - 4, rx: 7, fill: 'none', 'stroke-width': 2 }, g);
      outline.style.transition = 'stroke .6s';
      var cells = r.snps.map(function (snp, j) {
        var rect = S.el('rect', { x: snp * cw + 3, y: 5, width: cw - 6, height: LANE - 10, rx: 4 }, g);
        rect.style.transition = 'fill .6s, fill-opacity .6s';
        S.text(g, snp * cw + cw / 2, LANE / 2 + 5, String(r.values[j]),
          { 'font-size': digit, 'text-anchor': 'middle', 'font-weight': 600 });
        var ring = S.el('circle', { cx: snp * cw + cw / 2, cy: LANE / 2, r: Math.min(12, cw / 2 - 2),
          fill: 'none', stroke: C.violated, 'stroke-width': 2.5 }, g);
        ring.style.transition = FADE;
        return { rect: rect, ring: ring, value: r.values[j], bad: r.errors[j] };
      });
      S.el('title', {}, g).textContent = r.read_id + ': echipa ' + teamName(r.team);
      return { read: r, pos: p, g: g, outline: outline, cells: cells };
    });

    var counts = [parts[0].order.length, parts[1].order.length, parts[2].order.length];
    function apply(solved, animate) {
      if (!animate) { svg.classList.add('nq-instant'); }
      drawn.forEach(function (d) {
        var xy = solved ? solvedXY(d.pos) : startXY(d.pos), team = d.read.team;
        var suspect = team !== 0 && team !== 1;
        d.g.style.transform = 'translate(' + xy[0] + 'px,' + xy[1] + 'px)';
        d.outline.style.stroke = solved ? (suspect ? C.violated : C.team(team)) : C.suspect;
        d.outline.style.strokeDasharray = solved && suspect ? '6 4' : 'none';
        d.cells.forEach(function (cell) {
          cell.rect.style.fill = solved ? C.bit(cell.value) : C.suspect;
          cell.rect.style.fillOpacity = solved ? 1 : 0.3;
          cell.ring.style.opacity = solved && cell.bad ? 1 : 0;
        });
      });
      solvedEls.forEach(function (el) { el.style.opacity = solved ? 1 : 0; });
      startEls.forEach(function (el) { el.style.opacity = solved ? 0 : 1; });
      solve.disabled = solved;
      back.disabled = !solved;
      status.textContent = solved
        ? counts[0] + ' + ' + counts[1] + ' + ' + counts[2] + ' = ' + stage.metadata.n_reads + ' citiri'
        : stage.metadata.n_reads + ' citiri amestecate: nu știm din ce copie vine fiecare';
      if (!animate) { svg.getBoundingClientRect(); svg.classList.remove('nq-instant'); }
    }
    solve.addEventListener('click', function () { apply(true, true); });
    back.addEventListener('click', function () { apply(false, true); });

    container.appendChild(svg);
    S.html('p', { class: 'nq-note' }, container,
      'Bara unei citiri arată SNP-urile acoperite, nu lungimea reală a citirii.');
    apply(!!startSolved, false);
  }

  NeuroQ.views.assembly = {
    title: 'Asamblarea',
    lead: function (stage) {
      return 'Citirile sunt amestecate: nu știm din ce copie vine fiecare. Apasă Rezolvă ca să le vezi ' +
        'trecând în copia A sau în copia B, cu haplotipul reconstruit deasupra fiecărui panou; citirile ' +
        'care nu se potrivesc clar cu nicio copie rămân separat. ' +
        (stage.metadata.instance === 'nanopore' ? 'Aici sunt puține citiri lungi.'
          : stage.metadata.instance === 'illumina' ? 'Aici sunt multe citiri scurte.' : '');
    },
    render: function (container, stage, ctx) {
      draw(container, stage, ctx.frameEnd);
    }
  };
})();
