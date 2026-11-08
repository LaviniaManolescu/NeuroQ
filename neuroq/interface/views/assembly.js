// Etapa 7: asamblarea, ca într-un browser de cromozom.
// Citirile intră pe rând pe axa SNP-urilor, colorate după echipă; copiile de jos se
// completează coloană cu coloană prin votul majoritar al citirilor intrate.
// Aceeași animație pentru orice instanță: doar datele diferă.
(function () {
  var NeuroQ = window.NeuroQ = window.NeuroQ || {};
  NeuroQ.views = NeuroQ.views || {};
  var S = NeuroQ.svg, C = NeuroQ.colors;
  var LEFT = 120, LANE = 30, GAP = 26;

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

  function majority(votes) {
    return votes[1] > votes[0] ? 1 : votes[0] > votes[1] ? 0 : null;
  }

  // Copiile după primele k citiri: votul echipei A și al echipei B pe fiecare coloană.
  function votedCopies(order, k, nSnps) {
    var votes = [[], []];
    for (var c = 0; c < nSnps; c++) { votes[0].push([0, 0]); votes[1].push([0, 0]); }
    order.slice(0, k).forEach(function (item) {
      var team = item.read.team;
      if (team !== 0 && team !== 1) return; // citirile suspecte nu votează
      item.read.snps.forEach(function (snp, j) { votes[team][snp][item.read.values[j]] += 1; });
    });
    var a = [], b = [];
    for (c = 0; c < nSnps; c++) {
      var ma = majority(votes[0][c]), mb = majority(votes[1][c]);
      a.push(ma !== null ? ma : mb !== null ? 1 - mb : null);
      b.push(mb !== null ? mb : ma !== null ? 1 - ma : null);
    }
    return [a, b];
  }

  function draw(container, stage, startAtEnd) {
    var reads = stage.read_assignment, nSnps = stage.metadata.n_snps;
    var placed = layout(reads), labels = S.snpLabels(stage);
    var cw = Math.max(28, Math.min(58, Math.floor(980 / nSnps)));
    var rotate = labels.some(function (l) { return l.length > 3; });
    var top = 56 + (rotate ? 26 : 0);
    var tracksH = placed.lanes * LANE;
    var copiesY = top + tracksH + GAP;
    var W = LEFT + nSnps * cw + 20, H = copiesY + 2 * (cw + 6) + 10;
    var svg = S.root(W, H, 'Asamblarea copiilor din citiri');
    var x = function (c) { return LEFT + c * cw; };

    S.text(svg, LEFT, 20, S.snpAxisTitle(stage), { 'font-size': 15, class: 'muted' });
    labels.forEach(function (l, c) {
      var cx = x(c) + cw / 2, cy = top - 12;
      S.text(svg, cx, cy, l, rotate ? { 'font-size': 13, transform: 'rotate(-45 ' + cx + ' ' + cy + ')' }
        : { 'font-size': 16, 'text-anchor': 'middle' });
      S.el('line', { x1: cx, x2: cx, y1: top - 4, y2: top + tracksH, stroke: C.line,
                     'stroke-opacity': 0.35 }, svg);
    });
    S.text(svg, LEFT - 12, top + 18, 'citiri', { 'font-size': 15, 'text-anchor': 'end', class: 'muted' });

    var groups = placed.order.map(function (item) {
      var g = S.el('g', {}, svg), y = top + item.lane * LANE;
      var color = C.team(item.read.team);
      S.el('rect', { x: x(item.first) + 2, y: y + 4, width: (item.last - item.first + 1) * cw - 4,
                     height: LANE - 8, rx: 6, fill: color, 'fill-opacity': 0.22, stroke: color,
                     'stroke-width': 2, class: 'nq-read' }, g);
      item.read.snps.forEach(function (snp, j) {
        var bad = item.read.errors[j];
        S.text(g, x(snp) + cw / 2, y + LANE / 2 + 5, String(item.read.values[j]), {
          'font-size': 15, 'text-anchor': 'middle', 'font-weight': bad ? 800 : 600,
          fill: bad ? C.violated : C.ink
        });
        if (bad) {
          S.el('circle', { cx: x(snp) + cw / 2, cy: y + LANE / 2, r: Math.min(12, cw / 2 - 2),
                           fill: 'none', stroke: C.violated, 'stroke-width': 2.5 }, g);
        }
      });
      var title = S.el('title', {}, g);
      title.textContent = item.read.read_id + ': echipa ' +
        (item.read.team === 0 ? 'A' : item.read.team === 1 ? 'B' : 'suspectă');
      return g;
    });

    var copyCells = [0, 1].map(function (row) {
      var y = copiesY + row * (cw + 6);
      S.text(svg, LEFT - 12, y + cw / 2 + 6, row === 0 ? 'copia A' : 'copia B', {
        'font-size': 16, 'text-anchor': 'end', 'font-weight': 700, fill: row === 0 ? C.teamA : C.teamB
      });
      var cells = [];
      for (var c = 0; c < nSnps; c++) {
        var rect = S.el('rect', { x: x(c) + 2, y: y, width: cw - 4, height: cw, rx: 5,
                                  stroke: row === 0 ? C.teamA : C.teamB, 'stroke-width': 2 }, svg);
        var label = S.text(svg, x(c) + cw / 2, y + cw / 2 + 6, '', { 'font-size': 17,
                           'text-anchor': 'middle', 'font-weight': 700 });
        cells.push({ rect: rect, label: label });
      }
      return cells;
    });

    container.appendChild(svg);
    var legend = S.html('p', { class: 'nq-note' }, container);
    legend.textContent = 'Turcoaz = citire pusă în copia A, indigo = copia B, gri = suspectă (nu votează). ' +
      'Cerc roșu = valoare citită care diferă de copia decodată (eroare de citire). ' +
      'Bara unei citiri arată SNP-urile acoperite, nu lungimea reală a citirii.';

    NeuroQ.player(container, placed.order.length, function (k) {
      groups.forEach(function (g, i) {
        g.setAttribute('opacity', i < k ? 1 : 0);
        var bar = g.querySelector('.nq-read');
        bar.setAttribute('stroke-width', i === k - 1 ? 4 : 2);
      });
      var copies = votedCopies(placed.order, k, nSnps);
      copies.forEach(function (copy, row) {
        copy.forEach(function (v, c) {
          var cell = copyCells[row][c];
          cell.rect.setAttribute('fill', v === null ? 'none' : C.bit(v));
          cell.label.textContent = v === null ? '' : String(v);
          cell.label.setAttribute('fill', v === null ? C.ink : C.bitInk(v));
        });
      });
    }, { delay: 700, unit: 'citirea', startAtEnd: startAtEnd });
  }

  NeuroQ.views.assembly = {
    title: 'Asamblarea',
    lead: function (stage) {
      return 'Citirile intră pe rând, în ordinea poziției. Fiecare votează în copia echipei ei; ' +
        'sub ele, copiile A și B se completează prin vot majoritar pe fiecare SNP. ' +
        (stage.metadata.instance === 'nanopore' ? 'Aici sunt puține citiri lungi.'
          : stage.metadata.instance === 'illumina' ? 'Aici sunt multe citiri scurte.' : '');
    },
    render: function (container, stage, ctx) {
      draw(container, stage, ctx.frameEnd);
    }
  };
})();
