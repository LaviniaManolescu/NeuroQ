// Graful circular al cuplajelor: doar muchiile cu J nenul, colorate după semn.
// Cu opts.bitstring, nodurile se colorează după valoare, iar muchiile încălcate devin roșii.
(function () {
  var NeuroQ = window.NeuroQ = window.NeuroQ || {};
  NeuroQ.views = NeuroQ.views || {};
  var S = NeuroQ.svg, C = NeuroQ.colors;
  var SIZE = 520, RADIUS = 200, NODE = 22;

  function isViolated(edge, bits) {
    return bits !== null && ((edge.J > 0) === (bits[edge.i] === bits[edge.j]));
  }

  function nodeFill(stage, bit) {
    if (bit === null) return 'var(--bg2)';
    // pe citiri, bitul este echipa citirii; pe poziții, litera copiei A
    return stage.metadata.encoding === 'read' ? C.team(bit) : C.bit(bit);
  }

  NeuroQ.views.graph = {
    isViolated: isViolated,
    render: function (container, stage, opts) {
      var labels = S.variableLabels(stage), n = labels.length;
      var bits = opts.bitstring ? opts.bitstring.split('').map(Number) : null;
      var points = S.circleLayout(n, SIZE / 2, SIZE / 2, RADIUS);
      var max = Math.max.apply(null, [1e-9].concat(stage.ising.edges.map(function (e) {
        return Math.abs(e.J);
      })));
      var svg = S.root(SIZE, SIZE, 'Graful cuplajelor');
      var broken = 0;
      stage.ising.edges.forEach(function (e) {
        var bad = isViolated(e, bits);
        broken += bad ? 1 : 0;
        var a = points[e.i], b = points[e.j];
        S.el('line', {
          x1: a.x, y1: a.y, x2: b.x, y2: b.y,
          stroke: bad ? C.violated : C.coupling(e.J),
          'stroke-width': 1.5 + 5 * Math.abs(e.J) / max,
          'stroke-dasharray': e.J > 0 && !bad ? '8 6' : 'none',
          'stroke-opacity': bits && !bad ? 0.35 : 0.85
        }, svg);
        if (n <= 12) {
          var mx = (a.x + b.x) / 2, my = (a.y + b.y) / 2;
          S.el('rect', { x: mx - 17, y: my - 12, width: 34, height: 22, rx: 5, fill: 'var(--bg)',
                         'fill-opacity': 0.92 }, svg);
          S.text(svg, mx, my + 5, S.fmt(e.J), { 'font-size': 14, 'text-anchor': 'middle',
                 fill: bad ? C.violated : C.ink, 'font-weight': bad ? 700 : 400 });
        }
      });
      points.forEach(function (p, i) {
        var bit = bits ? bits[i] : null;
        S.el('circle', { cx: p.x, cy: p.y, r: NODE, fill: nodeFill(stage, bit), stroke: C.ink,
                         'stroke-width': 1.5 }, svg);
        S.text(svg, p.x, p.y + 5, labels[i].length > 4 ? 'S' + (i + 1) : labels[i], {
          'font-size': 14, 'text-anchor': 'middle', 'font-weight': 700,
          fill: bit === 1 && stage.metadata.encoding === 'position' ? C.bit1Ink : C.ink
        });
      });
      container.appendChild(svg);
      var legend = bits
        ? broken + (broken === 1 ? ' muchie încălcată' : ' muchii încălcate') + ' (roșu continuu)'
        : 'Albastru continuu: J < 0, atrage. Roșu întrerupt: J > 0, respinge. Grosimea arată |J|.';
      S.html('p', { class: 'nq-note' }, container, legend);
      return broken;
    }
  };
})();
