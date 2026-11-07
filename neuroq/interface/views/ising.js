// Etapa 3: cuplajele J cu semn și constanta; graful stă alături (views/graph.js).
(function () {
  var NeuroQ = window.NeuroQ = window.NeuroQ || {};
  NeuroQ.views = NeuroQ.views || {};
  var S = NeuroQ.svg, C = NeuroQ.colors;

  function couplingMatrix(stage) {
    var n = stage.a.length, J = [];
    for (var i = 0; i < n; i++) { J.push(new Array(n).fill(0)); }
    stage.ising.edges.forEach(function (e) { J[e.i][e.j] = e.J; });
    return J;
  }

  NeuroQ.views.ising = {
    title: 'Problema Ising',
    lead: function () {
      return 'Fiecare pereche primește J = (d − a) / 2. Roșu (J > 0) respinge: variabilele vor ' +
        'valori diferite. Albastru (J < 0) atrage: vor valori egale. Căutăm șirul cu energia minimă.';
    },
    couplingMatrix: couplingMatrix,
    render: function (container, stage) {
      var J = couplingMatrix(stage), n = J.length;
      var max = Math.max.apply(null, [1e-9].concat(stage.ising.edges.map(function (e) {
        return Math.abs(e.J);
      })));
      var row = S.html('div', { class: 'nq-row' }, container);
      var left = S.html('div', { class: 'nq-col' }, row);
      S.html('div', { class: 'nq-caption' }, left,
        'Cuplajele J, constanta = ' + S.fmt(stage.ising.offset));
      left.appendChild(S.heatmap({
        title: 'Matricea J', values: J, upperOnly: true,
        cell: Math.max(30, Math.min(46, Math.floor(520 / n))),
        rowLabels: S.variableLabels(stage), colLabels: S.variableLabels(stage),
        color: function (v) {
          return v !== 0 ? { fill: C.coupling(v), opacity: 0.2 + 0.6 * Math.abs(v) / max } : null;
        },
        label: function (v) { return v !== 0 ? S.fmt(v) : null; }
      }));
      S.html('p', { class: 'nq-note' }, left,
        'Energia unui șir = constanta + Σ J·zᵢ·zⱼ, cu z = +1 pentru bitul 0 și −1 pentru bitul 1.');
      var right = S.html('div', { class: 'nq-col' }, row);
      S.html('div', { class: 'nq-caption' }, right,
        'Graful: ' + n + ' variabile, ' + stage.ising.edges.length + ' cuplaje nenule');
      NeuroQ.views.graph.render(right, stage, {});
    }
  };
})();
