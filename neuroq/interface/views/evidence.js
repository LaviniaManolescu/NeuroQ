// Etapa 2: voturile a („la fel”) și d („opus”) pe perechi de variabile.
(function () {
  var NeuroQ = window.NeuroQ = window.NeuroQ || {};
  NeuroQ.views = NeuroQ.views || {};
  var S = NeuroQ.svg, C = NeuroQ.colors;

  function matrix(stage, values, color, caption) {
    var max = Math.max.apply(null, [1].concat(values.map(function (r) {
      return Math.max.apply(null, r);
    })));
    var n = values.length;
    var col = S.html('div', { class: 'nq-col' });
    S.html('div', { class: 'nq-caption' }, col, caption);
    col.appendChild(S.heatmap({
      title: caption, values: values, upperOnly: true,
      cell: Math.max(30, Math.min(46, Math.floor(520 / n))),
      rowLabels: S.variableLabels(stage), colLabels: S.variableLabels(stage),
      color: function (v) { return v > 0 ? { fill: color, opacity: 0.15 + 0.6 * v / max } : null; },
      label: function (v) { return v > 0 ? S.fmt(v) : null; }
    }));
    return col;
  }

  NeuroQ.views.evidence = {
    title: 'Dovezile',
    lead: function (stage) {
      var what = stage.metadata.encoding === 'position'
        ? 'pentru fiecare pereche de SNP-uri, câte citiri spun că au aceeași literă (a) sau litere diferite (d)'
        : 'pentru fiecare pereche de citiri, pe câte SNP-uri comune se potrivesc (a) sau se contrazic (d)';
      return 'Numărăm ' + what + '.';
    },
    render: function (container, stage) {
      var row = S.html('div', { class: 'nq-row' }, container);
      row.appendChild(matrix(stage, stage.a, C.same, 'a: voturi „la fel”'));
      row.appendChild(matrix(stage, stage.d, C.opposite, 'd: voturi „opus”'));
    }
  };
})();
