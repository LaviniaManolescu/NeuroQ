// Etapa 1: matricea citirilor × SNP-uri (0, 1 sau gol).
(function () {
  var NeuroQ = window.NeuroQ = window.NeuroQ || {};
  NeuroQ.views = NeuroQ.views || {};
  var S = NeuroQ.svg, C = NeuroQ.colors;

  NeuroQ.views.data = {
    title: 'Datele',
    lead: function (stage) {
      var m = stage.metadata;
      return m.n_reads + ' citiri pe ' + m.n_snps + ' SNP-uri. 0 este litera de referință, ' +
        '1 cea alternativă; o celulă goală înseamnă că citirea nu acoperă acel SNP. ' +
        'Nu știm din ce copie a cromozomului vine fiecare citire.';
    },
    render: function (container, stage) {
      var values = stage.evidence.values;
      var cols = values[0].length;
      var cell = Math.max(30, Math.min(46, Math.floor(900 / cols)));
      container.appendChild(S.heatmap({
        title: 'Matricea citirilor', values: values, cell: cell,
        rowLabels: S.readLabels(stage), colLabels: S.snpLabels(stage),
        colTitle: S.snpAxisTitle(stage),
        color: function (v) {
          return v === null ? null : { fill: C.bit(v), opacity: 1, ink: C.bitInk(v) };
        },
        label: function (v) { return v === null ? null : String(v); }
      }));
    }
  };
})();
