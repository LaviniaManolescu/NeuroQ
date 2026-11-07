// Helperi mici pentru desenele SVG și etichetele comune ale etapelor.
(function () {
  var NeuroQ = window.NeuroQ = window.NeuroQ || {};
  var NS = 'http://www.w3.org/2000/svg';

  function el(tag, attrs, parent) {
    var node = document.createElementNS(NS, tag);
    Object.keys(attrs || {}).forEach(function (k) { node.setAttribute(k, attrs[k]); });
    if (parent) parent.appendChild(node);
    return node;
  }

  function text(parent, x, y, value, attrs) {
    var node = el('text', Object.assign({ x: x, y: y }, attrs || {}), parent);
    node.textContent = value;
    return node;
  }

  // Figurile mici pot crește până la 1,6× (pentru proiector), dar nu peste 760 px.
  function root(width, height, label) {
    var maxWidth = Math.max(width, Math.min(width * 1.6, 760));
    return el('svg', {
      viewBox: '0 0 ' + width + ' ' + height, class: 'nq-fig', role: 'img',
      'aria-label': label || '', preserveAspectRatio: 'xMinYMin meet',
      style: 'max-width:' + maxWidth + 'px'
    });
  }

  function html(tag, attrs, parent, content) {
    var node = document.createElement(tag);
    Object.keys(attrs || {}).forEach(function (k) { node.setAttribute(k, attrs[k]); });
    if (content !== undefined) node.textContent = content;
    if (parent) parent.appendChild(node);
    return node;
  }

  function fmt(v) {
    return Number.isInteger(v) ? String(v) : String(Math.round(v * 100) / 100);
  }

  // Heatmap: opts = {values, rowLabels, colLabels, cell, color(v) -> {fill, opacity} | null,
  //                  label(v) -> text | null, colTitle, upperOnly}
  function heatmap(opts) {
    var rows = opts.values.length, cols = rows ? opts.values[0].length : 0;
    var cell = opts.cell || 40, left = opts.left || 110, top = opts.top || 70;
    var rotate = opts.colLabels.some(function (l) { return String(l).length > 3; });
    if (rotate) top += 24;
    var svg = root(left + cols * cell + 10, top + rows * cell + 10, opts.title);
    opts.colLabels.forEach(function (label, c) {
      var x = left + c * cell + cell / 2, y = top - 10;
      text(svg, x, y, label, rotate
        ? { 'font-size': 14, transform: 'rotate(-45 ' + x + ' ' + y + ')', 'text-anchor': 'start' }
        : { 'font-size': 16, 'text-anchor': 'middle' });
    });
    if (opts.colTitle) text(svg, left, 22, opts.colTitle, { 'font-size': 15, class: 'muted' });
    opts.rowLabels.forEach(function (label, r) {
      text(svg, left - 10, top + r * cell + cell / 2 + 6, label, { 'font-size': 15, 'text-anchor': 'end' });
    });
    opts.values.forEach(function (row, r) {
      row.forEach(function (v, c) {
        var x = left + c * cell, y = top + r * cell;
        var hidden = opts.upperOnly && c <= r;
        var style = hidden ? null : opts.color(v, r, c);
        el('rect', {
          x: x + 1, y: y + 1, width: cell - 2, height: cell - 2, rx: 4,
          fill: style ? style.fill : 'none', 'fill-opacity': style ? style.opacity : 1,
          stroke: NeuroQ.colors.line, 'stroke-opacity': hidden ? 0.25 : 0.6
        }, svg);
        var label = hidden ? null : opts.label(v, r, c);
        if (label !== null && label !== undefined) {
          text(svg, x + cell / 2, y + cell / 2 + 6, label, {
            'font-size': cell >= 36 ? 16 : 13, 'text-anchor': 'middle', 'font-weight': 600,
            fill: style && style.ink ? style.ink : NeuroQ.colors.ink
          });
        }
      });
    });
    return svg;
  }

  function circleLayout(n, cx, cy, radius) {
    var points = [];
    for (var i = 0; i < n; i++) {
      var angle = Math.PI / 2 - 2 * Math.PI * i / n;
      points.push({ x: cx + radius * Math.cos(angle), y: cy - radius * Math.sin(angle) });
    }
    return points;
  }

  // Exemplul din notebook are poziții implicite (1000, 2000, ...), nu coordonate reale.
  function usesIndexAxis(stage) {
    return stage.metadata.instance === 'notebook';
  }

  function snpLabels(stage) {
    return stage.evidence.snp_positions_bp.map(function (bp, i) {
      return usesIndexAxis(stage) ? String(i + 1) : String(bp);
    });
  }

  function snpAxisTitle(stage) {
    return usesIndexAxis(stage) ? 'SNP 1, 2, 3, …' : 'poziția SNP-ului (bp)';
  }

  // Etichetele variabilelor (qubiților): SNP-uri sau citiri, după codificare.
  function variableLabels(stage) {
    if (stage.metadata.encoding === 'position') return snpLabels(stage);
    return stage.evidence.read_ids.map(function (_, i) { return 'R' + (i + 1); });
  }

  function readLabels(stage) {
    return stage.evidence.read_ids.map(function (_, i) { return 'R' + (i + 1); });
  }

  NeuroQ.svg = {
    el: el, text: text, root: root, html: html, fmt: fmt, heatmap: heatmap,
    circleLayout: circleLayout, usesIndexAxis: usesIndexAxis, snpLabels: snpLabels,
    snpAxisTitle: snpAxisTitle, variableLabels: variableLabels, readLabels: readLabels
  };
})();
