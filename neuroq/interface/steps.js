// Controlul etapelor: pornește la evenimentul 'genometric:start' al paginii GenoMEtric,
// alege instanța din data.metadata.instance și afișează cele 7 etape, cu bara de pași,
// Înapoi/Înainte și săgețile de la tastatură. La etapele 4 și 5 apare și comutatorul de variante QAOA.
// Pentru capturi: ?instance=notebook&step=4&frame=end&variant=ws_qaoa
(function () {
  var NeuroQ = window.NeuroQ = window.NeuroQ || {};
  var S = NeuroQ.svg;
  var RESULTS_BASE = '../../results/'; // relativ la neuroq/interface/GenoMEtric.html
  var STEPS = ['data', 'evidence', 'ising', 'qaoa', 'top', 'copies', 'assembly'];
  var ENCODINGS = { position: 'pe poziții (un qubit pe SNP)', read: 'pe citiri (un qubit pe citire)' };

  // Variantele de QAOA din stages.js (stage.variants), în ordinea afișării.
  var VARIANT_ORDER = ['qaoa', 'qaoa_fixed', 'ws_qaoa'];
  var VARIANT_LABELS = { qaoa: 'QAOA standard', qaoa_fixed: 'QAOA cu un qubit fixat',
                         ws_qaoa: 'Warm-start QAOA' };
  var VARIANT_STEPS = ['qaoa', 'top'];  // etapele 4 și 5

  // variant: varianta aleasă, păstrată la trecerea între etape; settled: curbele apar deja complete
  var state = { stage: null, index: 0, frameEnd: false, variant: 'qaoa', settled: false, app: null };

  function app() {
    if (!state.app) {
      state.app = S.html('div', { id: 'nq-app', hidden: '' }, document.body);
    }
    return state.app;
  }

  function showUpload(visible) {
    var main = document.querySelector('main');
    if (main) main.hidden = !visible;
    app().hidden = visible;
  }

  // Acasă: înapoi la ecranul de încărcare al paginii, ca să se poată alege alt fișier.
  function goHome() {
    NeuroQ.player.stopAll();
    state.stage = null;
    showUpload(true);
    window.scrollTo(0, 0);
    var input = document.getElementById('file');
    if (input) input.focus();
  }

  function message(title, lines) {
    var root = app();
    root.textContent = '';
    var box = S.html('section', { class: 'nq-message', role: 'alert' }, root);
    S.html('h2', {}, box, title);
    lines.forEach(function (line) { S.html('p', {}, box, line); });
    var back = S.html('button', { type: 'button', class: 'primary' }, box, 'Încarcă alt fișier');
    back.addEventListener('click', goHome);
    showUpload(false);
    back.focus();
  }

  function available() {
    return window.STAGES ? Object.keys(window.STAGES) : [];
  }

  function variantKeys(stage) {
    var all = stage.variants || {};
    return VARIANT_ORDER.filter(function (key) { return all[key]; });
  }

  function currentVariant(stage) {
    var keys = variantKeys(stage);
    if (!keys.length) return null;
    var key = keys.indexOf(state.variant) !== -1 ? state.variant : keys[0];
    state.variant = key;
    return { key: key, data: stage.variants[key] };
  }

  function percent(value) {
    return (100 * value).toFixed(1).replace('.', ',') + '%';
  }

  // Comutatorul și rândul de sub el; apar doar dacă sunt cel puțin două variante.
  function variantBar(section, stage, current) {
    if (variantKeys(stage).length < 2) return;
    var group = S.html('div', { class: 'nq-variants', role: 'group', 'aria-label': 'Varianta de QAOA' },
      section);
    variantKeys(stage).forEach(function (key) {
      var button = S.html('button', { type: 'button', 'aria-pressed': String(key === current.key) },
        group, VARIANT_LABELS[key]);
      button.addEventListener('click', function () {
        if (key === state.variant) return;
        state.variant = key;
        go(state.index, true);
      });
    });
    var v = current.data, stats = S.html('div', { class: 'nq-vstats' }, section);
    [[Math.round(100 * v.metrics.accuracy) + '%', 'acuratețe'],
     [String(v.n_qubits), 'qubiți'],
     [String(v.n_rzz), 'porți RZZ pe strat'],
     [percent(v.confidence), 'încrederea variantei']].forEach(function (item) {
      var cell = S.html('div', {}, stats);
      S.html('span', { class: 'v' }, cell, item[0]);
      S.html('span', { class: 'k' }, cell, item[1]);
    });
  }

  function start(detail) {
    var data = detail && detail.data;
    var instance = data && data.metadata && data.metadata.instance;
    var known = available();
    if (!window.STAGES) {
      message('Etapele nu sunt încărcate', [
        'Fișierul results/stages.js lipsește sau nu s-a putut citi.',
        'Generează-l cu: python -m neuroq.experiments.export_stages']);
      return;
    }
    if (!instance) {
      message('Fișierul nu spune ce instanță este', [
        '„' + ((detail && detail.name) || 'fișierul') + '” nu are câmpul metadata.instance.',
        'Încarcă unul dintre fișierele din results/inputs/: ' + known.join(', ') + '.']);
      return;
    }
    if (!window.STAGES[instance]) {
      message('Instanța „' + instance + '” nu are etape calculate', [
        'Instanțe disponibile: ' + known.join(', ') + '.',
        'Adaug-o în export_stages.py și rulează din nou exportul.']);
      return;
    }
    state.stage = window.STAGES[instance];
    showUpload(false);
    go(state.index);
  }

  // settled: la schimbarea variantei etapa se redesenează cu curbele deja complete
  function go(index, settled) {
    NeuroQ.player.stopAll();
    state.settled = !!settled;
    state.index = Math.max(0, Math.min(STEPS.length - 1, index));
    var stage = state.stage, view = NeuroQ.views[STEPS[state.index]], root = app();
    var current = currentVariant(stage);
    var ctx = { resultsBase: RESULTS_BASE, frameEnd: state.frameEnd || state.settled,
                variant: current ? current.data : null, variantKey: current ? current.key : null };
    root.textContent = '';

    var head = S.html('header', { class: 'nq-head' }, root);
    var brand = S.html('div', { class: 'nq-brand' }, head, 'Geno');
    S.html('span', {}, brand, 'ME');
    brand.appendChild(document.createTextNode('tric'));
    S.html('div', { class: 'nq-meta' }, head, 'Instanța ' + stage.metadata.instance + ', codificarea ' +
      ENCODINGS[stage.metadata.encoding] + ', ' + stage.metadata.n_qubits + ' qubiți');
    var home = S.html('button', { type: 'button', class: 'nq-home',
                                  title: 'Înapoi la ecranul principal, pentru alt fișier' }, head, 'Acasă');
    home.addEventListener('click', goHome);

    var list = S.html('ol', { class: 'nq-steps', 'aria-label': 'Etapele' }, root);
    STEPS.forEach(function (key, i) {
      var button = S.html('button', { type: 'button' }, S.html('li', {}, list));
      if (i === state.index) button.setAttribute('aria-current', 'step');
      S.html('span', { class: 'num' }, button, String(i + 1));
      button.appendChild(document.createTextNode(NeuroQ.views[key].title));
      button.addEventListener('click', function () { go(i); });
    });

    var section = S.html('section', { class: 'nq-stage', 'aria-live': 'polite' }, root);
    S.html('h2', {}, section, (state.index + 1) + '. ' + view.title);
    S.html('p', { class: 'nq-lead' }, section, view.lead(stage, ctx));
    if (current && VARIANT_STEPS.indexOf(STEPS[state.index]) !== -1) variantBar(section, stage, current);
    try {
      view.render(section, stage, ctx);
    } catch (error) {
      S.html('p', { class: 'nq-note' }, section,
        'Etapa nu a putut fi desenată: datele din stages.js nu au forma așteptată (' + error.message + ').');
    }

    var nav = S.html('nav', { class: 'nq-nav' }, root);
    var prev = S.html('button', { type: 'button' }, nav, 'Înapoi');
    S.html('span', { class: 'nq-hint' }, nav, 'Săgețile ← și → schimbă etapa.');
    var next = S.html('button', { type: 'button', class: 'primary' }, nav,
      state.index === STEPS.length - 1 ? 'Acasă: încarcă alt fișier' : 'Înainte');
    prev.disabled = state.index === 0;
    prev.addEventListener('click', function () { go(state.index - 1); });
    next.addEventListener('click', function () {
      if (state.index === STEPS.length - 1) { goHome(); return; }
      go(state.index + 1);
    });
    var pressed = root.querySelector('.nq-variants [aria-pressed="true"]');
    if (settled && pressed) pressed.focus();   // după o schimbare de variantă, focusul rămâne pe comutator
  }

  document.addEventListener('genometric:start', function (event) {
    state.index = 0;
    state.variant = 'qaoa';
    start(event.detail);
  });

  document.addEventListener('keydown', function (event) {
    if (!state.stage || app().hidden) return;
    var tag = (event.target && event.target.tagName) || '';
    if (tag === 'SELECT' || tag === 'INPUT') return;
    if (event.key === 'ArrowRight') { go(state.index + 1); event.preventDefault(); }
    if (event.key === 'ArrowLeft') { go(state.index - 1); event.preventDefault(); }
  });

  // Link direct (folosit la capturi): ?instance=...&step=1..7&frame=end
  var params = new URLSearchParams(window.location.search);
  if (params.has('instance')) {
    state.index = Math.max(0, (Number(params.get('step')) || 1) - 1);
    state.frameEnd = params.get('frame') === 'end';
    if (params.has('variant')) state.variant = params.get('variant');
    start({ name: 'link', data: { metadata: { instance: params.get('instance') } } });
  }
})();
