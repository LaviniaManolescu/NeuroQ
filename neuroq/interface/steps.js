// Controlul etapelor: pornește la evenimentul 'genometric:start' al paginii GenoMEtric,
// alege instanța din data.metadata.instance și afișează cele 7 etape, cu bara de pași,
// Înapoi/Înainte și săgețile de la tastatură. Pentru capturi: ?instance=notebook&step=4&frame=end
(function () {
  var NeuroQ = window.NeuroQ = window.NeuroQ || {};
  var S = NeuroQ.svg;
  var RESULTS_BASE = '../../results/'; // relativ la neuroq/interface/GenoMEtric.html
  var STEPS = ['data', 'evidence', 'ising', 'qaoa', 'top', 'copies', 'assembly'];
  var ENCODINGS = { position: 'pe poziții (un qubit pe SNP)', read: 'pe citiri (un qubit pe citire)' };

  var state = { stage: null, index: 0, frameEnd: false, app: null };

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

  function message(title, lines) {
    var root = app();
    root.textContent = '';
    var box = S.html('section', { class: 'nq-message', role: 'alert' }, root);
    S.html('h2', {}, box, title);
    lines.forEach(function (line) { S.html('p', {}, box, line); });
    var back = S.html('button', { type: 'button', class: 'primary' }, box, 'Încarcă alt fișier');
    back.addEventListener('click', function () { showUpload(true); });
    showUpload(false);
    back.focus();
  }

  function available() {
    return window.STAGES ? Object.keys(window.STAGES) : [];
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

  function go(index) {
    NeuroQ.player.stopAll();
    state.index = Math.max(0, Math.min(STEPS.length - 1, index));
    var stage = state.stage, view = NeuroQ.views[STEPS[state.index]], root = app();
    root.textContent = '';

    var head = S.html('header', { class: 'nq-head' }, root);
    var brand = S.html('div', { class: 'nq-brand' }, head, 'Geno');
    S.html('span', {}, brand, 'ME');
    brand.appendChild(document.createTextNode('tric'));
    S.html('div', { class: 'nq-meta' }, head, 'Instanța ' + stage.metadata.instance + ', codificarea ' +
      ENCODINGS[stage.metadata.encoding] + ', ' + stage.metadata.n_qubits + ' qubiți');

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
    S.html('p', { class: 'nq-lead' }, section, view.lead(stage));
    try {
      view.render(section, stage, { resultsBase: RESULTS_BASE, frameEnd: state.frameEnd });
    } catch (error) {
      S.html('p', { class: 'nq-note' }, section,
        'Etapa nu a putut fi desenată: datele din stages.js nu au forma așteptată (' + error.message + ').');
    }

    var nav = S.html('nav', { class: 'nq-nav' }, root);
    var prev = S.html('button', { type: 'button' }, nav, 'Înapoi');
    S.html('span', { class: 'nq-hint' }, nav, 'Săgețile ← și → schimbă etapa.');
    var next = S.html('button', { type: 'button', class: 'primary' }, nav,
      state.index === STEPS.length - 1 ? 'Încarcă alt fișier' : 'Înainte');
    prev.disabled = state.index === 0;
    prev.addEventListener('click', function () { go(state.index - 1); });
    next.addEventListener('click', function () {
      if (state.index === STEPS.length - 1) { showUpload(true); return; }
      go(state.index + 1);
    });
  }

  document.addEventListener('genometric:start', function (event) {
    state.index = 0;
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
    start({ name: 'link', data: { metadata: { instance: params.get('instance') } } });
  }
})();
