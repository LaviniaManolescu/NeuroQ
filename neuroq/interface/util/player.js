// Redare pas cu pas pentru animații: Play, Pauză, Pas, viteză, Reset.
// player(total, draw, opts): draw(k) desenează starea după k pași (0..total).
(function () {
  var NeuroQ = window.NeuroQ = window.NeuroQ || {};
  var S = NeuroQ.svg;
  var SPEEDS = [['0.5', '0,5×'], ['1', '1×'], ['2', '2×'], ['4', '4×']];

  NeuroQ.player = function (container, total, draw, opts) {
    opts = opts || {};
    var baseDelay = opts.delay || 300, k = 0, timer = null, speed = 1;
    var bar = S.html('div', { class: 'nq-controls' }, container);
    var play = S.html('button', { type: 'button', class: 'primary' }, bar, 'Play');
    var step = opts.step === false ? null : S.html('button', { type: 'button' }, bar, 'Pas');
    var reset = S.html('button', { type: 'button' }, bar, 'Reset');
    var label = S.html('label', {}, bar, 'Viteza ');
    var select = S.html('select', { 'aria-label': 'Viteza animației' }, label);
    SPEEDS.forEach(function (s) {
      var o = S.html('option', { value: s[0] }, select, s[1]);
      if (s[0] === '1') o.selected = true;
    });
    var progress = S.html('span', { class: 'nq-progress', 'aria-live': 'polite' }, bar);

    function show() {
      draw(k);
      progress.textContent = (opts.unit || 'pasul') + ' ' + k + ' din ' + total;
      play.disabled = k >= total;
      if (step) step.disabled = k >= total;
    }
    function stop() {
      if (timer) { clearInterval(timer); timer = null; }
      play.textContent = 'Play';
    }
    function tick() {
      if (k >= total) { stop(); show(); return; }
      k += 1; show();
      if (k >= total) stop();
    }
    play.addEventListener('click', function () {
      if (timer) { stop(); return; }
      play.textContent = 'Pauză';
      timer = setInterval(tick, baseDelay / speed);
    });
    if (step) step.addEventListener('click', function () { stop(); tick(); });
    reset.addEventListener('click', function () { stop(); k = 0; show(); });
    select.addEventListener('change', function () {
      speed = Number(select.value);
      if (timer) { clearInterval(timer); timer = setInterval(tick, baseDelay / speed); }
    });

    k = opts.startAtEnd ? total : 0;
    show();
    var handle = { stop: stop };
    active.push(handle);
    return handle;
  };

  // Oprește toate animațiile când se schimbă etapa.
  var active = [];
  NeuroQ.player.stopAll = function () {
    active.forEach(function (p) { p.stop(); });
    active = [];
  };
})();
