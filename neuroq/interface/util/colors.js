// Culorile folosite la fel în toate etapele. Valorile vin din variabilele CSS ale paginii.
(function () {
  var NeuroQ = window.NeuroQ = window.NeuroQ || {};

  NeuroQ.colors = {
    teamA: 'var(--helix-a)',      // copia A = prima catenă a helixului
    teamB: 'var(--helix-b)',      // copia B = a doua catenă
    suspect: 'var(--nq-suspect)', // citire care nu se potrivește clar cu nicio copie
    same: 'var(--nq-same)',       // a („la fel”) și J < 0 (atrage)
    opposite: 'var(--nq-opp)',    // d („opus”) și J > 0 (respinge)
    violated: 'var(--err)',
    bit0: 'var(--nq-bit0)',
    bit1: 'var(--nq-bit1)',
    bit1Ink: 'var(--nq-bit1-ink)',
    ink: 'var(--ink)',
    muted: 'var(--muted)',
    line: 'var(--line)',
    accent: 'var(--accent)',

    team: function (team) {
      return team === 0 ? this.teamA : team === 1 ? this.teamB : this.suspect;
    },
    coupling: function (J) {
      return J > 0 ? this.opposite : this.same;
    },
    bit: function (value) {
      return value === 1 ? this.bit1 : this.bit0;
    },
    bitInk: function (value) {
      return value === 1 ? this.bit1Ink : this.ink;
    }
  };
})();
