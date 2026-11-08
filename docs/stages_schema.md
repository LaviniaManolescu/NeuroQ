# Schema `results/stages.js`

`stages.js` conține etapele demonstrației pentru interfață, calculate în Python, ca pagina să le poată afișa fără server. Îl scrie `python -m neuroq.experiments.export_stages`, care refolosește `run_stages` din `neuroq/experiments/stages.py`. Fișierul e scris atomic, în UTF-8, și conține doar JSON simplu: fără `NaN`, fără `Infinity`, fără tipuri numpy.

```js
window.STAGES = { "notebook": {...}, "illumina": {...}, "nanopore": {...} };
```

Fiecare instanță rulează pe codificarea care încape în simulator, fără fixarea unui qubit: `notebook` pe poziții, `illumina` pe poziții, `nanopore` pe citiri. Rezultatul e reproductibil (seed-uri fixe); între două exporturi diferă doar `metadata.qaoa_seconds`. Datele de intrare ale acelorași instanțe sunt în `results/inputs/<instanță>.json`, în formatul interfeței (`docs/source_format.md`), cu `metadata.instance` setat.

Indicii (SNP-uri, citiri, variabile) sunt de la 0 în tot fișierul, iar șirurile de biți sunt în ordinea proiectului: primul caracter din stânga este variabila 0.

## Câmpurile unei instanțe

| câmp | tip | ce înseamnă |
|---|---|---|
| `metadata.instance` | text | `"notebook"`, `"illumina"` sau `"nanopore"` |
| `metadata.encoding` | text | `"position"` (un qubit pe SNP) sau `"read"` (un qubit pe citire) |
| `metadata.n_qubits` | întreg | qubiții circuitului QAOA |
| `metadata.n_reads`, `metadata.n_snps` | întreg | dimensiunile matricei de citiri |
| `metadata.qaoa_seconds` | număr | durata antrenării și eșantionării QAOA, în secunde (singurul câmp care variază între rulări) |
| `evidence.values` | matrice citiri × SNP | `0`, `1` sau `null` (citirea nu acoperă SNP-ul) |
| `evidence.read_ids`, `evidence.snp_ids` | listă de text | identificatorii rândurilor și coloanelor |
| `evidence.snp_positions_bp` | listă de întregi | coordonata fiecărui SNP, în perechi de baze |
| `evidence.read_spans_bp` | listă de `[start, end]` | pentru fiecare citire, coordonatele primului și ultimului SNP observat. **Nu** sunt capetele reale ale citirii: generatorul nu păstrează lungimea citirilor |
| `a`, `d` | matrice n × n | voturile „la fel” și „opus” pe perechi de variabile (triunghiul de sus) |
| `ising.edges` | listă de `{i, j, J}` | cuplajele nenule, cu `i < j`; `J > 0` respinge, `J < 0` atrage |
| `ising.offset` | număr | constanta; energia = `offset + Σ J·z_i·z_j`, cu `z = +1` pentru bitul 0 și `−1` pentru 1 |
| `cost_history` | listă de liste | valoarea medie ⟨H⟩ la fiecare evaluare COBYLA, câte o listă pe fiecare start |
| `best_start` | întreg | indicele startului ales |
| `top_samples` | listă de `{bitstring, count, energy}` | cele mai frecvente 10 șiruri din 2048 de măsurători, descrescător după `count` |
| `result.bitstring`, `result.energy` | text, număr | șirul QAOA ales (cea mai mică energie observată) și energia lui |
| `haplotypes.decoded` | `[A, B]` | copiile decodate, câte un caracter pe SNP; `-` = SNP nefazat |
| `haplotypes.truth` | `[A, B]` sau `null` | perechea adevărată (doar la date sintetice) |
| `read_assignment` | listă, câte un obiect pe citire | `read_id`; `team` (`0` = copia A, `1` = copia B, `-1` = suspectă); `snps` = indicii SNP-urilor acoperite; `values` = valorile citite acolo; `errors[k]` = `true` dacă valoarea diferă de copia decodată a echipei (pentru o citire suspectă, de copia cea mai apropiată) |
| `metrics.accuracy` | număr între 0 și 1 | fracția de SNP-uri corecte față de adevăr, cu oglinda |
| `metrics.switch_errors` | număr | schimbările de fază greșite între SNP-uri consecutive |
| `metrics.mec_cost` | număr | câte valori trebuie corectate ca fiecare citire să se potrivească cu o copie |
| `circuit_png` | text sau `null` | calea imaginii circuitului QAOA, relativă la folderul `results/` (de exemplu `figures/stages/notebook_circuit.png`); `null` dacă circuitul are peste 12 qubiți |

## `variants`: trei variante de QAOA pe aceeași instanță

Câmpurile de la nivelul instanței (`cost_history`, `top_samples`, `result`, `haplotypes`, `circuit_png` etc.) sunt cele ale variantei `qaoa`. `variants` adaugă, pentru comutatorul de la etapele 4 și 5, cheile de mai jos; o variantă care nu se poate calcula pe o instanță lipsește (exportul tipărește motivul).

| cheie | ce este |
|---|---|
| `qaoa` | QAOA standard, identică cu câmpurile instanței |
| `qaoa_fixed` | primul qubit fixat la 0 cu `fix_variable`; problema are un qubit mai puțin, iar oglinda A↔B dispare |
| `ws_qaoa` | warm-start QAOA: circuitul pornește din relaxarea continuă a problemei |

Fiecare variantă are aceleași câmpuri (rezultatele ei, nu ale instanței):

| câmp | tip | ce înseamnă |
|---|---|---|
| `n_qubits` | întreg | qubiții problemei rezolvate (n − 1 la `qaoa_fixed`) |
| `n_rzz` | întreg | câte perechi au J ≠ 0, adică porțile RZZ ale unui strat |
| `ising` | `{edges, h, offset}` | J (muchiile nenule, cu indici ai problemei rezolvate), câmpurile locale h și constanta |
| `cost_history`, `best_start`, `qaoa_seconds` | ca la instanță | curba costului, startul ales și durata variantei |
| `circuit_png` | text sau `null` | imaginea circuitului variantei (`figures/stages/<instanță>_<variantă>_circuit.png`; `qaoa` păstrează `<instanță>_circuit.png`); `null` peste 12 qubiți |
| `all_samples` | listă de `{bitstring, count, energy}` | toate șirurile măsurate cel puțin o dată, în ordinea valorii binare a șirului; suma `count` este numărul de măsurători (2048) |
| `top_samples` | ca la instanță | cele mai frecvente 10 șiruri; prefixul lui `all_samples` sortat după `count` |
| `result` | `{bitstring, full_bitstring, energy}` | `bitstring` are lungimea problemei rezolvate (ca în `all_samples`); `full_bitstring` are lungimea problemei inițiale, cu bitul fixat reinserat (identic cu `bitstring` la celelalte variante) |
| `haplotypes`, `metrics` | ca la instanță | copiile decodate și metricile variantei |
| `read_counts` | `{a, b, suspect}` | câte citiri sunt în copia A, în copia B și suspecte; suma este `metadata.n_reads` |
| `confidence` | număr între 0 și 1 | frecvența soluției alese în cele 2048 de măsurători; la `qaoa` și `ws_qaoa` se adună și oglinda ei, la `qaoa_fixed` doar șirul ales |

Câmpuri specifice:

| câmp | doar la | ce înseamnă |
|---|---|---|
| `fixed_variables` | `qaoa_fixed` | `{indice original: bit}`; șirurile din variantă sunt mai scurte cu atâtea poziții |
| `c`, `c_regularized`, `epsilon` | `ws_qaoa` | probabilitatea de start a fiecărui qubit (P bit = 1) înainte și după regularizarea în [ε, 1 − ε] |

## Exemplu minuscul (2 SNP-uri, 2 citiri, codificare pe poziții)

```js
window.STAGES = {
 "mini": {
  "metadata": {"instance": "mini", "encoding": "position", "n_qubits": 2,
               "n_reads": 2, "n_snps": 2, "qaoa_seconds": 0.05},
  "evidence": {"values": [[0, 1], [1, null]], "read_ids": ["r1", "r2"],
               "snp_ids": ["1000", "2000"], "snp_positions_bp": [1000, 2000],
               "read_spans_bp": [[1000, 2000], [1000, 1000]]},
  "a": [[0, 0], [0, 0]],
  "d": [[0, 1], [0, 0]],
  "ising": {"edges": [{"i": 0, "j": 1, "J": 0.5}], "offset": 0.5},
  "cost_history": [[0.4, 0.1, 0.0], [0.3, 0.0]],
  "best_start": 1,
  "top_samples": [{"bitstring": "01", "count": 1030, "energy": 0.0},
                  {"bitstring": "10", "count": 1018, "energy": 0.0}],
  "result": {"bitstring": "01", "energy": 0.0},
  "haplotypes": {"decoded": ["01", "10"], "truth": ["01", "10"]},
  "read_assignment": [
   {"read_id": "r1", "team": 0, "snps": [0, 1], "values": [0, 1], "errors": [false, false]},
   {"read_id": "r2", "team": 1, "snps": [0], "values": [1], "errors": [false]}
  ],
  "metrics": {"accuracy": 1.0, "switch_errors": 0.0, "mec_cost": 0.0},
  "circuit_png": "figures/stages/mini_circuit.png"
 }
};
```
