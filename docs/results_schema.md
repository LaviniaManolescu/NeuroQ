# Schema results.json

`results.json` este singurul fișier pe care îl citește interfața. Îl scrie modulul `experiments`, câte un rând (`ResultRow`, din `neuroq/core/contracts.py`) pentru fiecare combinație instanță × metodă × codificare × transformări × backend. Un fișier de probă, cu date inventate, este în `examples/results_sample.json`; un test (`tests/test_core.py`) îl încarcă și construiește un `ResultRow` din fiecare rând, deci schema de aici, fișierul de probă și dataclass-ul nu se pot despărți fără ca testul să pice.

## Învelișul fișierului

| câmp | tip | ce înseamnă |
|---|---|---|
| `schema_version` | text | versiunea schemei, acum `"0.1"`; se mărește la orice schimbare de câmpuri |
| `is_sample` | bool | `true` în fișierul de probă, `false` în rezultatele reale; interfața afișează un avertisment când e `true` |
| `note` | text | explicație liberă (la probă: că datele sunt inventate) |
| `rows` | listă de rânduri | rândurile `ResultRow`, descrise mai jos |

## Un rând (`ResultRow`)

| câmp | tip | ce înseamnă |
|---|---|---|
| `instance` | obiect | ce date s-au folosit (cheile mai jos) |
| `method` | text | `"exact"`, `"qaoa"`, `"rqaoa"` sau `"greedy"` |
| `encoding` | text sau `null` | `"position"` (un qubit pe SNP), `"read"` (un qubit pe citire); `null` pentru `greedy`, care lucrează direct pe matrice |
| `transforms` | listă de text | transformările aplicate, în ordine: `"fix_variable"`, `"prune"`, `"stitch"`; listă goală dacă niciuna |
| `metrics` | obiect | **doar numere** (cheile mai jos); o metrică care nu se aplică lipsește |
| `resources` | obiect sau `null` | resursele cuantice (`ResourceReport`); `null` la `exact` și `greedy` |
| `details` | obiect | tot ce nu e număr: haplotipuri, eșantioane, verdict (cheile mai jos); poate fi gol |

### `instance`

| cheie | tip | ce înseamnă |
|---|---|---|
| `id` | text | identificatorul instanței; același `id` = aceleași date, deci rândurile se pot compara |
| `source` | text | `"synthetic"` sau `"giab"` |
| `technology` | text | `"illumina"` sau `"nanopore"` |
| `first_snp` | întreg | indicele (de la 0) primului SNP din fereastră |
| `n_snps` | întreg | câte SNP-uri are fereastra |
| `n_reads` | întreg | câte citiri au rămas după selecție și downsampling |
| `error_rate` | număr | rata de erori din citiri (la date reale: estimată sau `null`) |
| `seed` | întreg | seed-ul folosit pentru generare și downsampling |

### `metrics` (doar numere)

| cheie | tip | ce înseamnă |
|---|---|---|
| `pair_cost` | număr | energia soluției după offset = numărul de voturi pe perechi încălcate |
| `mec_cost` | număr | câte litere trebuie corectate ca fiecare citire să se potrivească cu copia ei cea mai apropiată |
| `accuracy` | număr între 0 și 1 | fracția de SNP-uri corecte față de adevăr, ținând cont de oglindă |
| `switch_errors` | întreg | câte salturi greșite între copii sunt față de adevăr |
| `optimum_frequency` | număr între 0 și 1 | cât de des apare optimul în eșantioane (doar metode cuantice) |
| `mean_sample_cost` | număr | costul mediu pe eșantioane (doar metode cuantice) |
| `cis_trans_correct` | 0 sau 1 | dacă verdictul cis/trans e corect (doar instanțe cu V1 și V2) |

### `resources` (`ResourceReport`)

| cheie | tip | ce înseamnă |
|---|---|---|
| `n_qubits` | întreg | qubiții problemei (după fixare, dacă e cazul) |
| `two_qubit_gates_before` | întreg | porți cu doi qubiți în circuitul logic, înainte de compilare |
| `two_qubit_gates_after` | întreg sau `null` | porți cu doi qubiți după compilarea pe cip |
| `depth` | întreg sau `null` | adâncimea circuitului după compilare |
| `backend_name` | text | numele simulatorului sau al cipului |
| `job_id` | text sau `null` | doar pentru hardware real |

### `details` (ce nu e număr)

| cheie | tip | ce înseamnă |
|---|---|---|
| `haplotypes.copy_a` | text de 0 și 1 | copia A, câte un caracter pe SNP, de la stânga la dreapta (SNP-ul 1 primul) |
| `haplotypes.copy_b` | text de 0 și 1 | copia B |
| `haplotypes.read_assignment` | listă de întregi | câte o valoare pe citire: `0` = copia A, `1` = copia B, `-1` = suspectă |
| `top_samples` | listă de `{bitstring, count}` | cele mai frecvente șiruri măsurate, descrescător după `count`, în ordinea proiectului (qubitul 0 în stânga); doar metode cuantice. Lungimea șirului = `resources.n_qubits` (deci mai scurt cu unu după fixare) |
| `cis_trans_verdict` | text | `"cis"` (V1 și V2 pe aceeași copie) sau `"trans"`; doar instanțe cu V1 și V2 |

## Reguli

- Șirurile de biți sunt mereu în ordinea proiectului, niciodată în ordinea Qiskit.
- `metrics` conține doar numere (fără text, fără liste), ca interfața să le poată pune direct pe grafice.
- O cheie nouă în `metrics` sau `details` nu strică interfața; o cheie redenumită sau ștearsă cere mărirea `schema_version`.
