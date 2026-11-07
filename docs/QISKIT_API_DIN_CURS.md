# Qiskit: tiparele de cod din QGSS 2026 și R2P

Acest fișier spune ce funcții și ce stil de apel s-au folosit efectiv în materialele de curs (QGSS 2026 Lab0, Lab1, Lab4a/b/c și R2P Unit 1–4). Claude Code îl citește înainte să scrie cod și îl folosește ca sursă de adevăr pentru API. Dacă un apel nu apare aici și nu ești sigur de el, se verifică în documentația oficială sau cu un test mic înainte să fie folosit în proiect.

Fiecare secțiune spune de unde vine tiparul: „din curs" înseamnă că apare în notebook-uri, „standard, neverificat în curs" înseamnă că e API normal Qiskit, dar nu l-am văzut în materiale.

## 1. Versiuni și ce nu se folosește

Materialele cer Qiskit SDK 2.x (2.0, 2.1 sau 2.5 după notebook), `qiskit-ibm-runtime` 0.40 sau mai nou, și lucrează exclusiv cu primitivele **V2**. Pachetele de bază sunt `qiskit[visualization]`, `qiskit-ibm-runtime`, `qiskit-aer`, `rustworkx`, `numpy`, `scipy`, `matplotlib`, iar pentru problemele de optimizare `qiskit-addon-opt-mapper` și `qopt-best-practices` (din GitHub: `qiskit-community/qopt-best-practices`).

**Mediul organizatorilor (`requirements.txt` din hackathon):** referința lor rulează cu Qiskit 2.1.0 și Aer 0.17.0, cu `qiskit>=2.1,<3`, `qiskit-aer>=0.17,<0.19`, numpy, scipy, pandas, matplotlib, networkx și Jupyter. Lista **nu conține `qiskit-ibm-runtime`**, `qiskit-addon-opt-mapper` sau `qopt-best-practices`. Consecințe: tot ce scriem trebuie să ruleze cu acest mediu pe simulator (Aer sau `qiskit.primitives`); partea de hardware IBM (`qiskit-ibm-runtime`) se adaugă separat și se importă doar în modulul de execuție, ca proiectul să meargă și fără ea; nu folosim opt-mapper sau qopt (Hamiltonianul se construiește direct din ponderi); `rustworkx` vine automat odată cu Qiskit, deci e disponibil fără să fie listat. Nu se cere și nu se instalează un Qiskit mai nou de 3.0, iar tiparele din acest fișier care vin din notebook-uri cu Qiskit 2.5 se verifică pe 2.1 înainte de folosire.

Nu se folosesc: `execute(...)`, `Aer.get_backend(...)`, `IBMQ`, primitivele V1 (`Sampler`/`Estimator` din vechiul `qiskit.primitives` cu `.run(circuits, observables)`), `BackendEstimator` și `BackendSampler`. Nu se scrie cod din memorie pentru ele.

**Atenție, corecție:** pentru construirea problemelor de optimizare cursul folosește `qiskit_addon_opt_mapper`, nu `qiskit-optimization` și nu `QuadraticProgram`. Vezi secțiunea 5.

## 2. Conectare și backend (din curs)

```python
from qiskit_ibm_runtime import QiskitRuntimeService
service = QiskitRuntimeService()            # token salvat o singură dată, nu în cod
backend = service.backend("ibm_torino")     # nume concret
backend = service.least_busy(min_num_qubits=127)   # sau cel mai liber
```

Tokenul IBM nu se pune niciodată în cod, notebook sau repo. În notebook-urile de curs apare `channel="ibm_quantum_platform"` și `instance=...` la salvarea contului; după salvare, `QiskitRuntimeService()` fără argumente e suficient. Tokenul se ține în variabilă de mediu sau în contul salvat local, iar repo-ul are `.env` în `.gitignore`.

Pentru dezvoltare fără hardware, cursul folosește trei variante. Atenție: `FakeTorino`, `FakeAlmadenV2` și `FakeManilaV2` vin din `qiskit_ibm_runtime.fake_provider`, deci cer `qiskit-ibm-runtime` (care nu e în `requirements.txt` al organizatorilor); în Qiskit 2.x nu mai există în `qiskit.providers.fake_provider`. Fără runtime, `GenericBackendV2(num_qubits=5, coupling_map=[[0,1],[1,2],[2,3],[3,4]], seed=...)` din `qiskit.providers.fake_provider` (testat pe Qiskit 2.1.0) ține locul unui cip fals.

```python
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime.fake_provider import FakeTorino, FakeAlmadenV2

backend = AerSimulator()                              # simulator ideal, fără zgomot
backend = FakeTorino()                                # backend fals, cu topologia și zgomotul unui cip real
backend = AerSimulator.from_backend(FakeTorino())     # simulator cu zgomotul acelui cip
```

În Unit 4 se mai folosește `StatevectorEstimator` din `qiskit.primitives` (calcul exact, fără shots) și `AerSimulator(method="matrix_product_state")` pentru circuite mai mari.

## 3. Transpilare (din curs)

```python
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

pm = generate_preset_pass_manager(backend=backend, optimization_level=3)
isa_circuit = pm.run(circuit)       # `pm.run([qc])` întoarce listă
```

În curs apare `optimization_level` 1, 2 și 3 (nivelul 3 în QAOA și în Lab4b, cu `seed_transpiler=seed` pentru reproductibilitate). Circuitul transpilat se numește „ISA circuit". Layout-ul ales de transpiler e în `isa_circuit.layout`.

**Observabilul se mută pe același layout** (apare în aproape fiecare notebook):

```python
isa_hamiltonian = hamiltonian.apply_layout(isa_circuit.layout)
```

Fără această linie, pe hardware se măsoară alți qubiți decât cei ai circuitului. Pe un simulator fără layout nu se simte, dar linia rămâne în cod ca să funcționeze la fel pe hardware.

## 4. Primitivele V2 (din curs)

Argumentul se numește `mode=` (nu `backend=` sau `session=`, în versiunile recente). Un **PUB** (primitive unified bloc) e un tuplu: pentru Estimator `(circuit, observable, parametri)`, pentru Sampler `(circuit, parametri)` sau doar `(circuit,)`.

```python
from qiskit_ibm_runtime import EstimatorV2 as Estimator, SamplerV2 as Sampler

estimator = Estimator(mode=backend)
pub = (isa_circuit, isa_hamiltonian, params)
result = estimator.run([pub]).result()
cost = result[0].data.evs            # valoarea medie; result[0].data.stds = eroarea standard

sampler = Sampler(mode=backend)
sampler.options.default_shots = 10000
result = sampler.run([(optimized_circuit,)], shots=10000).result()
counts = result[0].data.meas.get_counts()        # dict de tip {'0101': 523, ...}
counts_int = result[0].data.meas.get_int_counts()  # chei întregi
```

Detalii care contează:

- `data.meas` apare doar dacă circuitul are `measure_all()` (registrul se numește `meas`). Pentru circuitul QAOA, `measure_all()` se adaugă la copia folosită cu Sampler, nu la cea pentru Estimator.
- `data.evs` poate fi un array numpy; când funcția de cost e minimizată cu SciPy, întoarce un scalar sau un array cu un element, după caz.
- Un `.run` primește o listă de PUB-uri; se pot trimite mai multe odată.

## 5. Bucla variațională și QAOA (din curs)

Funcția de cost din notebook-urile QAOA și VQE, cu COBYLA din SciPy:

```python
from scipy.optimize import minimize
objective_func_vals = []

def cost_func_estimator(params, ansatz, hamiltonian, estimator):
    isa_hamiltonian = hamiltonian.apply_layout(ansatz.layout)
    pub = (ansatz, isa_hamiltonian, params)
    result = estimator.run([pub]).result()[0]
    cost = result.data.evs
    objective_func_vals.append(cost)       # pentru graficul de convergență
    return cost

result = minimize(
    cost_func_estimator, init_params,
    args=(candidate_circuit, cost_hamiltonian, estimator),
    method="COBYLA", tol=1e-2,
)
optimized_circuit = candidate_circuit.assign_parameters(result.x)
```

Ansatz-ul QAOA și Hamiltonianul de MaxCut din notebook-uri:

```python
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.library import QAOAAnsatz

cost_hamiltonian = SparsePauliOp.from_list(pauli_list)     # listă de (șir_pauli, coeficient)
circuit = QAOAAnsatz(cost_operator=cost_hamiltonian, reps=2)   # reps = p
circuit.measure_all()
```

Parametrii de pornire în curs sunt perechi (γ, β): `init_params = [np.pi, np.pi/2, np.pi, np.pi/2]` pentru p = 2. Alte ansatz-uri din curs: `efficient_su2` (VQE) din `qiskit.circuit.library`.

**Ordinea qubiților în șiruri Pauli (capcană):** `SparsePauliOp.from_list` citește șirul cu qubitul 0 pe poziția cea mai din dreapta. Funcția `build_max_cut_paulis` din curs construiește șirul cu qubitul i pe poziția i și apoi îl **inversează** (`"".join(paulis)[::-1]`) tocmai de aceea. La fel, la citirea rezultatului din Sampler, bitstring-ul cel mai probabil se inversează (`most_likely_bitstring.reverse()`), ca poziția i să corespundă nodului i.

**Conversia problemelor de optimizare (din Lab4b):**

```python
from qiskit_addon_opt_mapper.applications import Maxcut
from qiskit_addon_opt_mapper.translators import to_ising
from qopt_best_practices.circuit_library import annotated_qaoa_ansatz
from qopt_best_practices.transpilation import UnrollBoxes
from qopt_best_practices.transpilation.generate_preset_qaoa_pass_manager import generate_preset_qaoa_pass_manager
```

Lab4b spune explicit că `to_ising` primește un `OptimizationProblem`, **nu** obiectul `Maxcut` (se obține problema din `Maxcut`, apoi se convertește). Exercițiul în sine e lăsat necompletat în materialele pe care le avem, deci semnătura exactă și valoarea întoarsă (Hamiltonian și offset) se verifică în sursa `qiskit_addon_opt_mapper/translators/ising.py` sau în documentație înainte de folosire. Circuitul `annotated_qaoa_ansatz` se transpilează cu `pm.run(UnrollBoxes()(circuit))`. Lab4b arată și transpilarea cu `generate_preset_qaoa_pass_manager`, `SwapStrategy` și colorarea muchiilor pentru grafuri care se potrivesc direct pe topologia cipului.

Pentru Hamiltoniane scrise direct din ponderi, `SparsePauliOp.from_sparse_list([("ZZ", [i, j], w), ...], num_qubits=n)` este API standard Qiskit, dar **nu apare în curs**; se testează pe un graf mic înainte (de exemplu verificând energia unui bitstring cunoscut).

## 6. Rulare pe hardware: sesiuni, batch, opțiuni (din curs)

```python
from qiskit_ibm_runtime import Session, Batch

with Session(backend=backend) as session:
    estimator = Estimator(mode=session)
    estimator.options.default_shots = 1000
    estimator.options.dynamical_decoupling.enable = True
    estimator.options.dynamical_decoupling.sequence_type = "XY4"
    estimator.options.twirling.enable_gates = True
    estimator.options.twirling.num_randomizations = "auto"
    result = minimize(...)
```

Lab4b încearcă `Session` și, dacă nu merge (cont fără acces la sesiuni), **cade pe modul „job"**, adică `Estimator(mode=backend)`. Același reflex trebuie avut și în proiect: `try Session, except fallback la mode=backend`.

Pentru mitigare, notebook-ul Unit 3 (Error Mitigation) folosește `Batch` și compară, pe același circuit, `resilience_level = 0` cu opțiunile pornite pe rând: `dynamical_decoupling.enable` (cu `sequence_type` `"XY4"` sau `"XpXm"`), `resilience.measure_mitigation = True` (T-REx), `twirling.enable_gates = True` și `resilience.zne_mitigation = True` (ZNE). Rezultatul din curs: dynamical decoupling poate înrăutăți rezultatul când circuitul e foarte dens, deci se măsoară, nu se presupune. Pentru Sampler, opțiunile disponibile în curs sunt twirling și dynamical decoupling. `default_shots = 100000` apare în experimentele de mitigare; în bucla de optimizare, 1000.

## 7. Simulare locală, fără hardware

Cursul arată că primitivele din `qiskit_ibm_runtime` pot fi pornite pe un simulator Aer: `backend = AerSimulator(); sampler = Sampler(mode=backend)`. Același cod rulează deci pe simulator și pe hardware, schimbând doar `backend`. În mediul organizatorilor (fără runtime), simularea se face cu primitive locale, testate pe Qiskit 2.1.0: `StatevectorEstimator` și `StatevectorSampler` din `qiskit.primitives` (exacte, fără zgomot) și `SamplerV2`/`EstimatorV2` din `qiskit_aer.primitives`, create cu `SamplerV2.from_backend(AerSimulator.from_backend(backend_fals))` pentru simulare cu zgomot. `IdealSimulator` și bucla de antrenare folosesc aceste primitive; în `IBMHardware` se folosesc cele din runtime. Strategia proiectului: se optimizează pe simulator (sau pe `FakeTorino`/`AerSimulator.from_backend`), apoi un singur job final pe hardware.

## 8. Alte funcții folosite în curs, pentru recunoaștere

`QuantumCircuit`, `QuantumRegister`, `Parameter`, `ParameterVector`, `Statevector`, `plot_histogram`, `plot_bloch_multivector`, `plot_error_map`, `plot_gate_map`, `plot_circuit_layout`, `PauliEvolutionGate` și `LieTrotter` (simulare Trotter, Unit 4), `zz_feature_map` și `z_feature_map` cu `SVC` din scikit-learn (kerneluri cuantice), `mthree` (mitigarea erorilor de măsurare), `qiskit_addon_sqd` (SQD, Lab4c), `pyscf` (chimie), `rustworkx` și `networkx` (grafuri). `rustworkx` apare în notebook-urile de QAOA pentru grafuri (`rx.PyGraph`, `graph.edge_list()`).

## 9. Reguli pentru Claude Code

Se folosesc doar apelurile din acest fișier sau din documentația oficială Qiskit; dacă nu există certitudine, se scrie un test minuscul și se rulează înainte de a construi pe el. Orice apel de `apply_layout` rămâne în cod. Ordinea qubiților (little-endian) se tratează explicit și se testează pe un exemplu cu răspuns cunoscut. Tokenul IBM nu apare în cod, în loguri sau în commit-uri.
