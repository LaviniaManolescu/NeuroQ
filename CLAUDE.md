# NeuroQ: Quantum Genome Detectives (Qiskit Fall Fest 2026)

Proiect de phasing de haplotipuri cu QAOA. Citește acest fișier la începutul fiecărei sesiuni, apoi documentele de mai jos înainte de a scrie cod.

## Documente de citit

`docs/arhitectura-detaliata.md` descrie modulele, contractele dintre ele (pachetul `core`), codificările, experimentele și ordinea de construcție. `docs/pipeline.md` arată cum curg datele, iar `docs/module.md` arată ce module există și ce își transmit. Împreună sunt sursa de adevăr pentru structura codului. Dacă o cerere intră în conflict cu el, întreabă înainte să schimbi structura.

`docs/QISKIT_API_DIN_CURS.md` listează apelurile Qiskit folosite efectiv în materialele QGSS 2026 și R2P (primitive V2, transpilare, `apply_layout`, sesiuni, mitigare, simulatoare). Pentru orice cod Qiskit, folosește tiparele de acolo. Dacă ai nevoie de un apel care nu e în fișier și nu ești sigur de el, scrie un test minuscul și rulează-l înainte să construiești pe el; nu scrie API-uri Qiskit din memorie.

## Mediu

Organizatorii rulează cu Qiskit 2.1.0 și Aer 0.17.0 (`requirements.txt`: `qiskit>=2.1,<3`, `qiskit-aer>=0.17,<0.19`, numpy, scipy, pandas, matplotlib, networkx, Jupyter). Tot codul trebuie să ruleze pe acest mediu, pe simulator. `qiskit-ibm-runtime` nu e în listă: se instalează separat și se importă doar în modulul `execution`. Îl cer `IBMHardware` și backendurile false din `qiskit_ibm_runtime.fake_provider` (de exemplu `FakeManilaV2`, folosit în notebook), pentru că Qiskit 2.x nu mai conține backenduri false. Proiectul trebuie să pornească și fără el: `IdealSimulator` folosește doar Aer și primitivele din `qiskit.primitives` sau `qiskit_aer.primitives`. `FakeHardware` folosește runtime dacă e instalat, altfel `GenericBackendV2` din `qiskit.providers.fake_provider` (5 qubiți, hartă liniară) cu `AerSimulator.from_backend`. Nu adăuga dependențe noi fără să spui explicit care și de ce. Dacă un tipar din `docs/QISKIT_API_DIN_CURS.md` pare să ceară o versiune mai nouă, verifică-l pe 2.1 înainte.

## Reguli de cod

Fiecare modul face un singur lucru și comunică doar prin contractele din `core` (dataclass-uri). Qiskit se importă doar în modulele `solvers` și `execution`. Fără variabile globale și fără logică de experiment în solvere. Toate rulările folosesc seed fix. Clase mici, nume clare, type hints, câte un test mic pe modul. Funcțiile preluate din notebook-ul organizatorilor se marchează printr-un comentariu `# din notebook-ul organizatorilor`, pentru declarația de la predare.

Cheia IBM nu apare niciodată în cod, notebook, loguri, output sau commit. Se citește din variabilă de mediu; `.env` e în `.gitignore`. Fiecare job pe hardware își scrie backend-ul, numărul de shots și job ID-ul în `jobs.json` înainte de a aștepta rezultatul.

## Reguli de domeniu

Hamiltonianul are forma `offset + Σ h·z + Σ J·z·z`, cu `z = +1` pentru bitul 0 și `z = −1` pentru bitul 1. Un `J` pozitiv respinge (cele două variabile vor valori diferite), un `J` negativ atrage (vor valori egale). Căutăm energia minimă. Offsetul nu schimbă cine câștigă, doar valoarea afișată; energia unui șir se calculează într-un singur loc, în `IsingProblem`.

Soluția are simetrie de oglindă (inversarea tuturor biților dă aceeași împărțire, cu copiile A și B schimbate). Compararea cu adevărul trebuie să țină cont de ea, iar fixarea unui qubit o elimină.

Qiskit citește șirurile de biți cu qubitul 0 în dreapta. Conversia dintre ordinea Qiskit și ordinea proiectului se face o singură dată, în modulul `execution`, și se testează pe un exemplu cu răspuns cunoscut.

Estimator dă un număr (valoarea medie a Hamiltonianului) și se folosește în bucla de antrenare. Sampler dă bitstring-uri cu frecvențe și se folosește pentru răspunsul final și pentru încredere. Observabilul se mută pe layout cu `apply_layout` după transpilare. Pe hardware se rulează doar circuitul final, cu unghiurile antrenate pe simulator; mitigarea se măsoară, nu se presupune că ajută.

Verificarea de bază: pe exemplul din notebook (5 SNP-uri, 9 citiri), pipeline-ul trebuie să dea `01101` cu cost 1. Dacă un modul nou schimbă acest rezultat, e o greșeală de cod până se dovedește altfel.
