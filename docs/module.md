# Modulele NeuroQ

Codul este împărțit în module, fiecare cu o singură responsabilitate. Un modul primește de la cel dinainte rezultatul lui, îl prelucrează și îl dă mai departe, așa că oricine poate lucra pe modulul lui fără să cunoască restul. Săgețile arată ce se transmite între module.

```mermaid
flowchart TB
    SRCM["SOURCES<br/>aduce citirile:<br/>generator sintetic cu profiluri<br/>de tehnologie"]
    SELM["SELECTION<br/>alege fereastra de SNP-uri<br/>și, opțional, numărul<br/>de citiri păstrate"]
    ENCM["ENCODINGS<br/>transformă matricea<br/>într-o problemă cuantică:<br/>pe poziții sau pe citiri"]
    TRM["TRANSFORMS<br/>fixarea unui qubit,<br/>pruning, ferestre suprapuse"]
    SOLM["SOLVERS<br/>caută soluția:<br/>referință exactă, QAOA,<br/>RQAOA, baseline de lipire"]
    EXM["EXECUTION<br/>rulează circuitul:<br/>simulator, cip simulat,<br/>hardware IBM real"]
    EVM["EVALUATION<br/>măsoară calitatea soluției<br/>și resursele folosite,<br/>verdict cis/trans"]
    EXPM["EXPERIMENTS<br/>lista de experimente,<br/>rulate cu seed-uri fixe"]
    UIM["UI<br/>interfața web,<br/>afișează rezultatele"]

    SRCM -->|"citiri și SNP-uri"| SELM
    SELM -->|"matricea de citiri"| ENCM
    SELM -->|"matricea de citiri"| SOLM
    ENCM -->|"problemă Ising"| TRM
    TRM -->|"problemă Ising mai mică"| SOLM
    SOLM -->|"circuit"| EXM
    EXM -->|"măsurători"| SOLM
    SOLM -->|"soluția găsită"| EVM
    EXM -->|"resurse folosite"| EVM
    EVM -->|"rezultate"| EXPM
    EXPM -->|"results.json"| UIM

    subgraph LEG["Legendă"]
        direction TB
        L1["date"]
        L2["codificare și transformări"]
        L3["rezolvare și execuție"]
        L4["evaluare, experimente, interfață"]
    end

    classDef date fill:#eef2f7,stroke:#6b7c93,color:#1f2a37
    classDef enc fill:#dbe8fb,stroke:#3b6fc4,color:#14284b
    classDef solve fill:#dcf3e4,stroke:#2f8a55,color:#123b24
    classDef eval fill:#fdf0cf,stroke:#c48a12,color:#4a3305

    class SRCM,SELM,L1 date
    class ENCM,TRM,L2 enc
    class SOLM,EXM,L3 solve
    class EVM,EXPM,UIM,L4 eval
```

## Ce face fiecare modul

| modul | ce face | primește | dă mai departe |
|---|---|---|---|
| sources | Aduce citirile unei regiuni. Poate genera date sintetice care imită Illumina (citiri scurte) sau Nanopore (citiri lungi) și traduce literele A/C/G/T în 0 și 1. | nimic (începutul lanțului) | citirile și SNP-urile regiunii |
| selection | Alege câte SNP-uri intră în problemă și, opțional, păstrează aleator un număr de citiri (downsampling), de exemplu 10–15 la Nanopore. Fereastra și numărul de citiri sunt butoanele experimentelor. | citirile și SNP-urile | matricea de citiri |
| encodings | Transformă matricea în problemă cuantică, în două feluri: un qubit pe poziție (metoda organizatorilor) sau un qubit pe citire (metoda noastră). Cele două au funcții de pondere separate. | matricea de citiri | problema Ising |
| transforms | Simplifică problema: fixează un qubit, elimină legăturile slabe, taie probleme mari în ferestre suprapuse. | problema Ising | aceeași problemă, mai mică |
| solvers | Caută soluția: enumerare exactă (referința), QAOA, RQAOA și baseline-ul clasic de lipire a citirilor. | problema Ising (baseline-ul: matricea de citiri) | soluția găsită |
| execution | Rulează circuitele pe simulator, pe un cip IBM simulat cu zgomot sau pe hardware real. Poate și doar să le compileze, ca să numărăm porțile fără să le rulăm. | circuitul | măsurători și resursele folosite (qubiți, porți, adâncime) |
| evaluation | Măsoară cât de bună e soluția (contradicții rămase, acuratețe față de adevăr, salturi greșite între copii) și cât a costat. Dă verdictul cis/trans. | soluția și resursele | rezultate |
| experiments | Definește lista de experimente (ce date, ce codificare, ce metodă) și le rulează cu seed-uri fixe, ca rezultatele să se poată reproduce. | rezultatele | results.json |
| ui | Interfața web. Citește doar results.json, deci nu depinde de Qiskit. | results.json | — |
