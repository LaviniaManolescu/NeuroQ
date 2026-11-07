# NeuroQ · Arhitectura detaliată

Acest document explică amănunțit ce construim, ce face fiecare modul, ce își transmit modulele între ele și în ce ordine lucrăm. Este scris pentru echipă și pentru Claude Code, care îl poate citi ca să respecte arhitectura. Diagramele sunt în fișierele `module.md` (cine cu cine vorbește) și `pipeline.md` (cum curg datele). Părțile marcate „opțional" se fac doar dacă rămâne timp.

## Întrebarea investigației

Ce codificare a problemei de phasing folosește mai puține resurse cuantice (qubiți, porți cu doi qubiți, adâncime) în funcție de tipul de secvențiere, și cât ajută fixarea simetriei? Comparăm codificarea organizatorilor, cu un qubit pe poziție SNP, cu codificarea noastră, cu un qubit pe citire, pe aceleași ferestre de genom citite cu două tehnologii: citiri scurte de tip Illumina și citiri lungi de tip Nanopore. Ne așteptăm ca pozițiile să câștige la citiri scurte, iar citirile să câștige la citiri lungi, pe ferestre mari. În tabelul de direcții al organizatorilor, întrebarea intră la „Hardware/resources", iar fixarea și RQAOA la „QAOA/search". O a doua întrebare, tot măsurabilă, este dacă cele două codificări sunt la fel de exacte când crește rata de erori din citiri.

## Despre matrice

Aceeași fereastră de genom, adică aceleași SNP-uri, dă matrice diferite la tehnologii diferite: coloanele sunt aceleași, dar rândurile, adică citirile, diferă. Illumina dă multe rânduri, fiecare cu puține cifre completate, iar Nanopore dă puține rânduri, aproape pline. Pentru o tehnologie dată, ambele codificări pornesc din aceeași matrice, deci comparația dintre codificări se face corect, pe aceleași date. De aici tabelul de două ori două al experimentului: ambele codificări pe matricea Illumina și ambele pe matricea Nanopore.

O precizare de onestitate: profilurile sintetice sunt parametrizate după valori publicate: Illumina are citiri de aproximativ 150 de perechi de baze și o eroare de aproximativ 0,1%, Nanopore are citiri de 1–100 kb și o eroare de aproximativ 5%, iar în genomul uman apare un SNP heterozigot la 1.000–2.000 de baze. O citire Illumina reală acoperă deci de obicei 0–1 SNP, iar profilul Illumina cu 2–3 SNP-uri pe citire este deliberat optimist, ca în notebook-ul organizatorilor. Spunem asta când prezentăm profilurile.

## Cele două codificări

În matrice, 0 este litera de referință, 1 este litera alternativă, iar lipsa înseamnă că citirea nu acoperă poziția. Cele două codificări au funcții de pondere separate, pentru că se uită în direcții diferite. Amândouă produc un model Ising, adică un Hamiltonian H = c0 + Σ h·Z + Σ J·Z·Z, în care J negativ favorizează biți egali, iar J pozitiv favorizează biți diferiți.

Codificarea pe poziții este cea a organizatorilor. Pentru fiecare citire, funcția se uită în lungul rândului: pentru fiecare pereche de poziții pe care citirea o acoperă, dacă valorile sunt egale crește contorul A al perechii („la fel"), altfel crește contorul D („opus"). Cuplajul perechii este J = (D − A) / 2, iar offsetul este jumătate din numărul total de voturi, Σ (A + D) / 2. Un qubit reprezintă o poziție, iar valoarea lui este litera pe care o are copia A la acea poziție. Șirul măsurat este deci chiar copia A, iar copia B este complementul lui.

Codificarea pe citiri este a noastră. Funcția se uită în josul coloanelor: pentru fiecare pereche de citiri, parcurge coloanele acoperite de ambele și, dacă valorile sunt egale, crește contorul de potriviri a, altfel crește contorul de contradicții d. Cuplajul perechii de citiri este J = (d − a) / 2, iar offsetul este Σ (a + d) / 2. Un qubit reprezintă o citire, iar valoarea lui spune în ce echipă (copie) este pusă citirea. Citirile din aceeași copie se potrivesc, cele din copii diferite se contrazic pe fiecare poziție comună, deci un J pozitiv le separă în echipe diferite. Copiile se reconstruiesc la final prin majoritate pe fiecare coloană, după ce echipa B este întoarsă.

Folosim `/ 2` la ambele ca energia, după adăugarea offsetului, să fie exact numărul de voturi încălcate, ca în notebook-ul organizatorilor, iar energiile celor două codificări să fie comparabile. Împărțirea la doi nu schimbă soluția, doar scalează energiile; forma fără `/ 2` este și ea corectă.

Cele două funcții sunt oglinzi una alteia. Dacă transpunem matricea, funcția pe poziții aplicată pe matricea transpusă numără exact ce numără funcția pe citiri aplicată pe matricea originală, fiindcă fiecare „poziție" din matricea transpusă este de fapt o citire. Păstrăm totuși două funcții separate, pentru claritate, și folosim această proprietate ca test: contoarele celor două trebuie să coincidă pe matricea transpusă.

Cele două probleme nu sunt însă echivalente în general. Variabilele sunt diferite (poziții față de citiri), iar funcțiile de cost sunt diferite. Pe date curate dau aceleași copii, dar cu erori pot alege soluții diferite. De exemplu, o citire lungă cu o eroare produce multe voturi greșite pe perechile de poziții, pe când în codificarea pe citiri eroarea se diluează în ponderea dintre citiri. Acesta este motivul pentru care a doua întrebare a investigației are sens.

## Fereastra și downsampling-ul

Fereastra se alege peste coloane, adică peste SNP-uri, și o decidem noi. Citirile nu se aleg: fiecare fereastră ia toate citirile care acoperă cel puțin două dintre SNP-urile ei, păstrând doar coloanele ferestrei. Singurul control pe rânduri este downsampling-ul, adică păstrarea aleatoare, cu seed fix, a unui număr țintă de citiri.

Funcția de downsampling este aceeași pentru ambele tehnologii, dar o aplicăm diferit. La Nanopore, cu codificarea pe citiri, păstrăm 10–15 citiri, pentru că fiecare citire lungă acoperă oricum toată fereastra și altfel numărul de qubiți ar depăși numărul de poziții. La Illumina, cu codificarea pe poziții, nu facem downsampling: folosim toate citirile, iar numărul de qubiți este numărul de SNP-uri. La Illumina, cu codificarea pe citiri, nu facem downsampling, pentru că ar rămâne multe poziții nelegate; scopul acolo este să arătăm câți qubiți ar fi necesari, adică sute. Zidurile sunt deci în oglindă: pe poziții fereastra nu mai încape când are multe SNP-uri, iar pe citiri nu mai încape când are multe citiri.

Limitarea numărului de SNP-uri prin fereastră nu este gratuită. Pentru o regiune mare sunt necesare mai multe ferestre, iar ele trebuie lipite pe SNP-urile comune. Ferestrele suprapuse sunt soluția clasică pentru metoda pe poziții, nu apar în materialele organizatorilor și sunt opționale. Se taie cu pas fix, de exemplu SNP-urile 1–10, apoi 8–17, apoi 15–24, cu suprapunere de două-trei SNP-uri. Fiecare fereastră se rezolvă separat, apoi se lipesc pe SNP-urile comune, întorcând o fereastră dacă se potrivește doar în oglindă. Suprapunerea trebuie să fie legată de citiri; dacă nu este, acolo este o graniță naturală, impusă de date, și declarăm un bloc separat.

## Contractele dintre module

Contractele sunt obiecte simple de date (dataclass-uri în Python) pe care și le transmit modulele, ca niște fișe cu câmpuri fixe. Se scriu primele, înainte de orice altceva, ca fiecare coleg să poată lucra pe modulul lui fără să-l aștepte pe altul. Se află în pachetul `core`.

`SourceData` conține tot ce aduce o sursă dintr-o regiune mare: pozițiile SNP-urilor heterozigote, alelele de referință și alternative, citirile cu începutul lor și valorile 0, 1 sau lipsă pe fiecare SNP, opțional literele brute, etichetele speciale (de exemplu coloana V1 și coloana V2) și, dacă se cunoaște, perechea adevărată de haplotipuri.

`SelectionParams` descrie o tăietură: fereastra (primul SNP și câte SNP-uri), numărul țintă de citiri pentru downsampling, dacă se aplică, și seed-ul.

`ReadMatrix` este matricea citiri × SNP-uri cu valori 0, 1 și −1 pentru lipsă, identificatorii citirilor și ai SNP-urilor, tehnologia, etichetele speciale și adevărul, dacă există. Este intrarea comună pentru ambele codificări și pentru baseline.

`IsingProblem` conține cuplajele J (triunghiul de sus), câmpurile locale h (zero la început, apar după fixare), offsetul, ce reprezintă un qubit („poziție" sau „citire"), lista variabilelor fixate și o funcție de decodare care transformă un șir de biți în haplotipuri. Energia unui șir se calculează aici, ca offset plus Σ h·z plus Σ J·z·z.

`Solution` conține șirul ales, energia lui, distribuția completă de eșantioane, dacă metoda este cuantică, resursele folosite și metadatele rulării.

`Haplotypes` conține copia A, copia B și atribuirea fiecărei citiri la o copie sau marcarea ei ca suspectă.

`ResourceReport` conține numărul de qubiți, porțile cu doi qubiți înainte și după compilare, adâncimea, numele backend-ului și, pentru hardware, job ID-ul.

`ResultRow` este un rând din results.json: instanța, metoda, codificarea, transformările, metricile și resursele.

## Modulele, pe rând

### sources

`SyntheticSource` generează o pereche ascunsă de haplotipuri pe o regiune și citiri din ea, după un profil de tehnologie. Profilul spune câte SNP-uri acoperă în medie o citire, acoperirea și rata de erori. Avem un profil scurt, de tip Illumina, și unul lung, de tip Nanopore, la care o citire acoperă toată fereastra. Generatorul poate produce și litere A/C/G/T, iar un pas separat le traduce în 0 și 1 (0 este litera de referință, 1 litera alternativă), ca să arătăm drumul complet de la litere la matrice. Tot aici se construiesc cei doi pacienți pentru demonstrația cis/trans, prin alegerea copiei pe care stă fiecare dintre cele două mutații.

### selection

`WindowSelector` taie o fereastră de SNP-uri din `SourceData` și păstrează doar citirile care acoperă cel puțin două SNP-uri din fereastră, cu coloanele restrânse la fereastră. `Downsampler` păstrează aleator, cu seed fix, un număr țintă de citiri. Selecția este separată de sursă pentru că fereastra și numărul de citiri sunt butoanele principale ale experimentelor: le schimbăm fără să atingem restul codului.

### encodings

O clasă de bază `Encoding` cu metoda `encode`, care primește `ReadMatrix` și întoarce `IsingProblem`, și două implementări cu funcții de pondere separate. `PositionEncoding` este metoda organizatorilor. `ReadEncoding` este metoda noastră. Formulele și decodarea sunt descrise în secțiunea despre cele două codificări.

### transforms

`fix_variable` fixează un qubit (prima poziție la 0, sau prima citire în copia A). Elimină oglinda și scade un qubit, iar cuplajele cu qubitul fixat devin câmpuri locale h, care în circuit sunt porți rz pe un singur qubit. `prune` elimină cuplajele cu |J| sub un prag (opțional). `WindowStitcher` taie o fereastră mare în ferestre suprapuse, le rezolvă separat și le lipește prin majoritate pe SNP-urile comune (opțional, doar pentru cazul cu 50 de SNP-uri pe poziții).

### solvers

Două familii. Solverele pe Ising au metoda `solve`, care primește `IsingProblem` și întoarce `Solution`. `ExactSolver` enumeră toate șirurile, până la aproximativ 22 de qubiți, și este referința. `QAOASolver` construiește circuitul: H pe toți qubiții, apoi pe fiecare strat rz pentru câmpurile h, rzz pentru fiecare J nenul și rx pe fiecare qubit, cu unghiurile scalate la maximul valorilor absolute. Îl antrenează pe statevector cu COBYLA, din două puncte de start cu seed fix, și eșantionează 2048 de măsurători prin modulul de execuție. `RQAOASolver` (opțional) rulează QAOA, fixează perechea cea mai sigură, micșorează problema și repetă până o poate rezolva exact. Baseline-ul clasic `GreedyStitching` lucrează direct pe `ReadMatrix`: pornește de la o citire și lipește ce se potrivește. Există ca să arătăm unde cade metoda intuitivă la erori și la ordinea citirilor.

### execution

O clasă de bază `Executor` cu metoda `run`, care primește circuitul și numărul de măsurători și întoarce măsurătorile și un `ResourceReport`, cu trei implementări. `IdealSimulator` folosește Aer, fără zgomot. `FakeHardware` folosește un backend fals din Qiskit, cu harta cipului și zgomotul lui. `IBMHardware` folosește Runtime SamplerV2, cu cheia doar din variabilă de mediu, niciodată în cod, în notebook sau în output, și cu înregistrarea backend-ului, a numărului de măsurători și a job ID-ului într-un fișier `jobs.json` înainte de a aștepta rezultatul. Toate au și un mod „doar compilare", care transpilează circuitul pentru un cip și raportează resursele fără să ruleze nimic; așa obținem resursele pentru ferestrele mari. Atenție la ordinea biților din Qiskit: qubitul 0 este în dreapta, iar noi citim poziția 1 în stânga. Conversia se face o singură dată, aici.

### evaluation

Calculează pentru fiecare soluție costul pe perechi, costul MEC (câte litere trebuie corectate ca fiecare citire să se potrivească cu copia ei cea mai apropiată), acuratețea față de adevăr ținând cont de oglindă, numărul de switch errors, frecvența optimului și costul mediu din eșantioane, verdictul cis/trans pentru coloanele V1 și V2 și dacă este corect, plus resursele din `ResourceReport`.

### experiments

Definește instanțele și metodele ca simple configurații în Python și le rulează pe toate, cu seed-uri fixe, scriind rânduri `ResultRow` în `results.json`. Schema fișierului se fixează devreme, ca interfața să poată fi construită în paralel pe un fișier de probă.

### ui

Interfața este o pagină web statică (HTML și JavaScript, cu animații în D3) care citește doar `results.json`. QAOA nu rulează în browser, deci pentru metodele cuantice interfața arată combinațiile calculate dinainte. Tot ce este simplă numărătoare, de exemplu câți qubiți ar trebui pe poziții și câți pe citiri, se calculează în timp real.

## Experimentele

Experimentul principal compară cele două codificări pe ferestre de 10, 20 și 50 de SNP-uri, citite atât cu profilul Illumina, cât și cu profilul Nanopore. Pentru fiecare combinație măsurăm qubiții, porțile cu doi qubiți după compilare și adâncimea; acestea se pot obține pe tot intervalul, pentru că doar compilăm circuitul. Calitatea, adică dacă se găsesc copiile corecte, o măsurăm prin simulare doar unde încape, până la aproximativ 20 de qubiți. Fereastra de 10 este punctul de plecare, fereastra de 20 este comparația principală, iar fereastra de 50 arată zidul pentru codificarea pe poziții.

Fixarea qubitului se compară cu varianta fără fixare, pe aceleași instanțe: qubiți, porți, adâncime și frecvența optimului. Cu fixare, histograma are un singur vârf în loc de două.

Experimentul cu numărul de citiri păstrează aceeași fereastră și rulează cu 6, 10 și 15 citiri păstrate, ca să arate de la ce număr în jos metoda începe să greșească. Răspunde la întrebarea de ce folosim 10–15 citiri.

Experimentul despre erori crește rata de erori din citiri și compară acuratețea celor două codificări, ca să vedem dacă diferă (a doua întrebare a investigației).

Demonstrația cis/trans generează doi pacienți cu aceleași genotipuri, dar cu mutațiile V1 și V2 așezate diferit, pe aceeași copie sau pe copii diferite, și verifică dacă ambele codificări dau verdictul corect. Cu citiri scurte, legătura dintre V1 și V2 trece printr-un lanț de SNP-uri intermediare, iar o singură verigă greșită inversează verdictul; cu citiri lungi, o citire acoperă direct ambele mutații. Este ilustrativ, nu diagnostic.

Baseline-ul de lipire se rulează pe aceleași date, inclusiv cu citirile amestecate de mai multe ori, ca să se vadă că rezultatul lui depinde de ordine, iar cel al metodelor cu voturi nu.

Rularea pe hardware real se face pe fereastra de 10 în ambele codificări și, dacă se poate, pe fereastra de 20 în codificarea pe citiri, cu aceleași date, același număr de straturi, aceleași măsurători, aceleași opțiuni de mitigare și aceleași unghiuri antrenate. Raportăm backend-ul, numărul de măsurători, porțile după compilare, adâncimea, frecvența optimului și job ID-ul.

Opțional, dacă rămâne timp: RQAOA, pruning, lipirea pe ferestre suprapuse pentru 50 de SNP-uri pe poziții și mitigarea pe hardware (dynamical decoupling, twirling, corecția erorilor de citire).

## Ordinea de construcție

Faza 0 este scheletul: repo-ul, un fișier `CLAUDE.md` cu arhitectura și regulile de cod, pachetul `core` cu contractele și schema `results.json` cu un fișier de probă. Faza 1 este drumul complet pe exemplul din notebook: sursa sintetică, selecția, `PositionEncoding`, `ExactSolver`, `QAOASolver` pe simulator ideal și evaluarea de bază. Criteriul de reușită este că pipeline-ul reproduce notebook-ul: pe cele 5 SNP-uri și 9 citiri găsește `01101` cu cost 1. În paralel, interfața începe pe fișierul de probă. Faza 2 adaugă `ReadEncoding`, `GreedyStitching`, fixarea qubitului, profilurile Illumina și Nanopore, modul „doar compilare" și rularea experimentului principal pe date sintetice. Criteriul de reușită este că ambele codificări dau aceleași copii pe date curate. Faza 3 aduce rularea pe hardware real, demonstrația cis/trans și experimentul cu numărul de citiri. Faza 4 este opțională: RQAOA, pruning, ferestrele suprapuse și mitigarea.

## Împărțirea pe oameni

Propunere pentru cinci oameni. Unul se ocupă de prezentare de la început, pentru că povestea și diagramele sunt deja stabilite. Lavinia face faza 0, apoi `core`, `encodings` și `transforms`, și verifică fiecare merge. Al doilea coleg face `sources` și `selection`. Al treilea face `solvers` și `execution`, inclusiv contul IBM și hardware-ul. Al patrulea face `evaluation`, `experiments` și interfața, pornind de la fișierul de probă. Lucrăm pe branch-uri, câte unul pe modul, cu câte un test mic pe modul.

## Limite și direcții viitoare

În acest hackathon lucrăm doar pe date sintetice. Datele reale Genome in a Bottle (proba HG002, Illumina și Nanopore pe aceeași regiune), extrase cu pysam într-un script separat în Colab și citite de o sursă `GiabSource`, rămân o limită declarată și o direcție viitoare.

## Reguli de cod

Fiecare modul face un singur lucru și comunică doar prin contractele din `core`. Fără variabile globale, fără logică de experiment în solvere și fără Qiskit în afara modulelor `solvers` și `execution`. Cod curat, modular și fără structuri încâlcite, cu clase de bază acolo unde există mai multe implementări (`Encoding`, `Executor`, solvere). Toate rulările au seed fix. Cheia IBM nu apare niciodată în cod, în notebook sau în output. Funcțiile din notebook-ul organizatorilor pe care le refolosim sunt marcate ca atare, pentru declarația de la predare.
