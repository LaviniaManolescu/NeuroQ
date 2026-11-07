# Formatul JSON pentru surse

O sursă (`SourceData`) se poate salva și citi ca JSON rar: fiecare citire listează doar SNP-urile pe care le observă. Codul este în `neuroq/sources/json_io.py` (`save_source_json`, `load_source_json`, `save_reads_data_json`), într-un singur modul, ca formatul să se poată schimba ușor. **SNP-urile sunt numerotate de la 1** în fișier, ca în interfață și în documentație; în cod devin indici de la 0 o singură dată, la citire.

## Formatul NeuroQ (`format_version: 1`)

```json
{
  "format_version": 1,
  "n_snps": 5,
  "technology": "illumina",
  "snp_positions": [1000, 2000, 3000, 4000, 5000],
  "ref_alleles": ["A", "A", "A", "A", "A"],
  "alt_alleles": ["G", "G", "G", "G", "G"],
  "labels": {"V1": 2, "V2": 4},
  "reads": [
    {"id": "r1", "start": 1000, "snps": {"1": 0, "2": 1, "3": 1}},
    {"id": "r2", "snps": {"3": 1, "5": 0}}
  ],
  "truth": {"copy_a": "01101", "copy_b": "10010"},
  "meta": {"profile": "illumina", "seed": 7}
}
```

| câmp | obligatoriu | ce înseamnă |
|---|---|---|
| `format_version` | da | acum `1` |
| `n_snps` | da | numărul de SNP-uri (coloane) |
| `reads` | da, cel puțin una | fiecare citire: `snps` = `{indice SNP de la 1: 0 sau 1}`; ce lipsește înseamnă neobservat (−1 în cod). `id` și `start` sunt opționale |
| `technology` | nu | implicit `"unknown"` |
| `snp_positions` | nu | coordonatele genomice; implicit 1000·k pentru SNP-ul k |
| `ref_alleles`, `alt_alleles` | nu | implicit `"A"` și `"G"` |
| `labels` | nu | coloane speciale (V1, V2), tot de la 1 |
| `truth` | nu, doar la date sintetice | perechea adevărată, câte un caracter pe SNP |
| `meta` | nu | informații libere (profil, seed); nu intră în `SourceData` |

Implicit, `start` este poziția primului SNP observat. Literele brute (`raw_letters`) nu se salvează.

## Formatul interfeței (`reads_data.json`)

`load_source_json` citește și fișierul produs de interfață, iar `save_reads_data_json` îl scrie:

```json
{
  "metadata": {"num_snps": 9, "num_reads": 17, "description": "..."},
  "reads": [{"1": 1, "2": 0, "3": 0}, {"2": 0, "3": 0, "4": 1}]
}
```

Acest format păstrează doar citirile: fără identificatori, poziții sau adevăr.

## Erori

Un fișier gol, un JSON invalid, un indice SNP în afara lui 1..`n_snps`, o valoare diferită de 0/1, un fișier fără citiri sau o versiune necunoscută dau `ValueError` cu un mesaj care spune fișierul și citirea.
