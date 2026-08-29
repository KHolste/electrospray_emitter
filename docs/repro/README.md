# Reproduktionsfälle

## `avx_stack_alignment.cpp`

Minimaler Reproduktionsfall für den in
[../01_gap_analysis.md](../01_gap_analysis.md) Abschnitt 1.2 beschriebenen
Codegenerierungsfehler.

```sh
# stürzt in etwa der Hälfte der Läufe ab (MinGW-w64 GCC 16.1, Windows):
g++ -std=c++20 -O2 -mavx2 -Iinclude \
    src/elliptic.cpp src/linalg.cpp src/geometry.cpp src/bem.cpp src/meniscus.cpp \
    docs/repro/avx_stack_alignment.cpp -o repro.exe
for i in $(seq 30); do ./repro.exe >/dev/null 2>&1 || echo SIGSEGV; done

# läuft sauber:
g++ ... -mavx2 -mprefer-vector-width=128 ...
g++ ... -msse4.2 ...
```

Ursache: `main()` enthält alignment-pflichtige 256-Bit-Stackzugriffe
(`vmovdqa %ymm0,0x20(%rsp)`), ohne den Stack auf 32 Byte nachzurichten. Die
Win64-ABI garantiert nur 16 Byte. Ob es kracht, entscheidet die
Adressraum-Randomisierung.

Nicht reproduzierbar unter Linux/GCC 13.3. ASan und UBSan melden nichts.
`-mstackrealign` hilft nicht.

**Offen:** ein Reproduktionsfall ohne Abhängigkeit von der Projektbibliothek.
