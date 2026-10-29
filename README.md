# fft — быстрое преобразование Фурье (задание №8)

ДПФ в унитарной нормировке, алгоритм прореживания по частоте, всё in place:

$$
y_k = \frac{1}{\sqrt{N}} \sum_{j=0}^{N-1} x_j e^{\mp \frac{2\pi i}{N} jk},
\qquad \text{«−» — прямое, «+» — обратное.}
$$

Реализация — современный Fortran 2023 (`iso_fortran_env`, `merge`, сюры OpenMP),
двойная точность, размер $N = 2^p$ (нестепень двойки паддится нулями).

## Структура

- `src/precision_mod.f90` — kind-параметр `dp` (`iso_fortran_env: real64`);
- `src/io_mod.f90` — чтение `data.dat` и запись `result.dat`/`abs.dat`;
- `src/fft.f90` — ядро: `make_twiddles` (таблица $W_N^k$ из $N/2$ факторов),
  `fft_transform(x, inverse)` (прямое/обратное, ключ-параметр), внутренние
  `dif_stages` и публичный `bit_reverse_permute`;
- `app/main.f90` — CLI: `fft [<data.dat>] [<result.dat>] [<abs.dat>] [inverse|-i]`;
- `test/` — тесты (`fpm test` или `ctest`; файлы `test_*.f90`);
- `plot/` — demo: известный сигнал → forward → inverse (round trip) → `demo.png`.

## Форматы файлов (по заданию)

- `data.dat`: первая строка `# N`, дальше N пар `re im` (пробелы любые);
- если N не степень двойки — паддинг нулями до ближайшей;
- `result.dat` — тот же формат с заголовком `# N` (N — фактическое,
  т.е. дополненное число точек);
- `abs.dat` — модули `|y_k|` по одному в строке, без заголовка.

## Сборка и запуск

### fpm (автоматический захват src/, app/, test/)

```sh
fpm test                 # все тесты (OpenMP включён по умолчанию)
fpm run                  # data.dat -> result.dat, abs.dat + время
fpm run --profile perf   # с -O3 -march=native -ffast-math
```

OpenMP включён через `features`/`profiles` манифеста (fpm ≥ 0.13).
Fallback для старых fpm: `FPM_FFLAGS="-O3 -fopenmp" FPM_LDFLAGS="-fopenmp" fpm test`.

### CMake

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release   # Debug — с санитайзерами/fpe-трапами
cmake --build build -j
ctest --test-dir build                       # 4 теста
./build/fft [data.dat] [result.dat] [abs.dat] [inverse|-i]
```

`find_package(OpenMP REQUIRED COMPONENTS Fortran)` + `OpenMP::OpenMP_Fortran`.

### Параллелизм

Число потоков — переменная окружения `OMP_NUM_THREADS`:

```sh
OMP_NUM_THREADS=16 ./build/fft data.dat result.dat abs.dat
```

Параллельные участки: генерация поворачивающих множителей, бабочки каждой
стадии (`collapse(2)` по блокам × j — серия делится поровну и при одном
большом блоке, и при $N/2$ маленьких), битреверс-перестановка (пары
$(i, \mathrm{bitrev}(i))$ не пересекаются — обмены безопасны).

## Тесты

| файл | что проверяет |
|---|---|
| `test_make_twiddles.f90` | таблица $W_N^k$: значения, |w|=1, рекуррентность, w(1)=1, w(n/4)=±i |
| `test_fft_dft.f90` | прямое и обратное против наивного $\mathcal{O}(N^2)$ ДПФ задания, n = 2…256, псевдослучайные данные |
| `test_fft_identity.f90` | round-trip, дельта-импульс, константа, Парсеваль, линейность; большие размеры до $2^{20}$ |
| `test_bit_reverse.f90` | перестановка против независимой эталонной реверс-реализации, инволюция |

## Демо (gnuplot)

`plot/demo.sh` — round trip на известном сигнале
`0.7·cos(2π·32j/N) + 0.4·sin(2π·120j/N)`, N = 4096:

1. генерирует `plot/data.dat` (формат задания, мнимая часть = 0);
2. прогоняет программу вперёд: `data.dat → result.dat, abs.dat`;
3. подаёт `result.dat` назад с ключом `inverse`: `result.dat → round.dat`;
4. рисует `plot/demo.png`: сигнал, спектр в лог-шкале (пики ровно на k = 32 и
   k = 120, |X| = 0.7·√N/2 = 22.4 и 0.4·√N/2 = 12.8), реконструкция поверх
   оригинала и погрешность восстановления по всем точкам
   |x_j − reconstructed_j| (~10⁻¹⁴ при масштабе сигнала ~1).

```sh
sh plot/demo.sh          # нужен собранный бинарь (или FFT_BIN=<путь>) и gnuplot
```

## Производительность

Double precision, in place; компилятор gfortran 16.2,
`-O3 -march=native -ffast-math`; задача стриминговая, упор в пропускную
способность памяти:

| N | 1 поток | 16 потоков | ускорение |
|---|---|---|---|
| $2^{20} \approx 10^6$ | 0.125 s | 0.039 s | 3.2× |
| $2^{22} \approx 4\cdot10^6$ | 0.837 s | 0.225 s | 3.7× |
| $2^{24} \approx 1.7\cdot10^7$ | 3.32 s | 0.98 s | 3.4× |

(время без чтения/записи файлов; text-I/O остаётся основной накладкой для
$x \sim 10^7$ — десятки секунд на запись `result.dat`/`abs.dat`).
