# metaphaser: установка и запуск

[English documentation](README.md)

metaphaser получает сборку в FASTA и исходные нанопоровые риды в modBAM.
Он находит интервалы, классифицирует риды по стыкам L/R/P/B и сравнивает
модификации между состояниями в одинаковых позициях ДНК.
5mC, 5hmC и 6mA анализируются раздельно.

## Установка

Нужны Linux или macOS, Git и Conda/Mamba. Для Windows используйте Linux,
например через WSL2.

```bash
git clone --branch v0.2.0rc2 https://github.com/andreyeremin/metaphaser.git
cd metaphaser
conda env create -f environment.yml
conda activate metaphaser
python -m pip install -r requirements-tested.txt
python -m pip install --no-build-isolation --no-deps .
metaphaser --version
metaphaser doctor
```

Команды создают отдельное окружение с Python 3.12 и minimap2 2.30.
Последняя команда проверяет установку и выводит версии зависимостей.
Не нужно устанавливать пакет с похожим названием из PyPI. Исходный код
этой программы находится в репозитории выше.

Для сервера без доступа к PyPI предусмотрена
[установка из заранее скачанных wheel-файлов](docs/TROUBLESHOOTING.md#install-on-a-linux-server-without-internet-access).
Python и minimap2 должны быть установлены на сервере заранее.

## Быстрая проверка на примере

```bash
metaphaser demo --outdir demo_run --threads 2
```

Программа создаст 56 искусственных ридов и выполнит весь анализ.
Ожидается один интервал со всеми четырьмя типами стыков.
Отчёт: `demo_run/results/report.html`.

## Подготовка своих данных

Нужны два файла:

- сборка соответствующего образца, например `assembly.fasta`
- исходный BAM с полными последовательностями ридов и тегами MM/ML,
  например `reads.mod.bam`

Предварительное выравнивание BAM не обязательно. FASTQ не содержит нужных
вызовов модификаций. Для BAM, выровненного на другую сборку, добавьте
`--ignore-original-alignments` при проверке и запуске.

```bash
metaphaser check --assembly assembly.fasta --modbam reads.mod.bam
```

## Запуск одной командой

Для моделей с цитозиновыми вызовами во всех контекстах и вызовами 6mA:

```bash
metaphaser run \
  --assembly assembly.fasta \
  --modbam reads.mod.bam \
  --sample OI_SAMPLE \
  --cytosine-context all \
  --modifications 5mC,5hmC,6mA \
  --threads 4 --validation-workers 2 \
  --outdir results/OI_SAMPLE
```

Для данных SETUN с моделью `5mCG_5hmCG`:

```bash
metaphaser run \
  --assembly assembly.fasta \
  --modbam reads.mod.bam \
  --sample SETUN06 \
  --cytosine-context cpg \
  --modifications 5mC,5hmC \
  --threads 4 --validation-workers 2 \
  --outdir results/SETUN06
```

Подставьте свои пути вместо `assembly.fasta` и `reads.mod.bam`.
Для каждого образца используйте отдельное имя и новую папку результатов.
Число процессов проверки не должно превышать `--threads`.
Все сообщения программы и логи написаны на английском.

Контекст и модификации выбирают по модели basecalling.
`--modifications auto` выбирает доступные поддерживаемые каналы.
Отсутствующий канал не превращается в нулевое значение.
Значение `--cytosine-context` по умолчанию — `cpg`, поэтому для всех
контекстов необходимо явно указать `all`.

## Фоновый запуск и продолжение

Чтобы расчёт продолжался после закрытия SSH:

```bash
nohup metaphaser run \
  --assembly assembly.fasta --modbam reads.mod.bam \
  --sample SETUN06 --cytosine-context cpg --modifications 5mC,5hmC \
  --threads 4 --validation-workers 2 --outdir results/SETUN06 \
  > SETUN06.log 2>&1 < /dev/null &
echo $! > SETUN06.pid
tail -f SETUN06.log
```

Ctrl+C останавливает просмотр лога, если сейчас выполняется `tail`.
Сам фоновый расчёт продолжится. После прерывания расчёта повторите ту же
команду с добавлением `--resume`. Программа проверит входные файлы,
версию кода и параметры, затем использует завершённые этапы.
Обновлять программу посреди расчёта нельзя.

## Где смотреть результаты

| Файл | Содержание |
|---|---|
| `report.html` | Общий отчёт, открывается в браузере |
| `elements.tsv` | Отобранные интервалы и числа L/R/P/B |
| `read_classifications.tsv.gz` | Классификация каждого исходного рида |
| `read_locus_classifications.tsv.gz` | Назначения ридов внутри каждого интервала |
| `junction_calls.tsv.gz` | Стыки и признаки неоднозначности или сложности рида |
| `modifications/comparisons.tsv` | Средние уровни, разности и p-value |
| `modifications/read_profiles.tsv.gz` | Значения по отдельным ридам |
| `modifications/panels.tsv.gz` | Одинаковые позиции ДНК, использованные в сравнении |
| `modifications/epihaplotypes.tsv.gz` | Частоты профилей модификаций по состояниям |
| `manifest.json`, `summary.json` | Статус, параметры и итоговые числа |

L пересекает начало интервала, R — его конец по возрастанию координат
сборки. P пересекает соединение конца с началом, B — соединение внешних
фланков. Эти метки не обозначают ориентацию рида. Один рид может пересечь
и L, и R. В колонках `L_simple/R_simple/P_simple/B_simple` подсчитаны
простые однозначно назначенные риды после окончательных фильтров.

В `elements.tsv` попадают интервалы с минимум тремя P-ридами и минимум
тремя L- или R-ридами. Перекрывающиеся интервалы сохраняются отдельно.
Координаты отсчитываются от нуля, правая граница не включается.

В `comparisons.tsv` строка `strand=combined` содержит общий тест с
перестановками внутри ориентаций. Строки `+` и `-` содержат описательные
значения. Разность `delta_pp` выражена в процентных пунктах и считается
как P минус L или P минус R, согласно колонке `contrast`. Только статус
`tested` означает, что p-value рассчитано.

## Общая поправка BH

После завершения всех образцов передайте таблицы всей запланированной
серии. Команда сама выберет подходящие P–L/P–R сравнения.

```bash
metaphaser adjust \
  --comparisons results/OI_SAMPLE/modifications/comparisons.tsv \
                results/SETUN06/modifications/comparisons.tsv \
  --output results/primary_BH.tsv
```

В результате появятся `q_BH` и размер семейства тестов `BH_family_size`.
Список входных образцов нельзя выбирать по полученным p-value.

Подробности метода и ограничения: [README](README.md),
[проверки версии](docs/VALIDATION.md), [устранение проблем](docs/TROUBLESHOOTING.md).
