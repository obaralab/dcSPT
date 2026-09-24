# dcSPT — dual-colour single-particle tracking

Two colours, tracked independently, then compared. A separate toolkit from **SPTinMatlab**, not a
mode of it.

## Why separate

SPTinMatlab is built around one tracked species and its relationship to segmented organelles: a cell
is one movie, a channel is a mask, and the analyses ask "is this molecule near the ER". Dual colour
asks a different question with a different shape — two moving species, two clocks, and a partner that
may not have been imaged when you were looking. Bolting that onto the single-colour pipeline would
have meant a channel token threaded through every file name, a second time base through every
consumer, and a "which colour is this" question at every call site. The pair is the unit here
instead.

## What was taken from SPTinMatlab, and what was left

Copied and renamed `spt_* -> dc_*`, so both toolkits can sit on one MATLAB path:

| here | from | why |
|---|---|---|
| `core/dc_dog.m` | `spt_dog` | Difference-of-Gaussians spot filter, separable and bit-identical to the reference |
| `core/dc_detect.m` | `spt_detect` | local maxima, sub-pixel centroid, curvature and width gates |
| `core/dc_ridge.m` | `spt_ridge` | the curvature-ratio gate |
| `core/dc_measure.m` | `spt_measure` | per-spot intensity |
| `core/dc_tiff_pages.m` | `spt_tiff_pages` | one open handle per stack (3.7x faster than `imread` per page) |
| `core/dc_tiff_calib.m` | `spt_tiff_calib` | pixel size and frame interval from the file's own metadata |

Written here, because the single-colour pipeline had no need of it:

| | what it does |
|---|---|
| `core/dc_tiff_labels.m` | reads ImageJ's per-page slice labels — which channel and which timepoint the acquisition says each page is |

Deliberately **not** taken:

* **the ER link modes.** SPTinMatlab can bias or forbid a link by how much of it lies off a
  segmented ER. That is a statement about that experiment's biology and rests on a mask this toolkit
  does not have. `dc_track` is the plain LAP that remains.
* **the organelle-skeleton bleedthrough gate.** It rejects detections lying along a segmented
  filament. Here the thing bleeding through is the *other particle channel*, which is spot-shaped —
  a different problem needing a different answer, not a borrowed one.

## The model

A **colour** (`dc_channels`) says which pages of the stack it occupies and which timepoint of the
acquisition each of its frames belongs to. Everything in that description is an integer:

```
key  pages        tp           dt_s       nFrames
c2   1 3 5 7 ...  1 2 3 4 ...  0.0534237  2000     the odd pages
c4   2 4 6 8 ...  1 2 3 4 ...  0.0534237  2000     the even ones — the same timepoints, one page later
```

and it is **read from the file**, not assumed. ImageJ keeps a label per page and a hyperstack
acquisition writes its own indexing into it — `c:2/4 t:1/2000 - cellname #1` — so `dc_tiff_labels`
can be asked instead of deducing "odd pages are colour A" from a stride. That deduction assumes the
acquisition never dropped a page, never started on the other channel and never changed order, and
when it breaks the pages still alternate plausibly for a while: half of one colour's frames are handed
to the other, silently. A dropped page is in the test suite for exactly this reason.

**Frame numbers, not seconds, everywhere.** A second here would be a frame number multiplied by a
nominal interval: these files carry one `finterval` for the whole stack and **no per-page timestamps at
all**, so every second is derived and every derivation accumulates. The 1.1% between the `0.054`
typed into the first attempt's calibration and the `0.0534237` its builds carried is half a second
across a 1,791-frame track. A timepoint is what the microscope recorded — exact, and also the unit a
tracking gap is naturally counted in: frame k+2 after frame k is a gap of one, with nothing to compare
against a tolerance.

So the two colours are related by **set intersection on their timepoints**, not by matching times
within a tolerance. `c:2 t:5` and `c:4 t:5` are one moment imaged twice, and `dc_align` reports the
**page gap** between those two exposures rather than calling them simultaneous. A colour with no frame
at a timepoint contributes **nothing** to the merged view and is named — treating it as a distant
partner would turn "the other colour was not imaged" into "it was not there", which is the opposite
conclusion and invisible once it reaches a histogram.

`dt_s` is still per colour, because a diffusion coefficient has to come out in µm²/s eventually. It is
applied **once, at the end, to a frame count** and never indexed by — which is why a colour whose
interval is unknown still detects and tracks.

## The data model in one paragraph

Per-channel arrays stay per channel — the two colours do not share a frame axis, so `[frames x
tracks]` cannot hold both without resampling or lying about what a row means. The **merged** view is a
flat spot list where the shared moment is a column rather than an index: every spot carries a **channel
flag**, its own **frame**, and the **acquisition timepoint and page** it came from, all integers the
microscope supplied, checked for exact equality against the registry that holds each colour's maps.
Track ids are renumbered as a colour is added, because both trackers number from zero and "track 0"
otherwise exists twice. Masks are shared between the colours and looked up **by timepoint**, because
`page = frame + 1` gives two different answers for one moment. `dc_align` says which of four ways the
two colours' timepoints relate and `dc_merge` resolves any query into a timepoint — including the case
that looks fine and is not: the same frame count and the same timepoints, but two consecutive
exposures rather than one instant. Full argument in [docs/DATA_MODEL.md](docs/DATA_MODEL.md).

## Running it

```matlab
addpath('<this folder>/core', '<this folder>/drivers');
stack = '.../cellA.tif';
C   = dc_channels('fromStack', stack, 0.0267);  % ask the file which pages are which colour
cel = struct('stack', stack, 'base','cellA');
prm = struct('diamUm',0.4, 'thrAbs',60, 'pxUm',0.1, 'linkUm',1.5, 'gapUm',1.5, 'maxGap',1);
R   = dc_process_cell(cel, C, prm);             % R(1) and R(2): both colours, tracked independently

D = dc_dataset('new','cellA',C);
D = dc_dataset('addChannel', D, R(1));
D = dc_dataset('addChannel', D, R(2));
disp(dc_align(D).text)                          % how the two colours' timepoints relate
[S, info] = dc_merge(D, struct('tp', 6));       % both colours at timepoint 6
```

Where a stack carries no labels, `dc_channels('interleaved', nPages, dtPage)` assumes alternating
pages instead — an assumption, and it says so.

Tests: `run('tests/run_all.m')`.

## Where this is going

Tracking is done; relating the two colours is not. Next: chromatic **registration** (nothing here
corrects it yet, and a 100-300 nm offset is the same size as the distances being measured), then
partner distances with an explicit time-matching rule, then colocalization and encounter dwell times.
