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

Deliberately **not** taken:

* **the ER link modes.** SPTinMatlab can bias or forbid a link by how much of it lies off a
  segmented ER. That is a statement about that experiment's biology and rests on a mask this toolkit
  does not have. `dc_track` is the plain LAP that remains.
* **the organelle-skeleton bleedthrough gate.** It rejects detections lying along a segmented
  filament. Here the thing bleeding through is the *other particle channel*, which is spot-shaped —
  a different problem needing a different answer, not a borrowed one.

## The model

A **colour** (`dc_channels`) says where its frames are and when they happened:

```
key  stride  offset  dt_s   t0_s   file
a    2       0       0.02   0      ''     odd pages of the shared stack
b    2       1       0.02   0.01   ''     even pages — and it starts one page later
```

Two things this toolkit refuses to assume, both because they are silent when wrong:

* **`dt` is per colour.** Deriving it as page interval x stride assumes the only reason a colour has
  fewer frames is that it shares pages with another. A strobed colour has stride 1 and a dt with
  nothing to do with the page interval. Getting it wrong scales every diffusion coefficient, dwell
  time and rate by that factor.
* **`t0` is carried, not assumed to be zero.** Interleaved colours do not start together. Half a
  page sounds like nothing until it is set beside what is being measured: a molecule at
  0.1 µm²/s covers ~60 nm in 10 ms, which is the size of the colocalization distances this exists to
  report.

`dc_times` gives a colour's **acquisition** clock — when each frame happened, whether or not
anything was detected in it. That is a different question from when a molecule was seen, and both get
asked: `dc_time_match` answers "was the other colour even imaged then" against the first, and "where
was its molecule when mine was seen" against the second. A frame with no partner in time comes back
**unmatched, not far** — treating it as a large distance turns "B was not imaged" into "B was not
there", which is the opposite conclusion and invisible once it reaches a histogram.

## The data model in one paragraph

Per-channel arrays stay per channel — the two colours do not share a frame axis, so `[frames x
tracks]` cannot hold both without resampling or lying about what a row means. The **merged** view is
a flat spot list where time is a column rather than an index: every spot carries a **channel flag**
and a **time on the cell's clock**, and the clock itself lives once per colour in a registry, never
on a spot. Masks are shared between the colours and looked up **by time**, because `page = frame + 1`
gives two different answers for one instant. `dc_align` says which of four ways the two clocks relate
and `dc_merge` resolves any query into an instant — including the case that looks fine and is not:
same frame COUNT, same interval, offset by half a frame, never simultaneous. Full argument in
[docs/DATA_MODEL.md](docs/DATA_MODEL.md).

## Running it

```matlab
addpath('<this folder>/core', '<this folder>/drivers');
C   = dc_channels('interleaved', 0.01);        % two colours, one stack, 10 ms per page
cel = struct('stack','.../cellA.tif', 'base','cellA');
prm = struct('diamUm',0.4, 'thrAbs',60, 'pxUm',0.1, 'linkUm',1.5, 'gapUm',1.5, 'maxGap',1);
R   = dc_process_cell(cel, C, prm);            % R(1) and R(2): both colours, each on its own clock
```

Tests: `run('tests/run_all.m')`.

## Where this is going

Tracking is done; relating the two colours is not. Next: chromatic **registration** (nothing here
corrects it yet, and a 100-300 nm offset is the same size as the distances being measured), then
partner distances with an explicit time-matching rule, then colocalization and encounter dwell times.
