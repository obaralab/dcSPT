# What the first attempt shows

`IntegratedPipeline/DualColor/analysis_dualcolor/`, July 2025. Read-only; nothing there was changed.

## What it did

Ran the whole single-colour contact-site pipeline **twice** — once per colour, `ch24` and `ch3` —
and then tried to reconcile the two at the level of the **derived contact sites**:

* `cs_reconciliation.mat` pairs sites by centroid distance (≤ 0.30 µm) and IoU (≥ 0.10) on a 48×48
  grid.
* `cs_channel_offset.mat` holds a chromatic offset of **[−0.25, −0.45] µm** (0.515 µm total),
  estimated on one cell at 0.05 µm bins from 42,811 against 7,626 localizations.
* `cs_calib.mat` carries two frame intervals: `dt_s = 0.026712` and `dt2_s = 0.054`.

The result: **10 sites in ch24, 7 in ch3, 1 shared.** Nine A-only, six B-only.

## What went wrong — and what did not

**Not registration.** The obvious reading is that a 0.515 µm offset defeats a 0.30 µm pairing
threshold, and `R.shiftXY = [0 0]` says the offset was never applied. But applying it makes the
pairing *worse* (1 → 0), and a sweep of every shift within ±2 µm in 0.1 µm steps recovers at most
**2 of 10**. If the two colours marked the same structures with a rigid offset, a shift would find
them. The nearest ch3 site to each ch24 site is a **median 2.88 µm** away — a whole different scale
from any chromatic error.

So the sites genuinely do not correspond, and the 0.515 µm offset estimate is itself suspect: it was
fitted between two localization clouds of very different density (5.6× fewer in ch3), and it does not
improve agreement on the thing it was meant to fix.

**The level was wrong.** A contact site is the output of detection → density smoothing → a threshold,
run on each channel separately with 42,811 localizations in one and 7,626 in the other. Pairing those
outputs compounds every difference between the two channels into one yes/no per site, and when it
fails there is no way to see which stage caused it. Comparing **localizations** — each with a time,
a position and an error — keeps the question answerable.

**The two frame intervals disagree with the builds.** `cs_calib.mat` says the second channel is
`0.054` s; the builds say `frameInterval = 0.0534237` for ch3 against `0.0267118` for ch24 — an exact
2:1, which is what interleaved acquisition gives. The rounded 0.054 is 1.1% off, which is ~0.5 s of
drift across a 48 s track. A second, hand-typed copy of a number the data already carries is how that
gets in.

## What dcSPT does differently because of it

1. **Compare localizations, not derived objects.** Partner distance at the localization level, on one
   clock. Sites, if they are wanted, come after.
2. **A registration that cannot be shown to improve agreement is not applied.** Estimating an offset
   is easy; the residual check is the part that matters, and here it would have rejected the estimate.
3. **One dt per colour, from where the tracks came from.** `dc_channels` holds it, `dc_process_cell`
   refuses to run without it, and a declared value that disagrees with the stack's own metadata warns
   rather than being quietly believed.
4. **The density asymmetry is a fact about the experiment, not a nuisance.** 42,811 against 7,626
   localizations in the same cell means the two channels cannot be given the same thresholds and
   called comparable. Anything reported per channel has to carry its own n.
