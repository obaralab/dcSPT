# The data model

Decided before the tools, because every tool is a consequence of it.

## The one fact everything follows from

**The two colours do not share a frame axis.** In the real data here `ch24` runs at 26.7 ms and `ch3`
at 53.4 ms; interleaved colours are offset by half a frame; a strobed colour has its own rate again.
So "frame 7" is a different instant in each colour, and a `[frames x tracks]` matrix — the shape all
the per-track maths wants — cannot hold both. Forcing it means either resampling one colour onto the
other's grid, which invents localizations that were never measured, or padding with NaN, after which
row 7 silently means two different times.

Therefore:

* **per-channel arrays stay per-channel.** Tracking, MSD, D, step vectors — all of it is computed on
  one colour's own grid and never merged.
* **the merged view is a flat list of spots, where time is a COLUMN and not an index.** That
  representation has no frame-axis problem, which is exactly why it is the one that can be merged.

One is canonical, the other derived from it. Neither is optional.

## The structure

```
DCS                      one cell, both colours
  .cell                  char, the cell key — shared; the colours are of the SAME field of view
  .channels  1xC struct  THE REGISTRY — the only place a clock lives
      .key     'a'       short id; this is the flag value carried everywhere else
      .label   'Halo-Sec61B'
      .dt_s    0.0267118 this colour's frame interval
      .t0_s    0         when its frame 0 happened, on the cell's clock
      .stride/.offset    which pages of the stack it is
      .stack             the file it came from
      .pxUm              per channel: two cameras need not be binned alike
      .reg               registration onto the cell frame, and the residual that justifies it
  .tracks    1xC struct  canonical, per channel, on its own frame axis
      .matrix  [m x n x 3]  frame, x, y
      .lengths .trackIds .steps .msd ...
  .spots     table       DERIVED, both colours in one list
      .ch      categorical  the flag
      .frame   within that colour
      .t_s     the CELL's clock — what makes a merge mean anything
      .x .y .q .trackId .spotId
  .masks                 shared by both colours: they are properties of the cell
      .er .mito          stacks + their own page rate, looked up BY TIME
```

## Why the flag is a key and not a 1 or a 2

`ch` is `'a'`/`'b'`, stored as a categorical whose categories are the registry's keys. A bare 1/2 is
an index into something the file does not carry: subset a table, drop a colour, add a third, and the
integers still parse but mean something else. A key is self-describing, survives subsetting, and
fails loudly rather than quietly when it does not match the registry.

## Why dt is never on a spot

You said dt belongs "somewhere else", and the reason is worth stating: on the spot, two spots of the
same colour could disagree and nothing would notice. In the registry there is exactly one value per
colour, and `dc_validate` asserts that every spot's `t_s` equals `t0_s + frame*dt_s` for its own
channel. The spot carries `t_s` so that drawing and merging need no lookup — a denormalization, and
it is only safe because one function builds it and a test pins the identity.

The first attempt at dual colour shows the failure mode: `cs_calib.mat` carried a second, typed frame
interval of `0.054` while the builds carried `0.0534237`. A second copy of a number the data already
holds drifts by 1.1%, which is half a second across a 48-second track.

## Why masks are looked up by TIME

ER and mito belong to the cell, not to a colour — one segmentation serves both. But a mask page is
tied to a PAGE of the acquisition, and each colour samples pages differently, so `page = frame + 1`
gives two different answers for the same instant. The only rule that is right for both is
`page = f(t_s)`. (This is also the bug the single-colour pipeline still has on de-interleaved runs,
where every viewer is wrong by the stride.)

## What the merge panel must say

"Same frame number or not" is really "how do the two clocks relate", and there are four answers.
`dc_align` returns which one, with the numbers, and the panel says it in words:

| relation | what it means | what the panel does |
|---|---|---|
| `simultaneous` | same dt, same t0 | frames correspond 1:1; step by frame |
| `interleaved` | same dt, t0 differs by less than one frame | same COUNT but never simultaneous — step by frame, and show the constant lag |
| `subsampled` | one dt is ~k x the other | only 1 in k frames has a counterpart; step in the fast colour and show the partner's age |
| `unequal` | anything else | step in seconds; nearest partner within a stated tolerance |

The rule underneath: **the merged view is driven by time, never by frame index.** Pairing frame k
with frame k is right only in the first row of that table, and wrong in a way that looks fine in
every other.

Where a colour has no frame near the displayed instant, the panel says so rather than drawing its
last known position as if it were current — the same "unmeasured is not far" rule the time matcher
enforces.
