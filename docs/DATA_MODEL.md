# The data model

Decided before the tools, because every tool is a consequence of it.

## The one fact everything follows from

**The two colours do not share a frame axis.** In the real data here one colour of a `ch24` stack runs
at 26.7 ms per page and the other is the page next to it; a subsampled colour has only half the frames;
a strobed colour has its own rate again. So "frame 7" is a different moment in each colour, and a
`[frames x tracks]` matrix — the shape all the per-track maths wants — cannot hold both. Forcing it
means either resampling one colour onto the other's grid, which invents localizations that were never
measured, or padding with NaN, after which row 7 silently means two different moments.

Therefore:

* **per-channel arrays stay per-channel.** Tracking, MSD, D, step vectors — all of it is computed on
  one colour's own grid and never merged.
* **the merged view is a flat list of spots, where the shared moment is a COLUMN and not an index.**
  That representation has no frame-axis problem, which is exactly why it is the one that can be merged.

One is canonical, the other derived from it. Neither is optional.

## The shared index is an integer the microscope wrote

The column that makes a merge mean anything is the **acquisition timepoint** — not a time in seconds.

ImageJ writes a label per page (TIFF tag 50839) and a hyperstack acquisition leaves its own indexing
in it:

```
c:2/4 t:1/2000 - HVK-3C-Plate1-296-011 #1
c:4/4 t:1/2000 - HVK-3C-Plate1-296-011 #1
c:2/4 t:2/2000 - ...
```

Channel and timepoint, per page, as integers, **written by the microscope rather than inferred from
page parity**. `c:2 t:5` and `c:4 t:5` are one moment of the experiment imaged twice. `dc_tiff_labels`
reads them; `dc_channels('fromStack', ...)` turns them into, per colour, which pages it occupies and
which timepoint each of its frames belongs to.

There are **no per-page timestamps** in these files. Every page carries the same ImageJ header with a
single nominal `finterval`, so every second in the old model was a frame number multiplied by a
nominal interval — derived, and accumulating. The 1.1% between the `0.054` the first attempt typed
into `cs_calib.mat` and the `0.0534237` its builds carried is half a second across a 1,791-frame
track. The timepoint index has no such error: it is what the microscope recorded.

It is also the unit a **tracking gap** is naturally counted in. Frame k+2 after frame k is a gap of
one, with nothing to compare against a tolerance. The whole apparatus of matching two colours with a
tolerance — and the boundary bug that apparatus had, where a partner landing exactly on half the
spacing was rejected — disappears: two colours are compared by set intersection on their timepoints.

`dt_s` survives in the registry, per colour, because a diffusion coefficient has to come out in
µm²/s eventually. It is applied **once, at the end, to a frame count**. Nothing indexes by it, which
is why a colour whose interval is unknown still detects and tracks.

## The structure

```
DCS                      one cell, both colours
  .cell                  char, the cell key — shared; the colours are of the SAME field of view
  .channels  1xC struct  THE REGISTRY — the only place a colour's maps live
      .key     'c2'      from the acquisition's own channel number; the flag carried everywhere else
      .label   'channel 2'
      .pages   nFrames x 1   the 1-based page of the stack for each of this colour's frames
      .tp      nFrames x 1   the acquisition timepoint each of those frames belongs to
      .dt_s    0.0534237     seconds per frame — for REPORTING a result, never an index
      .nFrames
  .tracks    1xC struct  canonical, per channel, on its own frame axis
      .tracks  1xN cell   the chains as the tracker returned them
      .nTracks .idOffset  the colour's own 0-based numbering, and its offset into the cell's
  .spots     table       DERIVED, both colours in one list
      .ch      categorical  the flag
      .frame   0-based within that colour
      .tp      the acquisition timepoint — what makes a merge exact
      .page    which page of the stack this detection was read from
      .x .y .q
      .trackId    unique across the cell
      .trackLocal this colour's own index, which finds the track in .tracks
      .spotId
  .masks                 shared by both colours: they are properties of the cell
      .er .mito          stacks + the averaging WINDOW, in timepoints, looked up by timepoint
```

## Why the flag is a key and not a 1 or a 2

`ch` is `'c2'`/`'c4'`, stored as a categorical whose categories are the registry's keys. A bare 1/2 is
an index into something the file does not carry: subset a table, drop a colour, add a third, and the
integers still parse but mean something else. A key is self-describing, survives subsetting, and fails
loudly rather than quietly when it does not match the registry. Taking it from the acquisition's own
channel number means the name in the data is the name in the microscope's label.

## Why a spot carries its page and timepoint, and why that is checked exactly

The spot carries `tp` and `page` so that drawing, merging and going back to the pixels need no
lookup. That is a denormalization, and it is safe only because one function builds it and `validate`
asserts **integer equality** against the channel map: a spot's timepoint and page are exactly what the
map says for its frame.

Integer equality, not a tolerance. There is nothing approximate about an index, so a tolerance here
could only hide the bug it was meant to catch — which is precisely what the first attempt's second
copy of `dt` did.

## Why track ids are renumbered on the way in

Each colour is tracked on its own, so each colour's tracker numbers from zero and **"track 0" exists
in both** — the ids collide for the most ordinary reason there is. They are made global as a colour is
added, so one number identifies one track anywhere in the cell and a list, an export or a selection
never has to carry the colour beside it. The colour's own index stays as `trackLocal`, because a track
has to be findable in the cell array the tracker returned, and `idOffset` records the mapping instead
of leaving it to be re-derived. `validate` pins both: a global id is its local index plus the offset,
and a track id never appears in two colours.

## Why one set of masks serves both colours, and why the lookup takes a timepoint

ER and mito belong to the cell, not to a colour, and **one segmentation serves both**. The
justification is a separation of timescales: the organelles are **window-averaged** over many
timepoints, while the particles are tracked at every one. Both colours of a timepoint fall inside the
same averaging window, so the mask that is right for one is right for the other.

That is a statement about the biology rather than about the software, so it is **recorded and checked**
rather than assumed. `dc_masks('check', D)` reports the window and says either

> one mask serves both colours: each page averages 100 timepoints, so both colours of any timepoint
> fall well inside the same window and the organelle has not moved between them

or, if someone later images the organelle fast, that at this window the assumption is no longer safe.
It is the kind of assumption that stays true until an acquisition changes and then fails silently.

The lookup takes a **timepoint** for a separate reason. A mask page belongs to a moment of the
acquisition; each colour samples the acquisition differently, so `page = frame + 1` gives two different
answers for one moment and neither is the mask's own numbering. The timepoint is the only index the two
colours share. (This is also the bug the single-colour pipeline still has on de-interleaved runs, where
every viewer is wrong by the stride.)

A page is a **window of timepoints, not an instant**: page p covers `[tp0 + (p-1)*w, tp0 + p*w)`, and
which page a timepoint belongs to is integer division — containment, not nearest-centre, which would
be a half-window error at every boundary.

The window is **declared, not derived**. Deriving it from page counts needs an exact integer ratio
between the stacks; real acquisitions rarely give one, and when the ratio is not whole the fallback is
silently 1, which leaves every timepoint past the end of the organelle stack with no mask at all.

What the averaging costs is worth stating, because it is a floor on what any organelle distance can
mean: a distance to a window-averaged mask is a distance to where the organelle was **on average** over
that window. If it moved within the window the mask is blurred relative to any one moment. `check`
reports the window so that can be quoted rather than forgotten.

## What the merge panel must say

"Do the two colours have the same frame number or not" is really "do their frames fall on the same
timepoints", which is set comparison rather than a tolerance. `dc_align` returns which of four
relations holds, with the numbers, and the panel says it in words:

| relation | what it means | what the panel does |
|---|---|---|
| `matched` | both colours have a frame at every timepoint the other does | frame numbers correspond 1:1 — **and the page gap is shown, because the two exposures are consecutive pages, not one instant** |
| `subsampled` | one colour's timepoints are every k-th of the other's | only 1 in k timepoints has both; step in the fast colour and say the other has nothing here |
| `partial` | they overlap but neither contains the other | frame numbers mean nothing; only the shared timepoints can be compared |
| `disjoint` | no timepoint is shared | nothing can be compared frame to frame |

Note what `matched` does **not** mean. An interleaved acquisition shares every timepoint and has the
same frame count in both colours — the case a panel stepping by frame index gets silently wrong is
gone, but the two exposures are still consecutive pages. `pageGap` is reported so that difference is
visible rather than implied by the word.

The rule underneath: **the merged view is driven by the timepoint, never by frame index.** A frame
query is resolved through the reference colour's own map into a timepoint, and every other colour is
taken at THAT timepoint.

Where a colour has no frame at the displayed timepoint, the panel says so rather than drawing its last
known position as if it were current. Under `subsampled` that is the normal case for half the
timepoints, not an error — and a mark drawn from a neighbouring timepoint would sit beside a current
one and look like a coincidence that was never observed.
