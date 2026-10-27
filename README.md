# DLC_track_pharynx_egg

## Chemotaxis parameters (`initialize_load_files_high_res_universal_bootstrapping.py`)

This section describes every column the bootstrapping wrapper writes to
`chemotaxis_params.csv` / `chemotaxis_analysis_complete.h5`: how each one is calculated, over
which time window, and whether it is smoothed. The maths behind the core functions is in
[chemotaxis_analysis_high_res/chemotaxis_calculations.md](chemotaxis_analysis_high_res/chemotaxis_calculations.md).

Notation: `f` = `factor_px_to_mm`, `fps` = frame rate, `(cx, cy)` = `video_resolution / 2`,
`{c}` = center spline point (`n_spline_points / 2`), `{p}` = DLC body point (`nose`, `pharynx`, `tail`).

### Pipeline order
1. Load annotations, spline K, spline X/Y, and stage position. The stage position is linearly
   resampled to the number of video frames if its length differs.
2. Stage coordinate transform into the arena frame (`X_rel`, `Y_rel`), plus the odor position.
3. Spline-point positions (centroid, point 0, point `{c}`).
4. Centroid smoothing (Savitzky–Golay).
5. Time and border distances.
6. DLC body points (`--DLC_nose`, `--DLC_pharynx`, `--DLC_tail`; any subset).
7. Behaviour annotations and reversals.
8. Movement and speed.
9. Odor parameters.
10. Smoothing (`*_smoothed` columns).
11. Bootstrap.
12. Export: CSV, HDF5, PDF, and optional movie.

### Positions (per frame, mm)
| Column | Calculation |
|---|---|
| `X`, `Y`, `frame` | Stage position from `*TablePosRecord.txt` (assumed mm for zim01/zim06). `frame` = row index for zim01/zim06. |
| `X_rel`, `Y_rel` | Stage position in the arena frame. **zim06:** `X_rel = abs(Y − top_left_y)`, `Y_rel = abs(X − top_left_x)`. **zim01:** `X_rel = top_left_y − Y`, `Y_rel = top_left_x − X`. **crop:** `(px − top_left)·f`, rotated 90° with `abs`. The swap or rotation is the same one applied to the odor, so it only changes the viewing orientation. |
| `odor_x`, `odor_y` | Odor position, transformed exactly like the stage. |
| `X/Y_rel_skel_pos_centroid` | Stage + (mean of all spline points − centre)·f, then linear interpolation of gaps and a Savitzky–Golay filter (window = `fps` frames, made odd; polynomial order 3). |
| `X/Y_rel_skel_pos_centroid_raw` | The same before the Savitzky–Golay filter (after interpolation). |
| `X/Y_rel_skel_pos_0`, `X/Y_rel_skel_pos_{c}` | Stage + (spline point − centre)·f. zim06: `X_rel + (sX − cx)·f`, `Y_rel − (sY − cy)·f`. The spline CSVs store (image row, image column), so this lands on the correct stage axes. |
| `X/Y_rel_DLC_{p}` | Stage + (DLC point − centre)·f. zim06: `X_rel + (y − cy)·f`, `Y_rel − (x − cx)·f`. DLC is (column, row), so the axes are swapped relative to the spline; this matches imutils `bodyparts_coordinates.py`. The raw (unfiltered) DLC `.h5` is used, with no likelihood filtering. |
| `skel_pos_abs` group (CSV) | `X/Y_rel_skel_pos_{k}` for every spline point `k`. |
| `Spline_K` group (CSV) | Spline curvature, as read from the input. |

### Time, border, behaviour (per frame)
| Column | Calculation |
|---|---|
| `time_seconds` | `frame / fps`. |
| `distance_to_border_centroid`, `distance_to_border_nose`, `distance_to_border_DLC_{p}` | Minimum distance to the 4 walls of the arena (default 38 × 40.5 mm). `_nose` = spline point 0. |
| `behaviour_state` | Reversal annotation (−1 = reversal). |
| turn columns | As read from `turns_annotation.csv`. |
| `reversal_onset` / `reversal_end` | 1 on the frame where `behaviour_state` enters or leaves −1. |
| `reversal_frequency` | Trailing rolling **sum** of `reversal_onset` over `fps·60` frames (reversals per minute). |

### Movement (per frame)
| Column | Calculation |
|---|---|
| `centroid_dx_dt`, `centroid_dy_dt` | `np.gradient` of the smoothed centroid: centred difference over ±1 frame (mm/frame). |
| `centroid_displacement_vector_degrees` | `atan2(dy, dx)` (heading, degrees). |
| `centroid_displacement_magnitude` | `sqrt(dx² + dy²)` (mm/frame). |
| `{c}_dx_dt`, `{c}_dy_dt`, `{c}_displacement_vector_degrees`, `{c}_displacement_magnitude` | The same for spline point `{c}` (no Savitzky–Golay filter). |
| `speed_centroid`, `speed_center_{c}` | Displacement magnitude × `fps` (mm/s). |
| `speed_centroid_signed`, `speed_center_{c}_signed` | Speed × first column of the reversal annotation. |
| `curving_angle` | Change in centroid heading over 1 frame, wrapped to [−180, 180]. |

### Odor (per frame; only with `--odor_pos` and the diffusion arrays)
| Column | Calculation |
|---|---|
| `distance_to_odor_stage`, `_centroid`, `_0`, `_DLC_{p}` | Euclidean distance from that position to (`odor_x`, `odor_y`). |
| `conc_at_centroid`, `conc_at_0`, `conc_at_DLC_{p}` | Diffusion lookup. Uses row `round(time_seconds) + diffusion_time_offset` of the arrays and the nearest distance bin, so time resolution is 1 s. |
| `dC_centroid` | `conc_at_centroid.diff(fps)`: change over the last **fps frames (1 s)**. |
| `dC_0`, `dC_DLC_{p}` | `.diff(dC_lookback_frames)`: change over the last `dC_lookback_frames` frames (default 1). |
| `d_DLC_nose_pharynx_C`, `d_DLC_nose_tail_C`, `d_DLC_pharynx_tail_C` | `conc_at_DLC_a − conc_at_DLC_b` at the same frame (spatial difference along the body). |
| `radial_speed` | `−np.gradient(distance_to_odor_centroid) · fps` (mm/s; positive = towards the odor). |
| `bearing_angle` | `abs(angle from centroid to odor − centroid heading)`, wrapped, in [0, 180]. |
| `NI` | `radial_speed / speed_centroid` (NaN where speed = 0). |

### Smoothed vs per-frame columns
Every column above is saved **per frame and complete**. The columns listed below **also** get a
`<column>_smoothed` copy:
1. Outliers are masked to NaN, **in the copy used for smoothing only**. An outlier is a value
   outside mean ± 2.576·SD of that column over the whole recording.
2. A centred rolling mean (`min_periods=1`) is taken over `window × fps` frames.

| Window | Columns |
|---|---|
| 1 s | `speed_centroid`, `speed_center_{c}`, `radial_speed`, `reversal_frequency`, `bearing_angle`, `NI`, `curving_angle`, `distance_to_odor_centroid`, `conc_at_centroid`, `conc_at_0`, `dC_centroid`, `dC_0`, `distance_to_odor_DLC_{p}`, `conc_at_DLC_{p}`, `dC_DLC_{p}`, `d_DLC_*_C` |
| 10 s (default) | `distance_to_border_centroid`, `distance_to_border_nose`, `distance_to_odor_0` (not in the 1 s list, so it falls back to 10 s) |

All other columns are per frame only: positions, `time_seconds`, DLC border distances,
behaviour, dx/dt and heading, signed speeds, and `distance_to_odor_stage`.

### Bootstrap (`--bootstrap_iterations N`)
For each of N random odor positions inside the arena, every odor column is recomputed and
smoothed as above, for the centroid, point 0 and all DLC points. The `d_DLC_*_C` columns are
not saved for bootstrap iterations. Results are saved only in the HDF5 under `bootstrap_analysis/random_odor_XXX`.
The HDF5 also stores the real-odor means of NI, radial speed, bearing, distance and centroid concentration.

### Outputs (written next to `--reversal_annotation`)
- `chemotaxis_params.csv`: column groups `chemotaxis_parameter` (everything above),
  `Spline_K`, `skel_pos_abs`.
- `chemotaxis_analysis_complete.h5`: groups `worm_positions`, `odor_independent_parameters`,
  `real_odor_analysis` (`distances`, `concentrations`, `navigation`), `bootstrap_analysis`,
  `skeleton_data`, `metadata`.
- `chemotaxis_analysis.pdf`: overview plots and time series of the smoothed columns.
- `worm_movie.avi` (with `--generate_video`): every `fps/4`-th frame, with the
  `ratiometric` signal from `--signal_csv_path` and `conc_at_DLC_nose`.
