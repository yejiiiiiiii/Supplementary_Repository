#!/usr/bin/env python3
"""
Distance between two lines L(t) from video — Grey-spectrum two-line detector with LOCK-ON tracking

- Detects two nearly horizontal, dark lines in each video frame
- Turns the vertical distance between them into L_px, smooth it, convert to L_mm using L0, and save both a CSV and a monitoring video.

Row score:
  For each row, take the mean of the darkest p% pixels (keeps grey spectrum).
  Pencil lines give strong (dark) responses; seams/shadows are weaker.

Tracking:
  - Frame 1: pick two strongest peaks globally (with min separation)
  - Next frames: search only within ±TRACK_WIN_PX around previous y values
  - If a candidate jumps > MAX_JUMP_PX, ignore and keep previous
  - Enforce ordering and min separation

Outputs:
  CSV: t, L_mm, L_px, y_top_px, y_bot_px
  Monitoring video: <csv_stem>_monitor.mp4
"""

import sys, csv
from pathlib import Path
import numpy as np
import cv2 # open cv for video I/O and image processing

# ======== DEFAULTS (edit these so you can just run the file) ========
VIDEO_PATH    = r"C:\Users\yejihan\Downloads\3_3_1.mp4"
L0_MM         = 30.0 # initial gauge length (mm)
OUT_CSV       = r"C:\Users\yejihan\Downloads\L_vs_t.csv"
WRITE_MON_MP4 = True # save the video?
# ====================================================================

# -------- Detector tunables --------
P_DARK_FRAC   = 0.1   # darkest x% pixels per row
MIN_SEP_PX    = 30     # minimum vertical separation between the two lines (px)
GUARD_FRAC_TOP = 0.4 # top guard fraction. ignore this fraction of rows at top/bottom to avoid borders
GUARD_FRAC_BOT = 0.1 # bottom guard fraction
SMOOTH_GAUSS   = 3    # vertical smoothing of row scores, larger value = stronger tendency to merge nearby rows.
                       # Each row’s darkness score is recomputed as a weighted average of itself ±15 rows above and below, with Gaussian weights
# ------------ Tracking  ------------
TRACK_WIN_PX  = 2    # search window size (± px) around previous y
MAX_JUMP_PX   = 2  # ignore bigger-than-this jumps
# ----------- Smoothing  ------------
MEDIAN_K      = 5     # looks at the last 10 frames of L_px and outputs the median.
EMA_ALPHA     = 0.05   # EMA low-pass for L_px (smaller = smoother). L_ema = alpha * L(t) + (1-alpha) * L(t-1)
CALIB_FRAMES  = 5      # frames for mm/px calibration. mm/px = L_0 / avg(L_px[0:x])
# ------------------------------------

def preprocess_gray(gray):
    """contrast enhancement + denoising"""
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8)) # Contrast Limited Adaptive Histogram Equalization.
    g = clahe.apply(gray) # boost contrast locally so lines pop out
    g = cv2.bilateralFilter(g, d=7, sigmaColor=30, sigmaSpace=7) # reduce noise while keping edges sharp
    return g

def row_dark_score(gray, frac=P_DARK_FRAC, apply_guard=True):
    """
    Continuous row darkness score:
    For each row, average the darkest frac fraction of pixels.
    Return darkness so larger = darker line.
    """
    H, W = gray.shape
    k = max(1, int(round(frac * W)))
    part = np.partition(gray, kth=k-1, axis=1)[:, :k]
    darkest_mean = part.mean(axis=1).astype(np.float32)
    dark_score = 255.0 - darkest_mean
    dark_sm = cv2.GaussianBlur(dark_score.reshape(-1,1), (1, SMOOTH_GAUSS), 0).ravel()

    if apply_guard:
        g_top = int(round(GUARD_FRAC_TOP * H))
        g_bot = int(round(GUARD_FRAC_BOT * H))
        if g_top > 0:
            dark_sm[:g_top] = 0
        if g_bot > 0:
            dark_sm[H-g_bot:] = 0
    return dark_sm

def pick_two_peaks_global(score, min_sep=MIN_SEP_PX):
    """Pick two strongest peaks in 1D score with minimum separation (global init)."""
    idx = np.argsort(-score)
    chosen = []
    for y in idx:
        if score[y] <= 0:
            break
        if not chosen or all(abs(int(y)-int(z)) >= min_sep for z in chosen):
            chosen.append(int(y))
        if len(chosen) == 2:
            break
    if len(chosen) < 2:
        top = idx[:min(50, len(score))]
        if len(top) >= 2:
            return int(top.min()), int(top.max())
        else:
            return 0, max(0, len(score)-1)
    return tuple(sorted(chosen)) # [top bottom]

def pick_peak_near(score, y_prev, win_px=TRACK_WIN_PX):
    """Pick peak position within a vertical window around y_prev."""
    H = len(score) # height of the image
    a = max(0, int(round(y_prev)) - win_px)
    b = min(H, int(round(y_prev)) + win_px + 1)
    if b <= a + 1:
        return int(round(y_prev))
    y = int(np.argmax(score[a:b])) + a
    return y

def median_push(buf, k, x):
    '''
    keeps a sliding window of the last k values and returns their median. (oldest value is dumped as new value comes in)
    Smoothens noisy signals frame-by-frame -> makes the L_px measurement smoother by filtering out sudden noise.

    buf: a buffer (storage) of the most recent values.
    k: the window size (how many past values we want to keep).
    x: the new incoming value (for this frame).
    '''
    from collections import deque
    if buf is None:
        buf = deque(maxlen=k)
    buf.append(float(x))
    return buf, float(np.median(buf))

def process(video_path, L0_mm, out_csv, write_monitor=True):
    cap = cv2.VideoCapture(str(video_path))
    if not cap.isOpened():
        raise RuntimeError(f"Failed to open video: {video_path}")
    fps = cap.get(cv2.CAP_PROP_FPS) or 30.0
    W   = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    H   = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

    # monitoring writer
    mon = None
    if write_monitor:
        mon_path = str(Path(out_csv).with_name(Path(out_csv).stem + "_monitor.mp4")) # output file name
        fourcc = cv2.VideoWriter_fourcc(*"mp4v") # defines the codec for the output video
        mon = cv2.VideoWriter(mon_path, fourcc, fps, (W, H))
        if not mon.isOpened():
            mon = None

    rows = [] # [t, L_mm, ema_L, y_top, y_bot]
    median_buf = None
    ema_L = None
    y_top_prev = None
    y_bot_prev = None

    # main processing loop, runs frame-by-frame
    i = 0
    while True:
        ok, frame = cap.read() # reads one frame at a time from the video cap
        if not ok:
            break
        i += 1
        t = (i-1) / fps # time stamp

        # converts the frame to grayscale
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        g = preprocess_gray(gray) # CLAHE + denoising (bilateral filter)

        # --- detection with lock-on tracking ---
        if y_top_prev is None or y_bot_prev is None: # no history yet
            # first frame: apply guard, global init
            score = row_dark_score(g, frac=P_DARK_FRAC, apply_guard=True) # row scoring and Gaussian smoothing
            y_top, y_bot = pick_two_peaks_global(score, min_sep=MIN_SEP_PX)
        else:
            # later frames: no guard, local search near previous positions    
            score = row_dark_score(g, frac=P_DARK_FRAC, apply_guard=False)
            y_top = pick_peak_near(score, y_top_prev, TRACK_WIN_PX)
            y_bot = pick_peak_near(score, y_bot_prev, TRACK_WIN_PX)

            # outlier guard: ignore large jumps
            if abs(y_top - y_top_prev) > MAX_JUMP_PX:
                y_top = int(round(y_top_prev))
            if abs(y_bot - y_bot_prev) > MAX_JUMP_PX:
                y_bot = int(round(y_bot_prev))

            # maintain ordering and minimum separation
            if y_bot < y_top:
                y_top, y_bot = y_bot, y_top
            if (y_bot - y_top) < MIN_SEP_PX:
                mid = (y_top + y_bot) // 2
                y_top = max(0, mid - MIN_SEP_PX // 2)
                y_bot = min(H-1, y_top + MIN_SEP_PX)

        # update previous (for next frame)
        y_top_prev, y_bot_prev = float(y_top), float(y_bot)

        # length in pixels
        Lpx_raw = float(y_bot - y_top)

        # temporal smoothing of L
        median_buf, Lpx_med = median_push(median_buf, MEDIAN_K, Lpx_raw)
        if ema_L is None:
            ema_L = Lpx_med
        else:
            ema_L = EMA_ALPHA * Lpx_med + (1.0 - EMA_ALPHA) * ema_L # smoothed vertical distance between two lines

        rows.append([t, None, ema_L, float(y_top), float(y_bot)])

        # draw overlay
        vis = frame.copy() # makes a copy of the current video frame, acting like a per-frame buffer
        cv2.line(vis, (0, y_top), (W-1, y_top), (0, 0, 255), 2)
        cv2.line(vis, (0, y_bot), (W-1, y_bot), (0, 0, 255), 2)
        cv2.putText(vis, f"t={t:.3f}s  Lpx={ema_L:.1f}", (20, 40),
                    cv2.FONT_HERSHEY_SIMPLEX, 1.0, (255,255,255), 2)

        # tiny inset: normalized row score (for debugging)
        s = (score - score.min()) / (np.ptp(score) + 1e-6)
        bar = (255*(1.0 - s)).astype(np.uint8)
        bar = np.repeat(bar.reshape(-1,1), max(2, W//80), axis=1)
        bar = cv2.cvtColor(bar, cv2.COLOR_GRAY2BGR)
        vis[:, W-bar.shape[1]:W] = bar

        if mon is not None:
            mon.write(vis) # writing of frames

        cv2.imshow("Grey-spectrum two-line tracker", vis)
        if cv2.waitKey(1) & 0xFF == 27:
            break

    cap.release()
    cv2.destroyAllWindows()
    if mon is not None:
        mon.release()

    if not rows:
        raise RuntimeError("No frames processed.")

    # calibrate mm/px from first N frames
    N = min(CALIB_FRAMES, len(rows))
    L0px = float(np.mean([r[2] for r in rows[:N]])) # extracts the smoothed pixel length (ema_L) from the first N frames.
    mm_per_px = (L0_mm / L0px) if L0px > 0 else 1.0

    # fill L_mm and write CSV
    for r in rows:
        r[1] = r[2] * mm_per_px

    out_path = Path(out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t", "L_mm", "L_px", "y_top_px", "y_bot_px"])
        w.writerows(rows)

    print(f"Saved CSV: {out_path}")
    if mon is not None:
        print(f"Saved monitor video: {Path(out_csv).with_name(Path(out_csv).stem + '_monitor.mp4')}")
    print(f"mm/px = {mm_per_px:.6f}, frames = {len(rows)}")

if __name__ == "__main__":
    # allow optional positional overrides: script.py <video> <L0_mm> <out_csv>
    video = VIDEO_PATH
    L0    = L0_MM
    out   = OUT_CSV
    if len(sys.argv) >= 2: video = sys.argv[1]
    if len(sys.argv) >= 3: L0 = float(sys.argv[2])
    if len(sys.argv) >= 4: out = sys.argv[3]

    try:
        process(video, L0, out, write_monitor=WRITE_MON_MP4)
    except Exception as e:
        print(f"[Error] {e}")
        sys.exit(1)