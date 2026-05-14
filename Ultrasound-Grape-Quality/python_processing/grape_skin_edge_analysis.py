"""
Grape Skin Edge Blur Analysis from Ultrasound Images
Version 3 — Edge-Top Sampling with Horizontal Triplication

Methodology:
    - Identify the apex of the grape skin boundary (minimum Y on the largest contour)
    - Define three analysis positions: Apex, Left (−spacing), Right (+spacing)
    - At each position, extract three collinear horizontal points: L, C, R
    - Total: 3 positions × 3 points = 9 sampling points per image
    - Vertical intensity profiles (ESF) are extracted at each point
    - Edge Width (10–90% criterion) and Gradient Peak are derived from each ESF
    - Blur Index = Edge Width / Gradient Peak
"""

import numpy as np
import cv2
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from ipywidgets import FileUpload, VBox, Button, Output, Label, HBox, IntSlider, FloatSlider
from IPython.display import display
from PIL import Image
import io

# ============================================================================
# GLOBAL STATE
# ============================================================================

uploaded_images  = []
analysis_results = []

# Region of interest (absolute pixel coordinates)
X1, Y1 = 318, 110
X2, Y2 = 991, 300

# Horizontal analysis bounds within the ROI
X_MIN = 400
X_MAX = 800

# Physical scale
MM_PER_PIXEL = 0.056497

# ============================================================================
# WIDGETS
# ============================================================================

uploader = FileUpload(accept='image/*', multiple=True)

srad_iterations = IntSlider(
    value=15, min=5, max=30, step=1,
    description='SRAD Iterations:',
    style={'description_width': 'initial'}
)

srad_kappa = FloatSlider(
    value=30.0, min=10.0, max=50.0, step=1.0,
    description='SRAD Kappa:',
    style={'description_width': 'initial'}
)

percentile_threshold = IntSlider(
    value=95, min=90, max=99, step=1,
    description='Percentile Threshold:',
    style={'description_width': 'initial'}
)

position_spacing = IntSlider(
    value=3, min=3, max=10, step=1,
    description='Position Spacing (px):',
    style={'description_width': 'initial'}
)

horizontal_spacing = IntSlider(
    value=5, min=3, max=15, step=1,
    description='Horizontal Spacing (px):',
    style={'description_width': 'initial'}
)

profile_length = IntSlider(
    value=60, min=30, max=100, step=5,
    description='Profile Length (px):',
    style={'description_width': 'initial'}
)

analyze_button = Button(
    description='Analyze',
    button_style='success',
    layout={'width': '160px', 'height': '40px'}
)

out = Output()

# ============================================================================
# SPECKLE REDUCING ANISOTROPIC DIFFUSION (SRAD)
# ============================================================================

def srad_filter(img, num_iter=15, delta_t=0.15, kappa=30):
    """
    Speckle Reducing Anisotropic Diffusion filter.

    Parameters
    ----------
    img       : ndarray, uint8 grayscale input
    num_iter  : number of diffusion iterations
    delta_t   : time step
    kappa     : diffusion sensitivity (edge threshold)

    Returns
    -------
    ndarray, uint8 filtered image
    """
    img_out = img.astype(np.float64).copy()

    for _ in range(num_iter):
        nabla_N = np.roll(img_out, -1, axis=0) - img_out
        nabla_S = np.roll(img_out,  1, axis=0) - img_out
        nabla_E = np.roll(img_out, -1, axis=1) - img_out
        nabla_W = np.roll(img_out,  1, axis=1) - img_out

        c_N = np.exp(-(nabla_N / kappa) ** 2)
        c_S = np.exp(-(nabla_S / kappa) ** 2)
        c_E = np.exp(-(nabla_E / kappa) ** 2)
        c_W = np.exp(-(nabla_W / kappa) ** 2)

        img_out += delta_t * (
            c_N * nabla_N + c_S * nabla_S +
            c_E * nabla_E + c_W * nabla_W
        )

    return np.clip(img_out, 0, 255).astype(np.uint8)

# ============================================================================
# APEX DETECTION AND SAMPLING POINT GENERATION
# ============================================================================

def get_three_positions_with_horizontal_line(center_x, center_y,
                                              pos_spacing=3, h_spacing=5):
    """
    Generate 9 sampling points from the detected apex.

    Layout per position (L–C–R, horizontal only; T/B omitted to avoid
    sampling outside the boundary arc):

        L   C   R

    Parameters
    ----------
    center_x, center_y : apex coordinates (minimum Y on the boundary)
    pos_spacing        : lateral offset between the three main positions (px)
    h_spacing          : spacing between L, C, R within each position (px)

    Returns
    -------
    list of (x, y, position_label, point_label)
    """
    positions = [
        (center_x,               center_y, 'Apex'),
        (center_x - pos_spacing, center_y, 'Left'),
        (center_x + pos_spacing, center_y, 'Right'),
    ]

    all_points = []
    for pos_x, pos_y, pos_label in positions:
        for pt_x, pt_label in [
            (pos_x - h_spacing, 'L'),
            (pos_x,             'C'),
            (pos_x + h_spacing, 'R'),
        ]:
            if X_MIN <= pt_x <= X_MAX and Y1 <= pos_y <= Y2:
                all_points.append((pt_x, pos_y, pos_label, pt_label))

    return all_points

# ============================================================================
# PROFILE EXTRACTION AND EDGE METRICS
# ============================================================================

def extract_edge_profile(img, x, y, length=60):
    """
    Extract a vertical intensity profile (ESF) centred at (x, y).

    Returns
    -------
    positions : ndarray or None
    profile   : ndarray (Savitzky–Golay smoothed) or None
    """
    h, w = img.shape
    y_start = max(0, y - length // 2)
    y_end   = min(h, y + length // 2)

    if x < 0 or x >= w or (y_end - y_start) < 10:
        return None, None

    profile   = img[y_start:y_end, x].astype(float)
    positions = np.arange(len(profile))

    if len(profile) > 11:
        profile = savgol_filter(profile, window_length=11, polyorder=3)

    return positions, profile


def compute_edge_metrics(positions, profile):
    """
    Compute Edge Width (10–90% criterion) and Gradient Peak from an ESF.

    Returns
    -------
    edge_width    : float (pixels) or None
    gradient_peak : float or None
    """
    if positions is None or profile is None or len(profile) < 10:
        return None, None

    gradient      = np.gradient(profile)
    gradient_peak = np.max(np.abs(gradient))

    min_val   = np.min(profile)
    max_val   = np.max(profile)
    range_val = max_val - min_val

    if range_val < 1:
        return None, None

    idx_10     = np.argmin(np.abs(profile - (min_val + 0.1 * range_val)))
    idx_90     = np.argmin(np.abs(profile - (min_val + 0.9 * range_val)))
    edge_width = max(abs(positions[idx_90] - positions[idx_10]), 0.5)

    return edge_width, gradient_peak

# ============================================================================
# SINGLE-IMAGE ANALYSIS PIPELINE
# ============================================================================

def analyze_single_image(img):
    """
    Full analysis pipeline for one ultrasound image.

    Returns
    -------
    dict with keys: original, filtered, roi_filtered, center, results, all_points
    """
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY) if len(img.shape) == 3 else img
    gray = cv2.normalize(gray, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)

    roi         = gray[Y1:Y2, X1:X2].copy()
    roi_filtered = srad_filter(roi, num_iter=srad_iterations.value,
                               kappa=srad_kappa.value)

    img_filtered_full              = np.zeros_like(gray)
    img_filtered_full[Y1:Y2, X1:X2] = roi_filtered

    # Locate apex (minimum Y) on the dominant boundary contour
    x_min_roi  = max(0, X_MIN - X1)
    x_max_roi  = min(roi_filtered.shape[1], X_MAX - X1)
    roi_search = roi_filtered[:, x_min_roi:x_max_roi]

    threshold   = np.percentile(roi_search, percentile_threshold.value)
    high_signal = (roi_search >= threshold).astype(np.uint8)

    kernel      = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    high_signal = cv2.morphologyEx(high_signal, cv2.MORPH_CLOSE, kernel)
    high_signal = cv2.morphologyEx(high_signal, cv2.MORPH_OPEN,  kernel)

    contours, _ = cv2.findContours(high_signal, cv2.RETR_EXTERNAL,
                                   cv2.CHAIN_APPROX_SIMPLE)

    if not contours:
        loc            = np.unravel_index(roi_search.argmax(), roi_search.shape)
        center_x_roi   = x_min_roi + loc[1]
        center_y_roi   = loc[0]
    else:
        largest = max(contours, key=cv2.contourArea)
        pts     = largest[:, 0, :]          # shape (N, 2)
        top_idx = np.argmin(pts[:, 1])
        center_x_roi = x_min_roi + pts[top_idx, 0]
        center_y_roi = pts[top_idx, 1]

    center_x = center_x_roi + X1
    center_y = center_y_roi + Y1

    # Generate sampling points
    all_points = get_three_positions_with_horizontal_line(
        center_x, center_y,
        pos_spacing=position_spacing.value,
        h_spacing=horizontal_spacing.value,
    )

    # Extract profiles and compute metrics
    results = []
    for i, (x, y, pos_label, pt_label) in enumerate(all_points):
        positions, profile = extract_edge_profile(
            img_filtered_full, x, y, length=profile_length.value
        )
        edge_width, grad_peak = compute_edge_metrics(positions, profile)

        if edge_width is not None and grad_peak is not None:
            results.append({
                'point':          i,
                'x':              x,
                'y':              y,
                'position_label': pos_label,
                'point_label':    pt_label,
                'edge_width':     edge_width,
                'edge_width_mm':  edge_width * MM_PER_PIXEL,
                'gradient_peak':  grad_peak,
                'profile':        profile,
                'positions':      positions,
            })

    return {
        'original':     gray,
        'filtered':     img_filtered_full,
        'roi_filtered': roi_filtered,
        'center':       (center_x, center_y),
        'results':      results,
        'all_points':   all_points,
    }

# ============================================================================
# VISUALIZATION
# ============================================================================

def visualize_single_result(result, filename, idx):
    """Produce an 8-panel diagnostic figure for one image."""
    from matplotlib.patches import Rectangle

    n_pts = len(result['results'])
    fig   = plt.figure(figsize=(20, 12))
    fig.suptitle(
        f'Image {idx}: {filename}  —  {n_pts} valid sampling points',
        fontsize=13, fontweight='bold'
    )

    cx, cy = result['center']
    colors = {'Apex': '#d62728', 'Left': '#1f77b4', 'Right': '#2ca02c'}

    def _add_roi_overlay(ax):
        ax.add_patch(Rectangle(
            (X1, Y1), X2 - X1, Y2 - Y1,
            linewidth=1.5, edgecolor='lime', facecolor='none'
        ))
        ax.axvline(X_MIN, color='cyan', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.axvline(X_MAX, color='cyan', linestyle='--', linewidth=1.5, alpha=0.7)

    # --- Panel 1: Original image ---
    ax1 = plt.subplot(2, 4, 1)
    ax1.imshow(result['original'], cmap='gray')
    _add_roi_overlay(ax1)
    ax1.plot(cx, cy, 'r*', markersize=14, zorder=5, label='Apex')
    for r in result['results']:
        sz = 10 if r['point_label'] == 'C' else 7
        ax1.plot(r['x'], r['y'], 'o', color=colors[r['position_label']],
                 markersize=sz, alpha=0.85, zorder=4)
    ax1.set_title('Original image', fontsize=11)
    ax1.legend(loc='upper right', fontsize=8)
    ax1.axis('off')

    # --- Panel 2: SRAD-filtered image ---
    ax2 = plt.subplot(2, 4, 2)
    ax2.imshow(result['filtered'], cmap='gray')
    _add_roi_overlay(ax2)
    ax2.plot(cx, cy, 'r*', markersize=14, zorder=5)
    for r in result['results']:
        sz = 10 if r['point_label'] == 'C' else 7
        ax2.plot(r['x'], r['y'], 'o', color=colors[r['position_label']],
                 markersize=sz, alpha=0.85, zorder=4)
    ax2.set_title('After SRAD filtering', fontsize=11)
    ax2.axis('off')

    # --- Panel 3: Close-up sampling pattern ---
    ax3 = plt.subplot(2, 4, 3)
    ax3.imshow(result['filtered'], cmap='gray')
    ax3.set_xlim(cx - 50, cx + 50)
    ax3.set_ylim(cy + 30,  cy - 30)
    for r in result['results']:
        sz = 12 if r['point_label'] == 'C' else 8
        c  = colors[r['position_label']]
        ax3.plot(r['x'], r['y'], 'o', color=c, markersize=sz, alpha=0.9)
        ax3.text(r['x'] + 2, r['y'] - 2,
                 f"{r['position_label'][0]}{r['point_label']}",
                 color='white', fontsize=8, fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor=c, alpha=0.7))
    for pos_label in ['Apex', 'Left', 'Right']:
        pts = sorted([r for r in result['results'] if r['position_label'] == pos_label],
                     key=lambda r: r['x'])
        if len(pts) == 3:
            ax3.plot([p['x'] for p in pts], [p['y'] for p in pts],
                     '-', color=colors[pos_label], alpha=0.5, linewidth=1.8)
    ax3.set_title('Sampling pattern (L–C–R per position)', fontsize=11)
    ax3.axis('off')

    # --- Panel 4: Point count per position ---
    ax4 = plt.subplot(2, 4, 4)
    pos_counts = {k: sum(1 for r in result['results'] if r['position_label'] == k)
                  for k in ['Apex', 'Left', 'Right']}
    labels_bar = list(pos_counts.keys())
    counts_bar = list(pos_counts.values())
    bar_colors = [colors[k] for k in labels_bar]
    bars = ax4.bar(labels_bar, counts_bar, color=bar_colors, alpha=0.75)
    for bar, cnt in zip(bars, counts_bar):
        ax4.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                 f'{int(cnt)}/3', ha='center', va='bottom',
                 fontsize=11, fontweight='bold')
    ax4.set_ylabel('Valid points')
    ax4.set_title('Points per position', fontsize=11)
    ax4.set_ylim(0, 3.7)
    ax4.grid(True, alpha=0.3, axis='y')

    # --- Panel 5: ESF profiles ---
    ax5 = plt.subplot(2, 4, 5)
    for r in result['results']:
        if r['profile'] is not None:
            ls = '-' if r['point_label'] == 'C' else '--'
            lw = 2.0 if r['point_label'] == 'C' else 1.4
            ax5.plot(r['positions'], r['profile'],
                     color=colors[r['position_label']], linestyle=ls,
                     linewidth=lw, alpha=0.75,
                     label=f"{r['position_label'][0]}{r['point_label']}")
    ax5.set_xlabel('Position (px)')
    ax5.set_ylabel('Intensity')
    ax5.set_title('Edge Spread Functions (ESF)', fontsize=11)
    ax5.legend(fontsize=7, ncol=3)
    ax5.grid(True, alpha=0.3)

    # --- Panel 6: LSF (gradient of ESF) ---
    ax6 = plt.subplot(2, 4, 6)
    for r in result['results']:
        if r['profile'] is not None:
            ls = '-' if r['point_label'] == 'C' else '--'
            lw = 2.0 if r['point_label'] == 'C' else 1.4
            ax6.plot(r['positions'], np.gradient(r['profile']),
                     color=colors[r['position_label']], linestyle=ls,
                     linewidth=lw, alpha=0.75,
                     label=f"{r['position_label'][0]}{r['point_label']}")
    ax6.set_xlabel('Position (px)')
    ax6.set_ylabel('dI/dx')
    ax6.set_title('Line Spread Functions (LSF)', fontsize=11)
    ax6.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax6.legend(fontsize=7, ncol=3)
    ax6.grid(True, alpha=0.3)

    # --- Panel 7: Edge Width by position ---
    ax7 = plt.subplot(2, 4, 7)
    pos_ew = {k: [r['edge_width'] for r in result['results']
                  if r['position_label'] == k]
              for k in ['Apex', 'Left', 'Right']}
    valid_pos = [k for k in ['Apex', 'Left', 'Right'] if pos_ew[k]]
    means_ew  = [np.mean(pos_ew[k]) for k in valid_pos]
    stds_ew   = [np.std(pos_ew[k]) if len(pos_ew[k]) > 1 else 0 for k in valid_pos]
    x_pos     = np.arange(len(valid_pos))
    ax7.bar(x_pos, means_ew, yerr=stds_ew, capsize=5,
            color=[colors[k] for k in valid_pos], alpha=0.75)
    overall_ew = np.mean([r['edge_width'] for r in result['results']])
    ax7.axhline(overall_ew, color='k', linestyle='--', linewidth=1.8,
                label=f'Mean: {overall_ew:.2f} px')
    ax7.set_xticks(x_pos)
    ax7.set_xticklabels(valid_pos)
    ax7.set_ylabel('Edge Width (px)')
    ax7.set_title('Edge Width by position', fontsize=11)
    ax7.legend(fontsize=9)
    ax7.grid(True, alpha=0.3, axis='y')

    # --- Panel 8: Gradient Peak by position ---
    ax8 = plt.subplot(2, 4, 8)
    pos_gp   = {k: [r['gradient_peak'] for r in result['results']
                    if r['position_label'] == k]
                for k in ['Apex', 'Left', 'Right']}
    means_gp = [np.mean(pos_gp[k]) for k in valid_pos]
    stds_gp  = [np.std(pos_gp[k]) if len(pos_gp[k]) > 1 else 0 for k in valid_pos]
    ax8.bar(x_pos, means_gp, yerr=stds_gp, capsize=5,
            color=[colors[k] for k in valid_pos], alpha=0.75)
    overall_gp = np.mean([r['gradient_peak'] for r in result['results']])
    ax8.axhline(overall_gp, color='k', linestyle='--', linewidth=1.8,
                label=f'Mean: {overall_gp:.2f}')
    ax8.set_xticks(x_pos)
    ax8.set_xticklabels(valid_pos)
    ax8.set_ylabel('Gradient Peak')
    ax8.set_title('Gradient Peak by position', fontsize=11)
    ax8.legend(fontsize=9)
    ax8.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.show()


def display_summary_results():
    """Print aggregate statistics and produce a cross-sample comparison figure."""
    n = len(analysis_results)

    all_ew    = [r['edge_width_mean']    for r in analysis_results]
    all_ew_mm = [r['edge_width_mean_mm'] for r in analysis_results]
    all_gp    = [r['gradient_peak_mean'] for r in analysis_results]
    all_bi    = [r['blur_index']         for r in analysis_results]

    def _stat(vals):
        m = np.mean(vals)
        s = np.std(vals, ddof=1) if len(vals) > 1 else 0.0
        return m, s

    ew_m,  ew_s  = _stat(all_ew)
    ewmm_m, ewmm_s = _stat(all_ew_mm)
    gp_m,  gp_s  = _stat(all_gp)
    bi_m,  bi_s  = _stat(all_bi)

    print(f"\nSummary ({n} samples)\n{'─'*60}")
    print(f"  Edge Width     : {ew_m:.3f} ± {ew_s:.3f} px"
          f"  /  {ewmm_m:.3f} ± {ewmm_s:.3f} mm")
    print(f"  Gradient Peak  : {gp_m:.3f} ± {gp_s:.3f}")
    print(f"  Blur Index     : {bi_m:.4f} ± {bi_s:.4f}")
    print(f"\nPer-sample detail\n{'─'*60}")

    for i, r in enumerate(analysis_results, 1):
        print(f"  [{i}] {r['filename']}")
        print(f"      n = {r['n_points']}/9 points")
        print(f"      EW = {r['edge_width_mean']:.3f} ± {r['edge_width_std']:.3f} px"
              f"  /  {r['edge_width_mean_mm']:.3f} ± {r['edge_width_std_mm']:.3f} mm")
        print(f"      GP = {r['gradient_peak_mean']:.3f} ± {r['gradient_peak_std']:.3f}")
        print(f"      BI = {r['blur_index']:.4f}")

    if n <= 1:
        return

    labels = [r['filename'][:14] + '…' if len(r['filename']) > 14
              else r['filename'] for r in analysis_results]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle('Cross-sample Comparison', fontsize=13, fontweight='bold')

    for ax, data, mean_val, std_val, ylabel, title in [
        (axes[0], all_ew,    ew_m,  ew_s,
         'Edge Width (px)',  f'Edge Width\n({ewmm_m:.3f} ± {ewmm_s:.3f} mm)'),
        (axes[1], all_gp,    gp_m,  gp_s,
         'Gradient Peak',    'Gradient Peak'),
        (axes[2], all_bi,    bi_m,  bi_s,
         'Blur Index',       'Blur Index  (EW / GP)'),
    ]:
        ax.bar(range(n), data, alpha=0.72)
        ax.axhline(mean_val, color='r', linestyle='--', linewidth=1.8,
                   label=f'Mean: {mean_val:.3f}')
        ax.fill_between(range(n), mean_val - std_val, mean_val + std_val,
                        alpha=0.15, color='r')
        ax.set_xlabel('Sample')
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xticks(range(n))
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.show()

# ============================================================================
# ANALYSIS CALLBACK
# ============================================================================

def on_analyze_click(b):
    global uploaded_images, analysis_results

    with out:
        out.clear_output(wait=True)

        if not uploader.value:
            print('No images uploaded.')
            return

        print('=' * 70)
        print('Grape Skin Edge Blur Analysis  —  V3')
        print('=' * 70)
        print(f'Images          : {len(uploader.value)}')
        print(f'SRAD            : iterations={srad_iterations.value}, '
              f'kappa={srad_kappa.value}')
        print(f'Positions       : Apex + Left/Right (±{position_spacing.value} px)')
        print(f'Horizontal step : {horizontal_spacing.value} px  (L, C, R)')
        print(f'Sampling points : 3 positions × 3 points = 9 per image')
        print()

        uploaded_images  = []
        analysis_results = []

        for idx, (filename, file_info) in enumerate(uploader.value.items(), start=1):
            print(f'{'─'*70}')
            print(f'Image {idx}/{len(uploader.value)}: {filename}')

            try:
                image = Image.open(io.BytesIO(file_info['content'])).convert('RGB')
                img   = cv2.cvtColor(np.array(image), cv2.COLOR_RGB2BGR)
                uploaded_images.append((filename, img))

                result = analyze_single_image(img)

                if not result['results']:
                    print('  No valid boundary found.')
                    continue

                ews   = [r['edge_width']    for r in result['results']]
                ewmms = [r['edge_width_mm'] for r in result['results']]
                gps   = [r['gradient_peak'] for r in result['results']]

                ew_m   = np.mean(ews);   ew_s   = np.std(ews,   ddof=1) if len(ews)   > 1 else 0
                ewmm_m = np.mean(ewmms); ewmm_s = np.std(ewmms, ddof=1) if len(ewmms) > 1 else 0
                gp_m   = np.mean(gps);   gp_s   = np.std(gps,   ddof=1) if len(gps)   > 1 else 0
                bi     = ew_m / gp_m

                print(f'  Valid points   : {len(result["results"])}/9')
                print(f'  Apex           : {result["center"]}')
                print(f'  Edge Width     : {ew_m:.3f} ± {ew_s:.3f} px'
                      f'  /  {ewmm_m:.3f} ± {ewmm_s:.3f} mm')
                print(f'  Gradient Peak  : {gp_m:.3f} ± {gp_s:.3f}')
                print(f'  Blur Index     : {bi:.4f}')

                analysis_results.append({
                    'filename':           filename,
                    'edge_width_mean':    ew_m,
                    'edge_width_std':     ew_s,
                    'edge_width_mean_mm': ewmm_m,
                    'edge_width_std_mm':  ewmm_s,
                    'gradient_peak_mean': gp_m,
                    'gradient_peak_std':  gp_s,
                    'blur_index':         bi,
                    'n_points':           len(result['results']),
                    'details':            result,
                })

                visualize_single_result(result, filename, idx)

            except Exception as e:
                import traceback
                print(f'  Error: {e}')
                traceback.print_exc()

        if analysis_results:
            print('\n' + '=' * 70)
            display_summary_results()
        else:
            print('No images were successfully analysed.')


analyze_button.on_click(on_analyze_click)

# ============================================================================
# INTERFACE
# ============================================================================

print('Grape Skin Edge Blur Analysis  —  V3')
print('─' * 70)
print(f'ROI              : ({X1},{Y1}) – ({X2},{Y2})')
print(f'Analysis zone X  : {X_MIN} – {X_MAX} px  ({X_MAX - X_MIN} px wide)')
print(f'Scale            : {MM_PER_PIXEL} mm/px')
print(f'Metrics          : Edge Width (10–90%), Gradient Peak, Blur Index')
print('─' * 70)

display(VBox([
    Label('Step 1 — Upload ultrasound images'),
    uploader,
    Label('Step 2 — Parameters', layout={'margin': '16px 0 8px 0'}),
    VBox([
        srad_iterations,
        srad_kappa,
        percentile_threshold,
        position_spacing,
        horizontal_spacing,
        profile_length,
    ], layout={'padding': '10px', 'border': '1px solid #ccc', 'border_radius': '4px'}),
    HBox([analyze_button], layout={'margin': '16px 0', 'justify_content': 'center'}),
    out,
]))
