#!/usr/bin/env python3
"""Generate debug output image for the paper."""

import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Debug output text (sample)
debug_output = """$ pytest --fast-tcp --fast-tcp-debug tests/

[DEBUG] ═══════════════════════════════════════════
[DEBUG] Detected snapshot root at /project
[DEBUG] Snapshot diff resolved: old=a1b2c3d new=e4f5g6h
[DEBUG] Filesystem changes: +2 ~3 -1
[DEBUG]   [ADDED] tests/test_new_feature.py
[DEBUG]   [ADDED] tests/test_validation.py
[DEBUG]   [MODIFIED] tests/test_auth.py
[DEBUG]   [MODIFIED] tests/test_api.py
[DEBUG]   [MODIFIED] tests/test_utils.py
[DEBUG]   [DELETED] tests/test_deprecated.py
[DEBUG] Partition summary: 25 new / 42 unchanged / 5 deleted
[DEBUG] ═══════════════════════════════════════════
[DEBUG] Partition time:     0.0001s
[DEBUG] Preparation time:   0.0291s
[DEBUG] Prioritization time: 0.0086s
[DEBUG] Total FAST TCP time: 0.3782s
[DEBUG] ═══════════════════════════════════════════

========================= test session starts ==========================
platform darwin -- Python 3.11.0, pytest-7.4.0, pluggy-1.0.0
plugins: fast-tcp-0.1.0
collected 72 items (42 unchanged, 25 new, 5 deleted)

tests/test_new_feature.py::test_create ...
"""

# Create figure with terminal-like appearance
fig, ax = plt.subplots(figsize=(10, 7))
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.set_aspect('equal')

# Dark terminal background
terminal_bg = patches.FancyBboxPatch(
    (0.2, 0.2), 9.6, 9.6,
    boxstyle="round,pad=0.05,rounding_size=0.2",
    facecolor='#1e1e1e',
    edgecolor='#444444',
    linewidth=2
)
ax.add_patch(terminal_bg)

# Terminal title bar
title_bar = patches.FancyBboxPatch(
    (0.2, 9.3), 9.6, 0.5,
    boxstyle="round,pad=0.02,rounding_size=0.1",
    facecolor='#2d2d2d',
    edgecolor='#444444',
    linewidth=1
)
ax.add_patch(title_bar)

# Window buttons
for i, color in enumerate(['#ff5f56', '#ffbd2e', '#27c93f']):
    circle = patches.Circle((0.6 + i * 0.35, 9.55), 0.1, color=color)
    ax.add_patch(circle)

# Title text
ax.text(5, 9.55, 'Terminal - pytest --fast-tcp', fontsize=10, 
        color='#cccccc', ha='center', va='center', fontfamily='monospace')

# Terminal content
lines = debug_output.strip().split('\n')
y_pos = 9.0
line_height = 0.32

for line in lines[:25]:  # Limit to visible lines
    # Color coding
    if line.startswith('$'):
        color = '#27c93f'  # Green for prompt
    elif '[DEBUG]' in line:
        if '[ADDED]' in line:
            color = '#27c93f'  # Green for added
        elif '[MODIFIED]' in line:
            color = '#ffbd2e'  # Yellow for modified
        elif '[DELETED]' in line:
            color = '#ff5f56'  # Red for deleted
        elif 'time:' in line.lower() or 'Time:' in line:
            color = '#61afef'  # Blue for timing
        elif '═' in line:
            color = '#666666'  # Gray for separators
        else:
            color = '#98c379'  # Light green for debug
    elif 'collected' in line or 'platform' in line or 'plugins' in line:
        color = '#c678dd'  # Purple for pytest info
    elif '====' in line:
        color = '#61afef'  # Blue for pytest headers
    else:
        color = '#abb2bf'  # Default gray
    
    ax.text(0.4, y_pos, line[:80], fontsize=7, color=color,
            fontfamily='monospace', va='top')
    y_pos -= line_height

ax.axis('off')
plt.tight_layout()
plt.savefig('research/figures/debug_output.png', dpi=300, 
            bbox_inches='tight', facecolor='white', pad_inches=0.1)
plt.close()

print("Generated: research/figures/debug_output.png")
