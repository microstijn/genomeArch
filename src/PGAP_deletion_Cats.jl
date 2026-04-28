using CairoMakie

# 1. Apply the Custom Minimalist Theme
paper_theme = Theme(
    font = "Helvetica",
    fontsize = 14,
    Axis = (
        xgridvisible = false, ygridvisible = false,          
        topspinevisible = false, rightspinevisible = false,     
        spinewidth = 1.5, bottomspinecolor = :black,
        leftspinecolor = :black, xtickwidth = 1.5, ytickwidth = 1.5,
    )
)
set_theme!(paper_theme)

# Initialize Figure with a 2-Row Layout
fig = Figure(size = (1200, 900))

# ==============================================================================
# PANEL A: The Dual-Gated Survival Phase Space
# ==============================================================================
ax1 = Axis(fig[1, 1], 
    title = "A. The Dual-Gated Survival Logic (Truncation Phase)",
    xlabel = "Truncated Gene Length (L_Y')",
    ylabel = "Evidence Support Score (w_support)"
)

# Set axes limits
xlims!(ax1, 0, 250)
ylims!(ax1, 0, 45)

# Plot the three rigid geometric zones using Rect(x, y, width, height)
# 1. Deletion Void (< 60bp OR < 180bp with weak evidence)
poly!(ax1, Rect(0, 0, 60, 45), color = (:lightcoral, 0.4))
poly!(ax1, Rect(60, 0, 120, 25), color = (:lightcoral, 0.4))

# 2. High-Evidence Safe Zone (60-180bp AND w_support >= 25)
poly!(ax1, Rect(60, 25, 120, 20), color = (:skyblue, 0.5))

# 3. Universal Safe Zone (>= 180bp)
poly!(ax1, Rect(180, 0, 70, 45), color = (:lightgreen, 0.5))

# Add threshold lines
vlines!(ax1, [60, 180], color = :black, linestyle = :dash, linewidth = 1.5)
hlines!(ax1, [25], xmin=0.24, xmax=0.72, color = :black, linestyle = :dash, linewidth = 1.5)

# Annotate the zones
text!(ax1, 215, 22.5, text="SAFE\n(All Models)", align=(:center, :center), font=:bold)
text!(ax1, 120, 35, text="SAFE\n(High Evidence)", align=(:center, :center), font=:bold, color=:darkblue)
text!(ax1, 90, 12.5, text="DELETION\nVOID", align=(:center, :center), font=:bold, color=:darkred)
text!(ax1, 180, 2, text=" short_model_limit (180)", align=(:left, :bottom), fontsize=12)
text!(ax1, 60, 2, text=" abs_short (60)", align=(:left, :bottom), fontsize=12)
text!(ax1, 65, 25, text=" support-threshold (25.0)", align=(:left, :bottom), fontsize=12)


# ==============================================================================
# PANEL B: Hierarchical Evidence Collision Matrix
# ==============================================================================
ax2 = Axis(fig[1, 2], 
    title = "B. Hierarchical Evidence Resolution (Overlap > 120bp)",
    xlabel = "Gene B (Competitor)",
    ylabel = "Gene A (Anchor)",
    xticks = (1:3, ["rRNA / tRNA", "Best Alignment", "Ab Initio"]),
    yticks = (1:3, ["Ab Initio", "Best Alignment", "rRNA / tRNA"]),
    aspect = DataAspect()
)

# 3x3 Matrix representing Outcome for Gene A
# Values: 1 = Gene A Deleted (Red), 2 = Tie/Compare (Gray), 3 = Gene A Wins (Blue)
Z_plot = [
    2.0  1.0  1.0;  # X=1 (Gene B is rRNA)
    3.0  2.0  1.0;  # X=2 (Gene B is Align)
    3.0  3.0  2.0   # X=3 (Gene B is AbInitio)
]

heatmap!(ax2, 1:3, 1:3, Z_plot, colormap = [:lightcoral, :lightgray, :skyblue])

# Corrected Text Mapping (Reflecting the "Total Deletion" reality)
texts = [
    ("B Wins\n(A Deleted)", 1, 1), ("B Wins\n(A Deleted)", 2, 1), ("Compare\nScores", 3, 1),
    ("B Wins\n(A Deleted)", 1, 2), ("Compare\nScores", 2, 2),   ("A Wins\n(B Deleted)", 3, 2),
    ("Tie / Merge", 1, 3),           ("A Wins\n(B Deleted)", 2, 3), ("A Wins\n(B Deleted)", 3, 3)
]

for (t, x, y) in texts
    text!(ax2, x, y, text=t, align=(:center, :center), font=:bold)
end


# ==============================================================================
# PANEL C: The start_stop_allowance Optimization Track
# ==============================================================================
ax3 = Axis(fig[2, 1:2], 
    title = "C. Optimization via 'start_stop_allowance' (Resolving a 140bp Overlap)",
    xlabel = "Genomic Coordinate (bp)",
    limits = (0, 500, 0, 5)
)
hidespines!(ax3)
hidedecorations!(ax3, label=false)

# --- STATE 1: Original Conflict (Y = 3) ---
text!(ax3, 0, 4.2, text="1. Initial Conflict (Overlap = 140bp > 120bp Limit)", font=:bold)

poly!(ax3, Rect(0, 3, 300, 0.8), color = (:skyblue, 0.8))
poly!(ax3, Rect(160, 3, 300, 0.8), color = (:lightcoral, 0.8))
text!(ax3, 150, 3.4, text="Gene A", align=(:center, :center), font=:bold, color=:white)
text!(ax3, 310, 3.4, text="Gene B", align=(:center, :center), font=:bold, color=:white)

bracket!(ax3, 160, 2.8, 300, 2.8, text="140bp Overlap", offset=5, orientation=:down)

lines!(ax3, [240, 240, 300, 300, 240], [2.9, 3.9, 3.9, 2.9, 2.9], color=:black, linestyle=:dash)
lines!(ax3, [160, 160, 220, 220, 160], [2.9, 3.9, 3.9, 2.9, 2.9], color=:black, linestyle=:dash)
text!(ax3, 270, 4.0, text="60bp Allowance", align=(:center, :bottom), fontsize=12)
text!(ax3, 190, 4.0, text="60bp Allowance", align=(:center, :bottom), fontsize=12)

# --- STATE 2: Resolved State (Y = 1) ---
text!(ax3, 0, 2.2, text="2. Resolved State (Start/Stop motifs shifted within allowance)", font=:bold)

poly!(ax3, Rect(0, 1, 240, 0.8), color = (:dodgerblue, 0.9))
poly!(ax3, Rect(220, 1, 240, 0.8), color = (:firebrick, 0.9))

arrows!(ax3, [300], [2.8], [-60], [-0.8], color=:black, arrowsize=15)
arrows!(ax3, [160], [2.8], [60], [-0.8], color=:black, arrowsize=15)

bracket!(ax3, 220, 0.8, 240, 0.8, text="20bp Operon Overlap (< 120bp Limit)", offset=5, orientation=:down)

fig