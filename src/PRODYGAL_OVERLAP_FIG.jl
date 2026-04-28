using CairoMakie
using Distributions

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

# 2. Define Biological Parameters
μ, σ, β = 0.8, 1.5, 2.5
p_c, p_rbs = 0.05, 0.40
p = p_c * p_rbs
α = 0.5 # Softplus smoothing factor

# --- UNIFIED INDICATOR STYLES ---
ind_line = :dash
ind_color = (:black, 0.6)
ind_fill = (:black, 0.08)
ind_font = 12

# 3. Define Phase Spaces (Expanded Panel A to show the cap)
L_genes = 90:5:1500              
L_over_conv = 1:1:220            
L_over_uni = 1:1:300             

# 4. Calculate Matrices
prob_conv = zeros(length(L_genes), length(L_over_conv))
prob_uni_fail = zeros(length(L_genes), length(L_over_uni))
exp_loss_uni = zeros(length(L_genes), length(L_over_uni))

# Populate Convergent Deletion Probability
for (i, L) in enumerate(L_genes)
    dist = Normal(μ * L, σ * sqrt(L))
    for (j, Lo) in enumerate(L_over_conv)
        if Lo > 200
            prob_conv[i, j] = NaN # Stop contours in the forbidden zone
        else
            prob_conv[i, j] = cdf(dist, β * Lo)
        end
    end
end

# Populate Unidirectional Matrices 
for (i, L) in enumerate(L_genes)
    S_max = max(0, L - 90)
    for (j, Lo) in enumerate(L_over_uni)
        
        S_min_strict = max(0, Lo - 60)
        
        # 1. Check for the physical Deletion Void first
        if S_min_strict > S_max
            prob_uni_fail[i, j] = NaN 
            exp_loss_uni[i, j] = NaN
            continue
        end
        
        # 2. Panel B logic (Total Deletion Probability)
        if Lo <= 60
            prob_uni_fail[i, j] = 0.0  
        else
            k_strict = floor(Int, (S_max - S_min_strict) / 3) + 1
            prob_uni_fail[i, j] = (1 - p)^k_strict
        end
        
        # 3. Panel C logic (Expected Loss with Softplus Smooth)
        if α * (Lo - 60) > 50
            S_smooth = Lo - 60  
        else
            S_smooth = (1 / α) * log(1 + exp(α * (Lo - 60)))
        end
        
        k_smooth = max(0, floor(Int, (S_max - S_smooth) / 3))
        
        if k_smooth == 0 || (1 - (1 - p)^k_smooth) == 0
            E_penalty = 0.0
        else
            num = 1 - (k_smooth + 1)*(1 - p)^k_smooth + k_smooth*(1 - p)^(k_smooth + 1)
            den = p * (1 - (1 - p)^k_smooth)
            E_penalty = 3 * (num / den - 1)
        end
        
        P_trunc = cdf(Normal(40, 8), Lo)
        exp_loss_uni[i, j] = P_trunc * (S_smooth + E_penalty)
    end
end

# 5. Generate the 3-Panel Figure
fig = Figure(size = (900, 300))

# --- Panel A ---
ax1 = Axis(
    fig[1, 1],
    title = "Convergent overlap\n(probability of deletion)",
   # titlealign = :left,
    xlabel = "Gene Length (bp)",
    ylabel = "Overlap Length (bp)"
)
contour!(ax1, L_genes, L_over_conv, prob_conv, levels = 0.1:0.1:0.9, color = (:black, 0.3), linewidth = 1)
contour!(ax1, L_genes, L_over_conv, prob_conv, levels = [0.5], color = :darkred, linewidth = 2.5, linestyle = :solid)

# Panel A Geometries
#band!(ax1, L_genes, fill(200, length(L_genes)), fill(220, length(L_genes)), color=ind_fill)
hlines!(ax1, [200], color=ind_color, linestyle=ind_line, linewidth=1.5)
#text!(ax1, 800, 210, text="Deletion zone @ 200bp", color=ind_color, fontsize=ind_font, align=(:center, :center))
#text!(ax1, 1100, 160, text="50% Threshold", color=:darkred, fontsize=ind_font)

# --- Define the Diagonal Void Boundary for Panels B & C ---
# The exact mathematical limit where S_min > S_max simplifies to L_o > L - 30
void_boundary_y = max.(60, L_genes .- 30) 

# --- Panel B ---
ax2 = Axis(
    fig[1, 2],
    title = "Unidirectional overlap\n(Probability of deletion)",
     #   titlealign = :left,
    xlabel = "Original Gene Length (bp)",
    ylabel = "True Overlap Length (bp)"
)
contour!(ax2, L_genes, L_over_uni, prob_uni_fail, levels = 0.1:0.1:0.9, color = (:black, 0.3), linewidth = 1)
contour!(ax2, L_genes, L_over_uni, prob_uni_fail, levels = [0.5], color = :darkred, linewidth = 2.5, linestyle = :solid)

# Panel B Geometries
band!(ax2, L_genes, void_boundary_y, fill(300, length(L_genes)), color=ind_fill)
lines!(ax2, L_genes, void_boundary_y, color=ind_color, linestyle=ind_line, linewidth=1.5)
hlines!(ax2, [60], color=ind_color, linestyle=ind_line, linewidth=1.5)
#text!(ax2, 1000, 68, text="Mandatory Limit (60bp)", color=ind_color, fontsize=ind_font)
#text!(ax2, 350, 260, text="DELETION VOID\n(Gene strictly < 90bp)", color=ind_color, fontsize=14, align=(:center, :center))
#text!(ax2, 1100, 240, text="50% Threshold", color=:darkred, fontsize=ind_font)

# --- Panel C ---
ax3 = Axis(
    fig[1, 3],
    title = "Unidirectional overlap\n(Expected bp loss after truncation)",
    xlabel = "Original Gene Length (bp)",
    ylabel = "True Overlap Length (bp)"
)
contour!(ax3, L_genes, L_over_uni, exp_loss_uni, levels = 0:50:400, color = (:navy, 0.4), linewidth = 1, labels = true)
contour!(ax3, L_genes, L_over_uni, exp_loss_uni, levels = [200], color = :navy, linewidth = 2.5, linestyle = :solid)

# Panel C Geometries (Identical to B to homogenize)
band!(ax3, L_genes, void_boundary_y, fill(300, length(L_genes)), color=ind_fill)
lines!(ax3, L_genes, void_boundary_y, color=ind_color, linestyle=ind_line, linewidth=1.5)
hlines!(ax3, [60], color=ind_color, linestyle=ind_line, linewidth=1.5)

# Panel C Scheme Indicators
#text!(ax3, 1100, 70, text="Mandatory truncation zone", color=ind_color, fontsize=ind_font, align=(:center, :bottom))
#text!(ax3, 1100, 50, text="Competitive zone", color=ind_color, fontsize=ind_font, align=(:center, :top))
#text!(ax3, 350, 260, text="DELETION VOID\n(Truncation Impossible)", color=ind_color, fontsize=14, align=(:center, :center))
#text!(ax3, 1100, 115, text="200bp Loss", color=:navy, fontsize=ind_font)

for ax in [ax2, ax3]
    xlims!(ax, 90, 900)
    ylims!(ax, 0, 300)
end


for ax in [ax1]
    xlims!(ax, 90, 900)
end

Label(fig[1, 1, TopLeft()], "a", font=:bold, fontsize=20, halign=:right, padding=(0, 10, 5, 0))
Label(fig[1, 2, TopLeft()], "b", font=:bold, fontsize=20, halign=:right, padding=(0, 10, 5, 0))
Label(fig[1, 3, TopLeft()], "c", font=:bold, fontsize=20, halign=:right, padding=(0, 10, 5, 0))

fig