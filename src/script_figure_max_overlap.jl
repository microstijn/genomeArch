using CairoMakie

# 1. Setup data for the negative region (μ < 0)
# We avoid 0 exactly to prevent division by zero errors.
mu_neg = range(-0.1, stop=-0.001, length=1000)
L_limit_neg = sqrt.(-5 ./ mu_neg)

# 2. Create the Figure and Axis
fig = Figure(
    size = (300, 300),
    font = "DejaVu Sans",
    fontsize = 16
    )
ax = Axis(
    fig[1, 1],
    #title = L"Region satisfying $L_{max} \geq L_{actual}$ for $\mu \in [-0.1, 0.1]$ where $L_{max} = \frac{-10}{2 \mu L_{actual}}$",
    xlabel = L"\mu",
    xlabelsize = 18,
    ylabel = L"L_{actual}",
    ylabelsize = 24,
    limits = ((-0.1, 0.1), (0, 80))
) 

# 3. Plot the Boundary Line (Red Curve)
lines!(ax, mu_neg, L_limit_neg, color = :red, linewidth = 2, label = L"L_{actual} = \sqrt{-5 \mu}")

# 4. Shade the Valid Region (μ < 0)
# band! creates a filled area between two values (0 and the limit)
band!(ax, mu_neg, 0, L_limit_neg, color = (:skyblue, 0.4), label = L"\text{Satisfies}: L_{max} \geq L_{actual}")

# 5. Represent the No-Solution Zone (μ > 0)
# vspan! colors a vertical strip of the graph
vspan!(ax, 0, 0.1, color = (:gray, 0.2), label = L"No solution ($\mu > 0$)")

# 6. Add Reference Lines
vlines!(ax, [0], color = :black, linestyle = :dash, alpha = 0.5) # Zero line

# 7. Add Legend and Grid
axislegend(
    ax, position = :rt
)

ax.xgridstyle = :dash
ax.ygridstyle = :dash

# 8. Display or Save
save("my_graph.png", fig)
fig