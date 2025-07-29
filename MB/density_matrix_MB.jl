##this is the program for ploting the density matrix after one input photon in a MB theory

using Plots

# Parameters
λ = 1.0                     # Wavelength
nat = 1.0 / λ^3            # Atomic density
rb = 6.0 * λ               # Blockade radius


# Grid
x_vals = range(0, stop=20.0, length=500)
y_vals = range(0, stop=20.0, length=500)
Z = zeros(length(x_vals), length(y_vals))

# Function definition
function coherence(x, y, λ, nat, rb)
    δ = abs(x - y)
    if x > rb && y > rb
        L=exp(- (3/2π) * δ / λ * nat * λ^3)
    elseif x < rb && y < rb
        L=1
    elseif x <= rb
        #L=exp(- (3/(4π)) * δ / λ * nat * λ^3 - (3/(4π)) * (y - rb) / λ * nat * λ^3)
        L=exp(- (3/(2π)) * (y - rb) / λ * nat * λ^3)
    else # y <= rb
        #L=exp(- (3/(4π)) * δ / λ * nat * λ^3- (3/(4π)) * (x - rb) / λ * nat * λ^3)
        L=exp(- (3/(2π)) * (x - rb) / λ * nat * λ^3)
    end
    return L=exp(-1) *(1+L+1/2*L^2+1/6*L^3)
end


# Compute values
for (i, x) in enumerate(x_vals)
    for (j, y) in enumerate(y_vals)
        Z[i, j] = coherence(x, y, λ, nat, rb)
    end
end

# Plot
heatmap(x_vals, y_vals, Z',
    xlabel="x",
    ylabel="y",
    colorbar_title="Coherence",
    title="Spatial Coherence Map",
    aspect_ratio=:equal, xlims=(0,20), clims=(0.4, 1.0),
    c=:inferno)