nodes = Float64[
    0.0 0.0;
    1.0 0.0;
    2.0 0.0;
    3.0 0.0;
    1.0 1.0;
    2.0 1.0;
    3.0 1.0;
]

connectivity = Int64[
    1 2;
    1 5;
    2 3;
    2 5;
    2 6;
    3 4;
    3 6;
    3 7;
    4 7;
    5 6;
    6 7;
]

# Circular cross section.
D = 0.02     # Diameter (m)
A = π*D^2/4  # Area (m²)

# Young's modulus (Pa).
E = 205e9  # Steel

# Truss definition.
truss = Truss(nodes, connectivity, A, E);

# Boundary conditions.
fixed_nodes = [1, 4]
forces = Dict(2 => [0.0, -500.0], 3 => [0.0, -1000.0]);

# Solution.
d, f, σ = solve(truss, fixed_nodes, forces)
T = σ * A;
