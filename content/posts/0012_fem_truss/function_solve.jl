function solve(truss::Truss, fixed_nodes::Vector{Int64}, forces::Dict{Int64, Vector{Float64}})
    nds = truss.number_of_nodes
    nel = truss.number_of_elements

    # Displacement vector.
    d = zeros(Float64, 2*nds)

    # Force + reaction vector.
    f = zeros(Float64, 2*nds)

    # Stress vector.
    σ = zeros(Float64, nel)

    # Masks for essential and free nodes.
    en = fill(false, 2*nds)  # Essential nodes
    for i in fixed_nodes
        en[2*i-1:2*i] .= true
    end
    fn = .!en  # Free nodes

    # Stiffness submatrices.
    KEE = truss.K[en, en]
    KEF = truss.K[en, fn]
    KFE = truss.K[fn, en]
    KFF = truss.K[fn, fn]

    # Force vector at free nodes.
    for i in keys(forces)
        f[2*i-1:2*i] = forces[i]
    end
    fF = f[fn]  # Forces at free nodes

    # Displacement of free nodes.
    dF = KFF\fF
    d[fn] = dF
    
    # Reaction of fixed nodes.
    rE = KEF * dF
    f[en] = rE

    # Stress vector.
    for (i, e) in enumerate(truss.elements)
        i₁ = truss.connectivity[i, 1]
        i₂ = truss.connectivity[i, 2]
        
        de = zeros(Float64, 4)
        de[1:2] = d[2*i₁-1:2*i₁]
        de[3:4] = d[2*i₂-1:2*i₂]

        σ[i] = e.E/e.length * sum(e.R .* de)
    end

    # Reshaped displacement and force vectors
    dr = reshape(d, 2, nds)'
    fr = reshape(f, 2, nds)'
    
    return dr, fr, σ
end
