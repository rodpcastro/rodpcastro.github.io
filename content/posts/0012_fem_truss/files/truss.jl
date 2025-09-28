using LinearAlgebra


struct Element
    node1::Vector{Float64}
    node2::Vector{Float64}
    A::Float64
    E::Float64
    length::Float64
    R::Vector{Float64}
    K::Matrix{Float64}

    function Element(node1, node2, A, E)
        Δ = node2 - node1
        length = √sum(Δ.^2)
        sinϕ = Δ[2]/length
        cosϕ = Δ[1]/length
        s2 = sinϕ^2
        c2 = cosϕ^2
        sc = sinϕ*cosϕ

        R = Float64[-cosϕ, -sinϕ, cosϕ, sinϕ]
        
        K = A*E/length * Float64[
             c2  sc -c2 -sc;
             sc  s2 -sc -s2;
            -c2 -sc  c2  sc;
            -sc -s2  sc  s2;
        ]
        
        new(node1, node2, A, E, length, R, K)
    end
end


struct Truss
    nodes::Matrix{Float64}
    connectivity::Matrix{Int64}
    A::Float64
    E::Float64
    number_of_elements::Int64
    number_of_nodes::Int64
    elements::Vector{Element}
    K::Matrix{Float64}

    function Truss(nodes, connectivity, A, E)
        nel = size(connectivity, 1)
        nds = size(nodes, 1)
        elements = Array{Element}(undef, nel)
        K = zeros(Float64, 2*nds, 2*nds)
        for e in 1:nel
            i₁ = connectivity[e, 1]
            i₂ = connectivity[e, 2]

            element = Element(nodes[i₁,:], nodes[i₂,:], A, E)
            elements[e] = element

            L = zeros(Int8, 4, 2*nds)
            L[1:2, 2*i₁-1:2*i₁] = Matrix{Int8}(I, 2, 2)
            L[3:4, 2*i₂-1:2*i₂] = Matrix{Int8}(I, 2, 2)
            K += L' * element.K * L
        end
        new(nodes, connectivity, A, E, nel, nds, elements, K)
    end
end


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