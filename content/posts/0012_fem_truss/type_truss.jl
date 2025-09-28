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
