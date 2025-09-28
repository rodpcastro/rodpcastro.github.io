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
