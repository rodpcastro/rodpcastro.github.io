def get_bem_solution(nl, xv, yv):
    """Using 4*nl elements, provides the solution for given points."""
    
    # Pre-processing.
    xb, yb, bt, bv = define_boundary(nl)
    xm, ym, lm, nx, ny = elements_properties(xb, yb)

    # Processing.
    u, q = solve_bem(xb, yb, bt, bv, xm, ym, lm, nx, ny)

    # Post-processing.
    s = get_domain_values(u, q, xv, yv, xb, yb, nx, ny, lm)

    return s
