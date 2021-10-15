function get_epsilon_effective(mym, mat, f_hole, epsilon_encl::T) where T<:AbstractFloat
    # calculate epsilon effective
    # hole is part with epsilon effective / encl is enclosure
    f_encl = collect(1:size(mym.elements2faces,1))
    filter!(x->x != f_hole, f_encl)
    temp_encl = 0
    temp_hole = 1000 # Kelvin
    temp = zeros(size(mym.elements,1),1)
    set_bc_face!(mym, temp, f_hole, temp_hole)
    set_bc_face!(mym, temp, f_encl, temp_encl)
    epsilon_hole = 1.0
    epsilon = zeros(size(mym.elements,1),1)
    set_bc_face!(mym, epsilon, f_hole, epsilon_hole)
    set_bc_face!(mym, epsilon, f_encl, epsilon_encl)
    # calculation of Qp
    Qp, Ga = tempsolver_old_algorithm_for_epseff(mym, mat, temp, epsilon)
    # area of hole
    e1 = mym.elements2faces[f_hole,3]
    e2 = mym.elements2faces[f_hole,4]
    area_hole = sum(mym.area[e1:e2,1])
    G_hole = sum(Ga[e1:e2,1])
    # calculation of epsilon effective
    epsilon_eff = 1.0 - (G_hole / (area_hole * sigma * (temp_hole^4 - temp_encl^4)))
    return epsilon_eff
end

function get_epsilon_effective(mym, mat, p_hole; epsilon_ink = 0.1)
    # calculate epsilon effective
    # hole is part with epsilon effective / encl is enclosure
    p_encl = collect(1:m.npar)
    filter!(x->x != p_hole, p_encl)
    temp_encl = 0
    temp_hole = 1000 # Kelvin
    temp = zeros(m.nelem,1)
    set_bc_part!(m, temp, p_hole, temp_hole)
    set_bc_part!(m, temp, p_encl, temp_encl)
    epsilon_hole = 1.0
    n_epsilon_encl = round(Integer,1.0/epsilon_ink)
    epsilon_encl = Matrix{Float64}(undef,n_epsilon_encl,2)
    epsilon_encl[:,1] = [epsilon_ink*i for i =1:n_epsilon_encl]
    epsilon = zeros(m.nelem,1)
    set_bc_part!(m, epsilon, p_hole, epsilon_hole)
    for i = 1:n_epsilon_encl
        set_bc_part!(m, epsilon, p_encl, epsilon_encl[i,1])
        # calculation of Qp
        Qp, Ga = tempsolver(m, mat, temp, epsilon)
        area_hole = get_area_of_part(m, p_hole)
        G_hole = sum(Ga[m.elem2par[p_hole].first:m.elem2par[p_hole].last])
        Qp_hole = sum(Qp[m.elem2par[p_hole].first:m.elem2par[p_hole].last])
        # calculation of epsilon effective
        epsilon_eff = 1.0 - (G_hole / (area_hole * sigma * (temp_hole^4 - temp_encl^4)))
        epsilon_encl[i,2] = epsilon_eff
    end
    return epsilon_encl
end