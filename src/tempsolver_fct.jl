const sigma = 5.670374419E-8

function tempsolver_old_algorithm(m, vfmat, temp, epsilon, rounds_max)
    # old algorithm used in matlab (slightly optimized)
    n_elements = size(m.elements,1)
    # init matrix
    mat = zeros(n_elements, n_elements)
    Qp = zeros(n_elements, 1) 
    # prepare matrix for LGS
    mat = ((epsilon .- 1) .* vfmat)
    mat .= mat .+ Matrix{Int}(I, n_elements, n_elements)
    for rounds = 1:rounds_max
        # solve LGS
        B = (temp[:,1].^4) .* sigma .* epsilon
        J = mat \ B
        # Bestimmung der einfallenden Strahlung G
        G = zeros(Float64,n_elements,1)
        for i = 1:n_elements
            G[i,1] = sum(J[:,1] .* vfmat[i, :])
        end
        # Qp Radiation per element
        Qp = (J[:,1] - G[:,1]) .* m.area[:,1]
        # temp .= temp .* 1.1
    end 
    return Qp
end

function tempsolver_old_algorithm_for_epseff(m, vfmat, temp, epsilon)
    # old algorithm used in matlab (slightly optimized)
    n_elements = size(m.elements,1)
    # init matrix
    mat = zeros(n_elements, n_elements)
    Qp = zeros(n_elements, 1) 
    # prepare matrix for LGS
    mat = ((epsilon .- 1) .* vfmat)
    mat .= mat .+ Matrix{Int}(I, n_elements, n_elements)
    # solve LGS
    B = (temp[:,1].^4) .* sigma .* epsilon
    J = mat \ B
    # Bestimmung der einfallenden Strahlung G
    G = zeros(Float64,n_elements,1)
    for i = 1:n_elements
        G[i,1] = sum(J[:,1] .* vfmat[i, :])
    end
    # Qp Radiation per element
    Qp = (J[:,1] - G[:,1]) .* m.area[:,1]
    Ga = G[:,1] .* m.area[:,1]
    return Qp, Ga
end

function tempsolver_old_algorithm_optimized(m, vfmat, temp, epsilon, rounds_max)
    # old algorithm heavily optimized 
    # verändert nach howell s.222 (ansatz nach J direkt q berechnen (ohne Fläche))
    # does not work with epsilon = 1.0
    n_elements = size(m.elements,1)
    # init matrix
    mat = zeros(n_elements, n_elements)
    Qp = zeros(n_elements, 1)
    q = zeros(n_elements, 1)
    # prepare matrix for LGS
    mat = ((epsilon .- 1) .* vfmat)
    mat .= mat .+ Matrix{Int}(I, n_elements, n_elements)
    # solve LGS
    for rounds = 1:rounds_max
        B = (temp[:,1].^4) .* sigma .* epsilon
        J = mat \ B
        # q Radiation per element
        q = epsilon ./ (1 .- epsilon) .* ((sigma .* temp[:,1].^4) - J[:,1])
        Qp = q .* m.area[:,1]
        # temp .= temp .* 1.1
    end
    return Qp
end

function tempsolver_modest_algorithm(m, vfmat, temp, epsilon, rounds_max)
    # Modest algorithm optimized lean rounds
    n_elements = size(m.elements,1)
    mata = vfmat .* m.area
    matb = ((1 .- epsilon) ./ (epsilon .* m.area[:,1]))' .* mata[:,:]
    mata .= mata .* epsilon
    zusatz = 1 ./ epsilon .* Matrix{Int}(I, n_elements, n_elements)
    matb .= zusatz - matb
    matb .= inv(matb)
    mata = matb * mata
    # Qp Rad
    Qp = zeros(n_elements,1)
    for rounds = 1:rounds_max
        for i = 1:n_elements
            Qp[i,1] = sum(mata[i,:] .* sigma .* (temp[i,1].^4 .- temp[:,1].^4))
        end
        # temp .= temp .* 1.1
    end
    return Qp
end

function tempsolver_gebhart_algorithm(m, vfmat, temp, epsilon, rounds_max)
    # Gebhart algorithm optimized
    # does not work with epsilon = 1.0
    # TO DO: verbesserung der inversenbildung
    n_elements = size(m.elements,1)
    mat = zeros(n_elements, n_elements)
    zusatz = Matrix{Int}(I, n_elements, n_elements) .* (1 ./ (1 .- epsilon))
    mat = vfmat .- zusatz
    mat .= inv(mat)  
    mat .= mat .* (-1)
    mat = mat' * vfmat'
    for i = 1:n_elements
        mat[i,:] = mat[i,:] ./ (m.area[i,1] ./ m.area[:,1] .* ((1 .- epsilon) ./ epsilon))
    end
    Qp = zeros(n_elements,1)
    for rounds = 1:rounds_max
        for j = 1:n_elements
            Qsum = sum(epsilon .* m.area[:,1] .* mat[:,j] .* sigma .* temp[:,1].^4)
            Qp[j,1] = (epsilon[j,1] .* m.area[j,1] .* sigma .* temp[j,1].^4) - Qsum
        end
        # temp .= temp .* 1.1
    end
    return Qp
end

