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

function tempsolver_gebhart_getmatrix(m, vfmat, epsilon)
    # Gebhart algorithm optimized
    # aus Clark1974
    # does not work with epsilon = 1.0
    # klappt
    # TO DO: verbesserung der inversenbildung
    n_elements = size(m.elements,1)
    mat = zeros(n_elements, n_elements)
    # identity = zeros(n_elements, n_elements)
    # for i = 1:n_elements
    #     identity[i,i] = 1
    # end
    # zusatz = identity .* (1 ./ (1 .- epsilon))
    # mat = vfmat .- zusatz
    mat = vfmat .- Matrix(1I,n_elements,n_elements) .* (1 ./ (1 .- epsilon))
    mat .= inv(mat)  
    mat .= mat .* (-1)
    mat .= mat' * vfmat'
    for i = 1:n_elements
        mat[i,:] = mat[i,:] ./ (m.area[i,1] ./ m.area[:,1] .* ((1 .- epsilon[i,1]) ./ epsilon))
    end
    return mat
end

# function tempsolver_gebhart_getmatrix2(m, vfmat, epsilon)
#     # Gebhart algorithm optimized
#     # aus Bornside1990
#     # klappt NICHT
#     n_elements = size(m.elements,1)
#     mat = zeros(n_elements, n_elements)
#     for j = 1:n_elements
#         mat[:,j] = (epsilon[j,1] .- 1) .* vfmat[:,j]
#     end
#     mat .= mat .- Matrix{Int}(I, n_elements, n_elements)
#     gebhart = zeros(n_elements, n_elements)
#     G = zeros(n_elements, 1)
#     for j = 1:n_elements
#         b = epsilon[j,1] .* vfmat[:,j]
#         b .= b .* (-1)
#         G = mat \ b
#         gebhart[:,j] = G
#     end
#     return gebhart
# end

function tempsolver_gebhart_calcQ(m, gebhart, temp, epsilon, rounds_max)
    # Gebhart algorithm optimized
    # aus Clark1974
    # does not work with epsilon = 1.0
    n_elements = size(m.elements,1)
    Qp = zeros(n_elements,1)
    Qsum = zeros(n_elements,1)
    for rounds = 1:rounds_max
        for j = 1:n_elements
            Qsum[j,1] = sum(epsilon[:,1] .* m.area[:,1] .* gebhart[:,j] .* sigma .* temp[:,1].^4)
        end
        Qp[:,1] = (epsilon[:,1] .* m.area[:,1] .* sigma .* temp[:,1].^4) - Qsum[:,1]
        # temp .= temp .* 1.1
    end
    return Qp
end

function tempsolver_modest_getmatrix(m, vfmat, epsilon)
    # Modest2013 algorithm S.586ff
    # klappt
    n_elements = size(m.elements,1)
    ss = zeros(n_elements, n_elements)
    T = zeros(n_elements, n_elements)
    S = zeros(n_elements, n_elements)
    SS = zeros(n_elements, n_elements)
    for i = 1:n_elements
        ss[i,:] = vfmat[i,:] .* m.area[i,:]
    end
    for j = 1:n_elements
        S[:,j] = ss[:,j] .* epsilon[j,:]
    end
    for j = 1:n_elements
        T[:,j] = (1 .- epsilon[j,:]) .* ss[:,j] ./ (epsilon[j,:] .* m.area[j,1])
    end
    zusatz = 1 ./ epsilon .* Matrix(1I,n_elements,n_elements)
    T = zusatz .- T
    Tinv = inv(T)
    SS = Tinv * S
    return SS
end

function tempsolver_modest_calcQ(m, SS, temp, rounds_max)
    # Modest2013 algorithm S.586ff
    n_elements = size(m.elements,1)
    Qp = zeros(n_elements,1)
    for rounds = 1:rounds_max
        for i = 1:n_elements
            Qsum = 0
            for j = 1:n_elements
                Qsum = Qsum + SS[i,j] * sigma * (temp[i,1]^4 - temp[j,1]^4)
            end
            Qp[i,1] = Qsum
        end
        # temp .= temp .* 1.1
    end
    return Qp
end

function tempsolver_modest_getmatrix_opt(m, mat, epsilon)
    # Modest2013 algorithm S.586ff
    # optimiert nach eigenem ermessen
    # klappt
    n_elements = size(m.elements,1)
    mat2 = zeros(n_elements, n_elements)
    for i = 1:n_elements
        mat[i,:] .= mat[i,:] .* m.area[i,:]
    end
    for j = 1:n_elements
        mat[:,j] .= mat[:,j] .* epsilon[j,:]
    end
    for j = 1:n_elements
        mat2[:,j] = (1 .- epsilon[j,:]) .* (mat[:,j] ./ epsilon[j,:]) ./ (epsilon[j,:] .* m.area[j,1])
    end
    mat2 .= (1 ./ epsilon .* Matrix(1I,n_elements,n_elements)) .- mat2
    mat2 .= inv(mat2)
    mat2 .= mat2 * mat
    return mat2
end

function tempsolver_modest_calcQ_opt(m, SS, temp, rounds_max)
    # Modest2013 algorithm S.586ff
    n_elements = size(m.elements,1)
    Qp = zeros(n_elements,1)
    for rounds = 1:rounds_max
        for i = 1:n_elements
            Qp[i,1] = sum(SS[i,:] .* sigma .* (temp[i,1]^4 .- temp[:,1].^4))
        end
        # temp .= temp .* 1.1
    end
    return Qp
end

function tempsolver_old_algorithm_opt(m, vfmat, temp, epsilon, rounds_max)
    # old algorithm heavily optimized 
    # verändert nach howell s.222 (ansatz nach J direkt q berechnen (ohne Fläche))
    # does not work with epsilon = 1.0
    # klappt
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