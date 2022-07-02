module FARADS_QRAD

    using LinearAlgebra

    using FARADS_GEOM
    using FARADS_MESHING
    
    include("./tempsolver_fct.jl")
    export tempsolver_old_algorithm
    export sigma

    export tempsolver_gebhart_getmatrix
    export tempsolver_gebhart_calcQ


    include("./boundary_fct.jl")
    export set_bc_face!
    export set_bc_part!

    include("./epsilon_effective_fct.jl")
    export get_epsilon_effective

end