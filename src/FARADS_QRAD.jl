module FARADS_QRAD

    using LinearAlgebra

    using FARADS_GEOM
    using FARADS_MESHING
    
    include("./tempsolver_fct.jl")
    export tempsolver_old_algorithm
    export sigma

    export tempsolver_gebhart_getmatrix
    export tempsolver_gebhart_calcQ

    export tempsolver_modest_getmatrix
    export tempsolver_modest_calcQ

    export tempsolver_modest_getmatrix_opt
    export tempsolver_modest_calcQ_opt

    export tempsolver_old_algorithm_opt


    include("./boundary_fct.jl")
    export set_bc_face!
    export set_bc_part!

    export set_bc_face_by_com_linear!
    export set_bc_face_by_com_linear2!

    include("./epsilon_effective_fct.jl")
    export get_epsilon_effective

end