
function set_bc_face!(mym, bcvec, faces, bcface)
    # give all face elements one bc
    for f = faces
        e1 = mym.elements2faces[f,3]
        e2 = mym.elements2faces[f,4]
        bcvec[e1:e2,1] .= bcface
    end
end

function set_bc_part!(mym, bcvec, parts, bcpart)
    # give all part elements one bc
    for p = parts
        e1 = mym.elements2parts[p,3]
        e2 = mym.elements2parts[p,4]
        bcvec[e1:e2,1] .= bcpart
    end
end

function set_bc_face_by_com_linear!(mym, bcvec, faces, bc_min, bc_max, dim)
    # give face elements bc according to com

    bc_diff = bc_max - bc_min

    for f = faces

        e1 = mym.elements2faces[f,3]
        e2 = mym.elements2faces[f,4]

        min_XYZ = vec(minimum(mym.com[e1:e2,:], dims = 1))
        max_XYZ = vec(maximum(mym.com[e1:e2,:], dims = 1))
        dim_min = min_XYZ[dim]
        dim_max = max_XYZ[dim]
        # println("values dim=", dim, " between ", dim_min, " and ", dim_max)

        dim_diff = dim_max - dim_min

        for i = e1:e2
            bcvec[i,1] = bc_min + (bc_diff / dim_diff) * (mym.com[i,dim] - dim_min)
        end

    end
end


function set_bc_face_by_com_linear2!(mym, bcvec, faces, bc1, bc2, bc3, grenze, dim)
    # give face elements bc according to com

    bc_diff1 = bc2 - bc1
    bc_diff2 = bc3 - bc2

    for f = faces

        e1 = mym.elements2faces[f,3]
        e2 = mym.elements2faces[f,4]

        min_XYZ = vec(minimum(mym.com[e1:e2,:], dims = 1))
        max_XYZ = vec(maximum(mym.com[e1:e2,:], dims = 1))
        dim_min = min_XYZ[dim]
        dim_max = max_XYZ[dim]
        # println("values dim=", dim, " between ", dim_min, " and ", dim_max)

        dim_diff = dim_max - dim_min

        for i = e1:e2
            if (mym.com[i,dim] - dim_min) <= (grenze*dim_diff)
                bcvec[i,1] = bc1 + (bc_diff1 / (grenze*dim_diff)) * (mym.com[i,dim] - dim_min)
            else
                bcvec[i,1] = bc2 + (bc_diff2 / ((1-grenze)*dim_diff)) * (mym.com[i,dim] - dim_min - (grenze*dim_diff))
            end
        end

    end
end