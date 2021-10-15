
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