"""
  OrthonormalBasis{T}

Syntatically represents an orthonormal basis for ``\\mathbb{R}^2``. 
"""
Basis{T} = SMatrix{2,2,T,4}

"""
  orthonormal_basis(x)

Get an orthonormal basis from a choice of axis ``e_1``.
"""
function orthonormal_basis(x)
    x̂ = SVector{2}(normalize(x)...)
    ŷ = SVector(-x̂[2], x̂[1])
    return hcat(x̂, ŷ)
end

# helper function 
function _hcat_and_normalize(e...)
    return mapreduce(normalize, hcat, e)
end

"""
  change_of_basis_matrix(old_basis, new_basis)

Compute the change of basis matrix to go from `old_basis` to `new_basis`. 
"""
function change_of_basis_matrix(old_basis, new_basis)
    return inv(new_basis) * old_basis
end

change_of_basis_matrix(::UniformScaling, ::UniformScaling) = I

"""
  change_basis(p::Vec, v::Basis, w::Basis)

Change from coordinates `p` in `v` to new coordinates in `w`.
"""
change_basis(p, v, w) = change_of_basis_matrix(v, w) * p

"""
  change_basis(p::Vec, w::Basis)

Change from coordinates `p` in the standard basis to new coordinates in `w`. 
"""
change_basis(p, w) = change_basis(p, I, w)
