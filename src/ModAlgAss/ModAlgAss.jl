type ModAlgAss{S, T}
  base_ring::S
  action::Vector{T}
  dimension::Int
  isirreducible::Int
  dimension_splitting_field::Int
  
  function ModAlgAss{S, T}(action::Vector{T}) where {S, T}
    z = new{S, T}()
    z.action = action
    z.dimension = cols(action[1])
    z.base_ring = base_ring(action[1])
    if z.dimension == 1
      z.isirreducible = 1
      z.dimension_splitting_field = 1
    else 
      z.isirreducible = 0
      z.dimension_splitting_field = 0
    end
    return z
  end
end

function ModAlgAss(action::Vector{T}) where {T}
  @assert length(action) > 0
  S = typeof(base_ring(action[1]))
  return ModAlgAss{S, T}(action)
end

function isirreducible_known(M::ModAlgAss)
  return M.isirreducible != 0
end

function isirreducible(M::ModAlgAss)
  if M.isirreducible != 0
    return M.isirreducible == 1
  else
    error("Not implemented yet")
  end
end

function dimension(M::ModAlgAss)
  return M.dimension
end

function base_ring(M::ModAlgAss)
  return M.base_ring
end

# We a basis for the morphisms M -> N, that is,
# xTM.i = xN.iT for all i
function _homspace(M::ModAlgAss{S, T}, N::ModAlgAss{S, T}) where {S, T}
  base_ring(M) != base_ring(N) && error("Base rings must coincide")
  length(M.action) != length(N.action) && error("Number of action generators must be equal")
  if dimension(M) != dimension(N)
    return false
  end
  d = dimension(M)
  R = base_ring(M)
  r = length(M.action)
  X = zero_matrix(R, d^2 * r, d^2)
  for s in 1:r
    MM = M.action[s]
    NN = N.action[s]
    for k in 1:d
      for l in 1:d
        row_index = d^2 * (s - 1) + d * (k - 1) + l
        for j in 1:d
          X[row_index, (k - 1)* d + j] += MM[j, l]
          X[row_index, (j - 1)* d + l] -= NN[k, j]
        end
      end
    end
  end
  l, NS = nullspace(X)
  return NS
end

function homspace(M::ModAlgAss{S, T}, N::ModAlgAss{S, T}) where {S, T}
  d = dimension(M)
  R = base_ring(M)
  NS = _homspace(M, N)
  l = cols(NS)
  res = Vector{T}(l)
  for s in 1:l
    A = zero_matrix(R, d, d)
    for i in 1:d
      for j in 1:d
        A[i, j] = NS[(i - 1) * d + j, s]
      end
    end
    res[s] = A
  end
  return res
end

################################################################################
#
#  Lattices
#
################################################################################

type FLattice
  mod::ModAlgAss{FlintRationalField, fmpq_mat}
  basis_mat::FakeFmpqMat
  basis_mat_inv::FakeFmpqMat
  induced_action::Vector{fmpz_mat}

  function FLattice(mod::ModAlgAss{FlintRationalField, fmpq_mat}, basis_mat::FakeFmpqMat)
    return new(mod, basis_mat, inv(basis_mat))
  end

  function FLattice(mod::ModAlgAss{FlintRationalField, fmpq_mat}, basis_mat::fmpz_mat)
    return new(mod, FakeFmpqMat(basis_mat), inv(FakeFmpqMat(basis_mat)))
  end
end

function FLattice(mod::ModAlgAss{FlintRationalField, fmpq_mat})
  if all(isone(denominator(M)) for M in mod.action)
    return FLattice(mod, identity_matrix(FlintZZ, dimension(mod)))
  else
    error("No canonical lattice, since action is not integral")
  end
end

rank(L::FLattice) = dimension(L.mod)

function induced_action(L::FLattice)
  if isdefined(L, :induced_action)
    return L.induced_action
  else
    L.induced_action = fmpz_mat[to_fmpz_mat(L.basis_mat * FakeFmpqMat(A) * L.basis_mat_inv) for A in L.mod.action]
    return L.induced_action
  end
end

function _homspace(M::FLattice, N::FLattice)
  H = homspace(M.mod, N.mod)
  m = rank(M)
  n = rank(N)
  # We need to get the matrices with respect to the bases of L and M
  HLtoM = [ M.basis_mat * FakeFmpqMat(T) * N.basis_mat_inv for T in H ]
  X = zero_matrix(FlintZZ, length(HLtoM), m * n)
  for t in 1:length(HLtoM)
    T = HLtoM[t]
    for i in 1:m
      for j in 1:n
        X[t, (i - 1) * m + j] = T.num[i, j]
      end
    end
  end
  S = saturate(X)
  XX = Vector{fmpz_mat}(length(HLtoM))
  for t in 1:length(XX)
    A = zero_matrix(FlintZZ, m, n)
    for i in 1:m
      for j in 1:n
        A[i, j] = X[t, (i - 1) * m + j]
      end
    end
    XX[t] = A
  end
  return XX
end

function change_ring(L::FLattice, R::Ring)
  M = induced_action(L)
  return ModAlgAss([ change_ring(T, R) for T in M ])
end

function change_ring(M::MatElem, R::Ring)
  N = zero_matrix(R, rows(M), cols(M))
  for i in 1:rows(M)
    for j in 1:cols(M)
      N[i, j] = R(M[i, j])
    end
  end
  return N
end
