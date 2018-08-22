###############################################################################
#
#  NfRelOrd
#
###############################################################################

# I don't like that we have to drag S around
# Still hoping for julia/#18466

mutable struct NfRelOrdSet{T}
  nf::RelativeExtension{T}

  function NfRelOrdSet{T}(K::RelativeExtension{T}) where {T}
    a = new(K)
    return a
  end
end

mutable struct NfRelOrd{T, S} <: Ring
  nf::RelativeExtension{T}
  basis_nf::Vector{RelativeElement{T}}
  basis_mat::Generic.Mat{T}
  basis_mat_inv::Generic.Mat{T}
  basis_pmat::PMat{T, S}
  pseudo_basis::Vector{Tuple{RelativeElement{T}, S}}

  disc_abs::NfOrdIdl # used if T == nf_elem
  disc_rel#::NfRelOrdIdl{T} # used otherwise; is a forward declaration
  parent::NfRelOrdSet{T}

  isequation_order::Bool

  ismaximal::Int                   # 0 Not known
                                   # 1 Known to be maximal
                                   # 2 Known to not be maximal

  trace_mat::Generic.Mat{T}

  inv_coeff_ideals::Vector{S}

  function NfRelOrd{T, S}(K::RelativeExtension{T}) where {T, S}
    z = new{T, S}()
    z.nf = K
    z.parent = NfRelOrdSet{T}(K)
    z.isequation_order = false
    z.ismaximal = 0
    return z
  end

  function NfRelOrd{T, S}(K::RelativeExtension{T}, M::PMat{T, S}) where {T, S}
    z = NfRelOrd{T, S}(K)
    z.nf = K
    z.parent = NfRelOrdSet{T}(K)
    z.basis_pmat = M
    z.basis_mat = M.matrix
    return z
  end

  function NfRelOrd{T, S}(K::RelativeExtension{T}, M::Generic.Mat{T}) where {T, S}
    z = NfRelOrd{T, S}(K)
    z.nf = K
    z.parent = NfRelOrdSet{T}(K)
    z.basis_mat = M
    z.basis_pmat = pseudo_matrix(M)
    return z
  end
end

###############################################################################
#
#  NfRelOrdElem
#
###############################################################################

mutable struct NfRelOrdElem{T} <: RingElem
  parent#::NfRelOrd{T, S} # I don't want to drag the S around
  elem_in_nf::RelativeElement{T}
  elem_in_basis::Vector{T}
  has_coord::Bool

  function NfRelOrdElem{T}(O::NfRelOrd{T}) where {T}
    z = new{T}()
    z.parent = O
    z.elem_in_nf = zero(nf(O))
    z.elem_in_basis = Vector{T}(undef, degree(O))
    z.has_coord = false
    return z
  end

  function NfRelOrdElem{T}(O::NfRelOrd{T}, a::RelativeElement{T}) where {T}
    z = new{T}()
    z.parent = O
    z.elem_in_nf = a
    z.elem_in_basis = Vector{T}(undef, degree(O))
    z.has_coord = false
    return z
  end

  function NfRelOrdElem{T}(O::NfRelOrd{T}, a::RelativeElement{T}, arr::Vector{T}) where {T}
    z = new{T}()
    z.parent = O
    z.elem_in_nf = a
    z.elem_in_basis = arr
    z.has_coord = true
    return z
  end
end

###############################################################################
#
#  NfRelOrdFracIdl
#
###############################################################################

mutable struct NfRelOrdFracIdlSet{T, S}
  order::NfRelOrd{T, S}

  function NfRelOrdFracIdlSet{T, S}(O::NfRelOrd{T, S}) where {T, S}
    a = new(O)
    return a
  end
end

mutable struct NfRelOrdFracIdl{T, S}
  order::NfRelOrd{T, S}
  parent::NfRelOrdFracIdlSet{T, S}
  basis_pmat::PMat{T, S}
  pseudo_basis::Vector{Tuple{RelativeElement{T}, S}}
  basis_mat::Generic.Mat{T}
  basis_mat_inv::Generic.Mat{T}
  den::fmpz

  norm
  has_norm::Bool

  function NfRelOrdFracIdl{T, S}(O::NfRelOrd{T, S}) where {T, S}
    z = new{T, S}()
    z.order = O
    z.parent = NfRelOrdFracIdlSet{T, S}(O)
    z.has_norm = false
    return z
  end

  function NfRelOrdFracIdl{T, S}(O::NfRelOrd{T, S}, M::PMat{T, S}) where {T, S}
    z = NfRelOrdFracIdl{T, S}(O)
    z.basis_pmat = M
    z.basis_mat = M.matrix
    return z
  end

  function NfRelOrdFracIdl{T, S}(O::NfRelOrd{T, S}, a::Array{Tuple{T1, S}}) where {T1 <: RelativeElement{T}, S} where T
    z = NfRelOrdFracIdl{T, S}(O)
    z.pseudo_basis = a
    return z
  end
end

###############################################################################
#
#  NfRelOrdIdl
#
###############################################################################

mutable struct NfRelOrdIdlSet{T, S}
  order::NfRelOrd{T, S}

  function NfRelOrdIdlSet{T, S}(O::NfRelOrd{T, S}) where {T, S}
    a = new(O)
    return a
  end
end

mutable struct NfRelOrdIdl{T, S}
  order::NfRelOrd{T, S}
  parent::NfRelOrdIdlSet{T, S}
  basis_pmat::PMat{T, S}
  pseudo_basis::Vector{Tuple{RelativeElement{T}, S}}
  basis_mat::Generic.Mat{T}
  basis_mat_inv::Generic.Mat{T}

  norm
  has_norm::Bool
  is_prime::Int            # 0: don't know
                           # 1 known to be prime
                           # 2 known to be not prime
  splitting_type::Tuple{Int, Int}

  minimum
  non_index_div_poly::fq_poly # only used if the ideal is a prime ideal not dividing the index
  p_uniformizer::NfRelOrdElem{T}
  anti_uniformizer::RelativeElement{T}

  function NfRelOrdIdl{T, S}(O::NfRelOrd{T, S}) where {T, S}
    z = new{T, S}()
    z.order = O
    z.parent = NfRelOrdIdlSet{T, S}(O)
    z.has_norm = false
    z.is_prime = 0
    z.splitting_type = (0,0)
    return z
  end

  function NfRelOrdIdl{T, S}(O::NfRelOrd{T, S}, M::PMat{T, S}) where {T, S}
    z = NfRelOrdIdl{T, S}(O)
    z.basis_pmat = M
    z.basis_mat = M.matrix
    return z
  end

  function NfRelOrdIdl{T, S}(O::NfRelOrd{T, S}, a::Array{Tuple{T1, S}}) where {T1 <: RelativeElement{T}, S} where T
    z = NfRelOrdIdl{T, S}(O)
    z.pseudo_basis = a
    return z
  end
end
