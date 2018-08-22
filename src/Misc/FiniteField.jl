##
## rand for Flint-Finite fields
##
## fq_nmod has no base ring
function rand(R::FqNmodFiniteField)
  #gen is not a generator for the group!
  r = R(0)
  for i=0:degree(R)
    r *= gen(R)
    r += rand(1:characteristic(R))
  end

  return r
end

function rand(R::FqFiniteField)
  #gen is not a generator for the group!
  r = R(0)
  for i=0:degree(R)
    r *= gen(R)
    r += rand(1:characteristic(R))
  end

  return r
end

################################################################################
#
#  Iterators for finite fields
#
################################################################################

# FqNmodFiniteField

function Base.iterate(F::FqNmodFiniteField)
  return zero(F), zeros(UInt, degree(F))
end

function Base.iterate(F::FqNmodFiniteField, st::Vector{UInt})
  done = true
  for j in 1:length(st)
    if st[j] != F.n - 1
      done = false
    end
  end

  if done
    return nothing
  end

  st[1] = st[1] + 1
  for j in 1:(length(st) - 1)
    if st[j] == F.n
      st[j] = 0
      st[j + 1] = st[j + 1] + 1
    end
  end
  
  d = F()
  ccall((:fq_nmod_init2, :libflint), Nothing, 
        (Ref{fq_nmod}, Ref{FqNmodFiniteField}), d, F)
  for j in 1:length(st)
         ccall((:nmod_poly_set_coeff_ui, :libflint), Nothing,
              (Ref{fq_nmod}, Int, UInt), d, j - 1, st[j])
  end

  return d, st
end

Base.length(F::FqNmodFiniteField) = BigInt(F.n)^degree(F)

Base.IteratorSize(::Type{FqNmodFiniteField}) = Base.HasLength()

Base.eltype(::FqNmodFiniteField) = fq_nmod

# FqFiniteField

function Base.iterate(F::FqFiniteField)
  return zero(F), zeros(FlintZZ, degree(F))
end

function Base.iterate(F::FqFiniteField, st::Vector{fmpz})
  done = true
  for j in 1:length(st)
    if st[j] != characteristic(F) - 1
      done = false
    end
  end

  if done
    return nothing
  end

  st[1] = st[1] + 1
  for j in 1:(length(st) - 1)
    if st[j] == characteristic(F)
      st[j] = 0
      st[j + 1] = st[j + 1] + 1
    end
  end

  d = F()
  ccall((:fq_init2, :libflint), Nothing, 
        (Ref{fq}, Ref{FqFiniteField}), d, F)
  for j in 1:length(st)
         ccall((:fmpz_mod_poly_set_coeff_fmpz, :libflint), Nothing,
               (Ref{fq}, Int, Ref{fmpz}), d, j - 1, st[j])
  end

  return d, st
end

Base.eltype(::FqFiniteField) = fq

Base.length(F::FqFiniteField) = BigInt(characteristic(F)^degree(F))

Base.IteratorSize(::Type{FqFiniteField}) = Base.HasLength()

sub!(z::T, x::T, y::T) where {T} = x - y

function _reduce(a::fq_nmod)
  A = parent(a)
  if a.length < 2*degree(A)
    ccall((:fq_nmod_reduce, :libflint), Nothing, (Ref{fq_nmod}, Ref{FqNmodFiniteField}), a, A)
  else
    ccall((:nmod_poly_rem, :libflint), Nothing, (Ref{fq_nmod}, Ref{fq_nmod}, Ref{Nothing}, Ref{Nothing}), a, a, pointer_from_objref(A)+6*sizeof(Int) + 2*sizeof(Ptr{Nothing}), pointer_from_objref(A)+sizeof(fmpz))
  end
end

function (R::FqFiniteField)(x::fmpz_mod_poly)
  z = R()
  ccall((:fq_set, :libflint), Nothing, (Ref{Nemo.fq}, Ref{Nemo.fmpz_mod_poly}, Ref{Nemo.FqFiniteField}), z, x, R)
  ccall((:fq_reduce, :libflint), Nothing, (Ref{Nemo.fq}, Ref{Nemo.FqFiniteField}), z, R)
  return z
end

#TODO: move elsewhere - and use. There are more calls to nmod_set/reduce
function (A::FqNmodFiniteField)(x::nmod_poly)
  u = A()
  ccall((:fq_nmod_set, :libflint), Nothing,
                     (Ref{fq_nmod}, Ref{nmod_poly}, Ref{FqNmodFiniteField}),
                                     u, x, A)
  _reduce(u)
#  ccall((:fq_nmod_reduce, :libflint), Nothing, (Ref{fq_nmod}, Ref{FqNmodFiniteField}), u, A)                                   
  return u
end

function _nf_to_fq!(a::fq_nmod, b::nf_elem, K::FqNmodFiniteField, a_tmp::nmod_poly)
  nf_elem_to_nmod_poly!(a_tmp, b)
  ccall((:fq_nmod_set, :libflint), Nothing,
                     (Ref{fq_nmod}, Ref{nmod_poly}, Ref{FqNmodFiniteField}),
                                     a, a_tmp, K)
  _reduce(a)
#  ccall((:fq_nmod_reduce, :libflint), Nothing, (Ref{fq_nmod}, Ref{FqNmodFiniteField}), a, K)                                   
end
  

