module MultDep

using Hecke
import Base.*

function multiplicative_group_mod_units_fac_elem(A::Array{nf_elem, 1}; use_max_ord::Bool = false)
  k = parent(A[1])
  @assert all(i->parent(i) === k, A)
  if use_max_ord
    zk = maximal_order(k)
    cp = coprime_base(A, zk)
  else
    cp = coprime_base(A)
  end
  sort!(cp, lt = (a,b) -> norm(a) > norm(b))
  M = sparse_matrix(FlintZZ)
  for a = A
    T = Tuple{Int, fmpz}[]
    for i = 1:length(cp)
      p = cp[i]
      v = valuation(a, p)
      if v != 0
        push!(T, (i, fmpz(v)))
      end
#      isone(I) && break
    end
    push!(M, sparse_row(FlintZZ, T))
  end
  h, t = Hecke.hnf_kannan_bachem(M, Val{true}, truncate = true)
  return h, t, cp
end

function units(h::SMat, t, b::Array{nf_elem, 1})
  u = FacElem{nf_elem, AnticNumberField}[]
  for i=nrows(h)+1:length(b)
    k = [fmpz(0) for i=b]
    k[i] = 1
    for i in length(t):-1:1
      Hecke.apply_right!(k, t[i])
    end
    push!(u, FacElem(b, k))
  end
  return u
end

function unit_group_mod_torsion_fac_elem(O::NfAbsOrd, u::Array{FacElem{nf_elem, AnticNumberField}, 1})
  U = Hecke._unit_group_init(O)
  s = signature(O)
  r = s[1] + s[2] - 1
  for y = u
    is_tors, p = istorsion_unit(y, false, U.tors_prec)
    U.tors_prec = max(p, U.tors_prec)
    if is_tors
      continue
    end
    Hecke._add_unit(U, y)
    if length(U.units) >= r
      break
    end
  end
  if length(U.units) < r
    # maybe use pAdic stuff here...
    error("not complete yet")
  end

  U.full_rank = true

  U.units = Hecke.reduce(U.units, U.tors_prec)

  for y = u
    is_tors, p = istorsion_unit(y, false, U.tors_prec)
    U.tors_prec = max(p, U.tors_prec)
    if is_tors
      continue
    end
    x = Hecke.reduce_mod_units([y], U)[1]
    is_tors, p = istorsion_unit(x, false, U.tors_prec)
    U.tors_prec = max(p, U.tors_prec)
    if is_tors
      continue
    end
    Hecke._add_dependent_unit(U, x)
  end
  return U
end

function *(O1::NfAbsOrd, O2::NfAbsOrd)
  k = nf(O1)
  @assert k === nf(O2)
  b1 = basis(O1, k)
  n = length(b1)
  b2 = basis(O2, k)
  p = [x*y for (x,y) in Base.Iterators.ProductIterator((b1, b2))]
  d = reduce(lcm, [denominator(x) for x = p])
  M = zero_matrix(FlintZZ, n*n, n)
  z = fmpz()
  for i = 1:n*n
    a = p[i]*d
    Hecke.elem_to_mat_row!(M, i, z, a)
  end
  h = hnf(M)
  b = nf_elem[]
  for i=1:n
    push!(b, Hecke.elem_from_mat_row(k, h, i, d))
  end
  return Hecke.Order(k, b)
end

mutable struct GeIdeal
  a::NfAbsOrdIdl
  function GeIdeal(a::NfAbsOrdIdl)
    o =order(a)
    if o.ismaximal == 1
      return new(a)
    end
    #keep track of known maximality and use this to speed things up
    o = ring_of_multipliers(a)
    if o === order(a)
      return new(a)
    else
      return new(extend(a, o))
    end
  end
  function GeIdeal(a::nf_elem)
    o = Hecke.any_order(parent(a))
    d = denominator(a, o)
    return new(o(a*d)*o), new(o(d)*o)
  end
end

import Hecke.gcd, Hecke.isone, Hecke.*, Hecke.gcd_into!, Hecke.copy, Hecke.divexact, 
       Hecke.isunit, Hecke.coprime_base, Hecke.valuation

function make_compatible!(a::GeIdeal, b::GeIdeal)
  o1 = order(a.a)
  o2 = order(b.a)
  if o1 === o2
    return
  end
  o = o1*o2
  a.a = extend(a.a, o)
  b.a = extend(b.a, o)
end

#isunit, divexact, gcd, isone, *, gcd_into!, copy
isone(a::GeIdeal) = isone(a.a)
isunit(a::GeIdeal) = isone(a)
copy(a::GeIdeal) = GeIdeal(a.a)

function gcd(a::GeIdeal, b::GeIdeal)
  make_compatible!(a, b)
  return GeIdeal(a.a + b.a)
end

function gcd_into!(a::GeIdeal, b::GeIdeal, c::GeIdeal)
  make_compatible!(b, c)
  a.a = b.a + c.a
  return a
end

function *(a::GeIdeal, b::GeIdeal)
  make_compatible!(a, b)
  return GeIdeal(a.a * b.a)
end

function divexact(a::GeIdeal, b::GeIdeal)
  make_compatible!(a, b)
  return GeIdeal(divexact(a.a, b.a))
end

Hecke.norm(a::GeIdeal) = norm(a.a)

function coprime_base(A::Array{nf_elem, 1})
  c = Array{GeIdeal, 1}()
  for a = A
    n,d = GeIdeal(a)
    isone(n) || push!(c, n)
    isone(d) || push!(c, d)
  end
  return coprime_base(c)
end

function coprime_base(A::Array{nf_elem, 1}, O::NfAbsOrd)
  c = Array{NfAbsOrdIdl{AnticNumberField, nf_elem}, 1}()
  for a = A
    n,d = integral_split(a*O)
    isone(n) || push!(c, n)
    isone(d) || push!(c, d)
  end
  return coprime_base(c)
end


function valuation(a::nf_elem, p::GeIdeal)
  return valuation(a, p.a)
end

end
