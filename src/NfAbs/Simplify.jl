export simplify

Markdown.doc"""
    simplify(K::AnticNumberField; canonical::Bool = false) -> AnticNumberField, NfToNfMor
 > Tries to find an isomorphic field $L$ given by a "nicer" defining polynomial.
 > By default, "nice" is defined to be of smaller index, testing is done only using
 > a LLL-basis of the maximal order.
 > If \texttt{canonical} is set to {{{true}}}, then a canonical defining
 > polynomial is found, where canonical is using the pari-definition of {{{polredabs}}}
 > in http://beta.lmfdb.org/knowledge/show/nf.polredabs.
 > Both version require a LLL reduced basis for the maximal order.
"""
function simplify(K::AnticNumberField; canonical::Bool = false)
  if canonical
    a, f = polredabs(K)
  else
    ZK = lll(maximal_order(K))
    I = index(ZK)^2
    D = discriminant(ZK)
    B = basis(ZK)
    b = ZK(gen(K))
    f = K.pol
    for i=1:length(B)
      ff = minpoly(B[i])
      if degree(ff) < degree(K)
        continue
      end
      id = div(discriminant(ff), D)
      if id<I
        b = B[i]
        I = id
        f = ff
      end
    end
    a = b.elem_in_nf
  end
  Qx,x=PolynomialRing(FlintQQ)
  L = number_field(Qx(f), cached=false)[1]
  m = NfToNfMor(L, K, a)
  return L, m
end

 #a block is a partition of 1:n
 #given by the subfield of parent(a) defined by a
 #the embeddings used are in R
 #K = parent(a)
 # then K has embeddings into the finite field (parent of R[1])
 # given by the roots (in R) of the minpoly of K
 #integers in 1:n are in the same block iff a(R[i]) == a(R[j])
 #the length of such a block (system) is the degree of Q(a):Q, the length
 # of a block is the degree K:Q(a)
 # a is primitive iff the block system has length n
function _block(a::nf_elem, R::Array{fq_nmod, 1}, ap)
  c = FlintZZ()
  for i=0:a.elem_length
    Nemo.num_coeff!(c, a, i)
    setcoeff!(ap, i, c)
  end
#  ap = Ft(Zx(a*denominator(a)))
  s = [ap(x) for x = R]
  b = []
  a = BitSet()
  i = 0
  n = length(R)
  while i < n
    i += 1
    if i in a
      continue
    end
    z = s[i]
    push!(b, findall(x->s[x] == z, 1:n))
    for j in b[end]
      push!(a, j)
    end
  end
  return b
end

#given 2 block systems b1, b2 for elements a1, a2, this computes the
#system for Q(a1, a2), the compositum of Q(a1) and Q(a2) as subfields of K
function _meet(b1, b2)
  b = []
  for i=b1
    for j = i
      for h = b2
        if j in h
          s = intersect(i, h)
          if ! (s in b)
            push!(b, s)
          end
        end
      end
    end
  end
  return b
end

function polredabs(K::AnticNumberField)
  #intended to implement 
  # http://beta.lmfdb.org/knowledge/show/nf.polredabs
  #as in pari
  #TODO: figure out the separation of T2-norms....
  ZK = lll(maximal_order(K))
  I = index(ZK)^2
  D = discriminant(ZK)
  B = basis(ZK)
  b = gen(K)
  f = K.pol
  
  p = 2^20
  d = 1
  while true
    p = next_prime(p)
    R = ResidueRing(FlintZZ, p, cached=false)
    lp = factor(K.pol, R)
    if any(t->t>1, values(lp.fac))
      continue
    end
    d = Base.reduce(lcm, [degree(x) for x = keys(lp.fac)], init = 1)
    if d < degree(f)^2
      break
    end
  end

  F, w = FiniteField(p, d, "w", cached=false)
  Ft, t = PolynomialRing(F, "t", cached=false)
  ap = Ft()
  R = roots(K.pol, F)
  Zx = FlintZZ["x"][1]
  n = degree(K)

  b = _block(B[1].elem_in_nf, R, ap)
  i = 2
  while length(b) < degree(K)
    bb = _block(B[i].elem_in_nf, R, ap)
    b = _meet(b, bb)
    i += 1
  end
  i -= 1
#  println("need to use at least the first $i basis elements...")
  pr = 100
  old = precision(BigFloat)
  E = 1
  while true
    setprecision(BigFloat, pr)
    try
      E = enum_ctx_from_ideal(ideal(ZK, 1), zero_matrix(FlintZZ, 1, 1), prec = pr, TU = BigFloat, TC = BigFloat)
      if E.C[end] + 0.0001 == E.C[end]  # very very crude...
        pr *= 2
        continue
      end
      break
    catch e
      if isa(e, InexactError)
        pr *= 2
        continue
      end
      rethrow(e)
    end
  end

  l = zeros(FlintZZ, n)
  l[i] = 1

  scale = 1.0
  enum_ctx_start(E, matrix(FlintZZ, 1, n, l), eps = 1.01)

  a = gen(K)
  all_a = [a]
  la = length(a)*BigFloat(E.t_den^2)
  Ec = BigFloat(E.c//E.d)
  eps = BigFloat(E.d)^(1//2)

  found_pe = false
  while !found_pe
    while enum_ctx_next(E)
#      @show E.x
      M = E.x*E.t
      q = elem_from_mat_row(K, M, 1, E.t_den)
      bb = _block(q, R, ap)
      if length(bb) < n
        continue
      end
      found_pe = true
#  @show    llq = length(q)
#  @show sum(E.C[i,i]*(BigFloat(E.x[1,i]) + E.tail[i])^2 for i=1:E.limit)/BigInt(E.t_den^2)
      lq = Ec - (E.l[1] - E.C[1, 1]*(BigFloat(E.x[1,1]) + E.tail[1])^2) #wrong, but where?
#      @show lq/E.t_den^2

      if lq < la + eps
        if lq > la - eps
          push!(all_a, q)
  #        @show "new one"
        else
          a = q
          all_a = [a]
          if lq/la < 0.8
  #          @show "re-init"
            enum_ctx_start(E, E.x, eps = 1.01)  #update upperbound
            Ec = BigFloat(E.c//E.d)
          end
          la = lq
  #        @show Float64(la/E.t_den^2)
        end  
      end
    end
    scale *= 2
    enum_ctx_start(E, matrix(FlintZZ, 1, n, l), eps = scale)
    Ec = BigFloat(E.c//E.d)
  end

  setprecision(BigFloat, old)
  all_f = [(x, minpoly(x)) for x=all_a]
  all_d = [abs(discriminant(x[2])) for x= all_f]
  m = minimum(all_d)

  L1 = all_f[findall(i->all_d[i] == m, 1:length(all_d))]

  function Q1Q2(f::PolyElem)
    q1 = parent(f)()
    q2 = parent(f)()
    g = gen(parent(f))
    for j=0:degree(f)
      if isodd(j)
        q2 += coeff(f, j)*g^div(j, 2)
      else
        q1 += coeff(f, j)*g^div(j, 2)
      end
    end
    return q1, q2
  end
  function minQ(A::Tuple)
    a = A[1]
    f = A[2]
    q1, q2 = Q1Q2(f)
    if lead(q1)>0 && lead(q2) > 0
      return (-A[1], f(-gen(parent(f)))*(-1)^degree(f))
    else
      return (A[1], f)
    end
  end

  L2 = [minQ(x) for x=L1]

  function int_cmp(a, b)
    if a==b
      return 0
    end
    if abs(a) == abs(b)
      if a>b
        return 1
      else
        return -1
      end
    end
    return cmp(abs(a), abs(b))
  end

  function il(F, G)
    f = F[2]
    g = G[2]
    i = degree(f)
    while i>0 && int_cmp(coeff(f, i), coeff(g, i))==0 
      i -= 1
    end
    return int_cmp(coeff(f, i), coeff(g, i))<0
  end

  L3 = sort(L2, lt = il)

  return L3[1]
end


