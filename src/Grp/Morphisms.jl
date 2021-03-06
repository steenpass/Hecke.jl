export image, preimage, has_preimage, isinjective,
issurjective, isbijective, automorphisms, morphisms, find_small_group,
inv, *, id_hom, domain, codomain, divisors, multiples

################################################################################
#
#  Morphism functionality
#
################################################################################

@doc Markdown.doc"""
    image(f::GrpGenToGrpGenMor, g::GrpGenElem) -> h::GrpGenElem
Returns the image of $g$ under $f$.
"""
image(f::GrpGenToGrpGenMor, g::GrpGenElem) = f.img[g.i]

@doc Markdown.doc"""
    preimage(f::GrpGenToGrpGenMor, g::GrpGenElem) -> h::GrpGenElem
Returns one element of the preimage of $g$ under $f$.
"""
function preimage(f::GrpGenToGrpGenMor, g::GrpGenElem)
  h = findfirst(x -> f(x) == g, collect(f.domain))
   if h == nothing
     error("$g has no preimage under $f")
   end
   return f.domain[h]
end

@doc Markdown.doc"""
    has_preimage(f::GrpGenToGrpGenMor, g::GrpGenElem) -> (b::Bool, h::GrpGenElem)
Returns whether $f$ has a preimage. If so, the second return value is an
element $h$ with $f(h) = g$.
"""
function has_preimage(f::GrpGenToGrpGenMor, g::GrpGenElem)
  h = findfirst(x -> f(x) == g, collect(f.domain))
   if h == nothing
     return false
   end
   return (true, f.domain[h])
end

@doc Markdown.doc"""
    *(f::GrpGenToGrpGenMor, g::GrpGenToGrpGenMor) -> h::GrpGenToGrpGenMor
Returns the composition (f * g) = g(f)$.
"""
function *(f::GrpGenToGrpGenMor, g::GrpGenToGrpGenMor)
  return GrpGenToGrpGenMor(f.domain, g.codomain, [g(f(x)) for x in collect(f.domain)])
end

@doc Markdown.doc"""
    inv(f::GrpGenToGrpGenMor) -> h::GrpGenToGrpGenMor
Assumes that $f$ is an isomorphism. Returns the inverse of $f$.
"""
function inv(f::GrpGenToGrpGenMor)
  return GrpGenToGrpGenMor(f.codomain, f.domain, collect(f.domain)[sortperm(getindex.(f.img))])
end

function Base.show(io::IO, f::GrpGenToGrpGenMor)
  println(io, "Morphism from group\n", f.domain, " to\n", f.codomain)
end

domain(f::GrpGenToGrpGenMor) = f.domain

codomain(f::GrpGenToGrpGenMor) = f.codomain

id_hom(G::GrpGen) = GrpGenToGrpGenMor(G, G, collect(G))

image(GtoH::GrpGenToGrpGenMor) = subgroup(GtoH.codomain, unique(GtoH.img))

function kernel(GtoH::GrpGenToGrpGenMor)
  G = GtoH.domain
  H = GtoH.codomain
  return subgroup(G, getindex.(Ref(G), findall(x-> GtoH(x) == id(H), collect(G))))
end

function issurjective(GtoH::GrpGenToGrpGenMor)
  return order(GtoH.codomain) == length(unique(GtoH.img))
end
#finite groups
isinjective(GtoH::GrpGenToGrpGenMor) = issurjective(GtoH)

isbijective(GtoH::GrpGenToGrpGenMor) = issurjective(GtoH)

################################################################################
#
#  Find Morphisms
#
################################################################################

# TODO: Cache the orders of the generators of the small_groups.
#       Do not recompute it here.
function find_small_group(G::GrpGen)
  l = order(G)

  elements_by_orders = Dict{Int, Array{GrpGenElem, 1}}()

  for i in 1:l
    g = G[i]
    o = order(g)
    if haskey(elements_by_orders, o)
      push!(elements_by_orders[o], g)
    else
      elements_by_orders[o] = [g]
    end
  end

  candidates = Int[]

  local ordershape

  for j in 1:length(small_groups_1_63[order(G)])
    ordershape = small_groups_1_63[order(G)][j][4]

    candidate = true
    for (o, no) in ordershape
      if !haskey(elements_by_orders, o)
        candidate = false
        break
      else
        elts = elements_by_orders[o]
        if length(elts) != no
          candidate = false
          break
        end
      end
     end

     if candidate
        push!(candidates, j)
     end
  end

  @assert length(candidates) > 0


  sort!(candidates, rev = true)

  idG = id(G)

  for j in candidates
    H = small_groups_1_63[order(G)][j]

    elbyord = [elements_by_orders[order(o)] for o in H[1]]

    it = Iterators.product(elbyord...)

    words = H[2]

    for poss in it
      is_hom = true
      for w in words
        if eval_word(collect(poss), w) != idG
          is_hom = false
          break
        end
      end

      if is_hom
        if length(closure(collect(poss), *, idG)) == order(G)
          # Found it!
          H = small_group(order(G), j)
          return (order(G), j), H, _spin_up_morphism(gens(H), collect(poss))
        end
      end
    end
  end
  throw(error("Could not identify group"))
end

function eval_word(S, w::Vector{Int})
  g = id(parent(S[1]))
  for i in 1:length(w)
    if w[i] > 0
      g = g * S[w[i]]
    else
      g = g * inv(S[-w[i]])
    end
  end
  return g
end

@doc Markdown.doc"""
    automorphisms(G::GrpGen) -> A::Array{GrpGenToGrpGenMor,1}
Returns all group isomorphisms from $G$ to $G$.
"""
function automorphisms(G::GrpGen)
  Gn, GntoG = find_small_group(G)[2:3]
  auts = _automorphisms(Gn)
  return [inv(GntoG) * aut * GntoG for aut in auts]
end

function _automorphisms(G::GrpGen)
  @assert isfrom_db(G)
  i, j = G.small_group_id
  Gdata = small_groups_1_63[i][j]

  l = order(G)

  elements_by_orders = Dict{Int, Array{GrpGenElem, 1}}()

  # TODO: I think the following is cached somewhere (in the database)
  for i in 1:l
    g = G[i]
    o = order(g)
    if haskey(elements_by_orders, o)
      push!(elements_by_orders[o], g)
    else
      elements_by_orders[o] = [g]
    end
  end


  elbyord = [elements_by_orders[order(o)] for o in Gdata[1]]

  it = Iterators.product(elbyord...)

  words::Vector{Vector{Int}} = Gdata[2]

  idG = id(G)

  auts = _aut_group(it, words, idG, order(G))::Vector{Vector{GrpGenElem}}

  # Any element A of auts determines an isomorphism by mapping gens(G)[i] to A[i]

  Ggens = gens(G)

  # TODO: preallocate
  return [_spin_up_morphism(Ggens, a) for a in auts]
end

function _spin_up_morphism(domain::Vector{GrpGenElem}, codomain::Vector{GrpGenElem})
  @assert length(domain) > 0
  @assert length(domain) == length(codomain)
  G = parent(domain[1])
  H = parent(codomain[1])
  pairs = [(domain[i], codomain[i]) for i in 1:length(domain)]
  cl = closure(pairs, (x, y) -> (x[1]*y[1], x[2]*y[2]), (id(G), id(H)))
  img = Vector{GrpGenElem}(undef, length(G))
  for i in 1:length(G)
    img[cl[i][1][]] = cl[i][2]
  end
  phi = GrpGenToGrpGenMor(G, H, img)

  # TODO: Remove this assertion once this is battle tested
  for g in G
    for h in G
      @assert phi(g * h) == phi(g) * phi(h)
    end
  end
  return phi
end

@noinline function _aut_group(it, words, idG, n)
  auts = Vector{GrpGenElem}[]
  for poss in it
    is_hom = true
    for w in words
      if eval_word(poss, w) != idG
        is_hom = false
        break
      end
    end

    if is_hom
      cposs = collect(poss)
      if length(closure(cposs, *, idG)) == n
        push!(auts, cposs)
      end
    end
  end

  return auts
end

function _morphisms_with_gens(G::GrpGen, H::GrpGen, Gens::Array{GrpGenElem,1}, Rels::Array{Array{Int64,1},1})

  l = order(H)

  elements_by_orders = Dict{Int, Array{GrpGenElem, 1}}()

  for i in 1:l
    h = H[i]
    o = order(h)
    for k in multiples(o, order(G))
      if haskey(elements_by_orders, k)
        push!(elements_by_orders[k], h)
      else
        elements_by_orders[k] = [h]
      end
    end
  end


  elbyord = [elements_by_orders[order(o)] for o in Gens]

  it = Iterators.product(elbyord...)

  words::Vector{Vector{Int}} = Rels

  idH = id(H)

  homs = _hom_group(it, words, idH)::Vector{Vector{GrpGenElem}}

  Ggens = Gens

  return [_spin_up_morphism(Ggens, a) for a in homs]
end

@doc Markdown.doc"""
    morphisms(G::GrpGen, H::GrpGen) -> A::Array{GrpGenToGrpGenMor,1}
Returns all group homomorphisms from $G$ to $H$.
"""
function morphisms(G::GrpGen, H::GrpGen)
  Gn, isom = find_small_group(G)[2:3]
  return [inv(isom) * mor for mor in _morphisms(Gn,H)]
end

function _morphisms(G::GrpGen, H::GrpGen)
  @assert isfrom_db(G)
  i, j = G.small_group_id
  Gdata = small_groups_1_63[i][j]

  l = order(H)

  elements_by_orders = Dict{Int, Array{GrpGenElem, 1}}()

  for i in 1:l
    h = H[i]
    o = order(h)
    for k in multiples(o, order(G))
      if haskey(elements_by_orders, k)
        push!(elements_by_orders[k], h)
      else
        elements_by_orders[k] = [h]
      end
    end
  end


  elbyord = [elements_by_orders[order(o)] for o in Gdata[1]]

  it = Iterators.product(elbyord...)

  words::Vector{Vector{Int}} = Gdata[2]

  idH = id(H)

  homs = _hom_group(it, words, idH)::Vector{Vector{GrpGenElem}}

  # Any element a of homs determines an homomorphism by mapping gens(G)[i] to A[i]

  Ggens = gens(G)

  # TODO: preallocate
  return [_spin_up_morphism(Ggens, a) for a in homs]
end

@noinline function _hom_group(it, words, idH)
  homs = Vector{GrpGenElem}[]
  for poss in it
    is_hom = true
    for w in words
      if eval_word(poss, w) != idH
        is_hom = false
        break
      end
    end

    if is_hom
      cposs = collect(poss)
      push!(homs, cposs)
    end
  end
  return homs
end

divisors(n::Int64) =  findall(x -> mod(n,x) == 0, 1:n)

multiples(n::Int64, b::Int64) =  [i * n for i in 1:Int64(floor(b/n))]
