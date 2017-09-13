################################################################################
#
#           GrpAb/Lattice.jl : Lattice of abelian groups
#
# This file is part of Hecke.
#
# Copyright (c) 2015, 2016, 2017: Claus Fieker, Tommy Hofmann
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
#  Copyright (C) 2017 Tommy Hofmann
#
################################################################################

################################################################################
#
#  Graph type with primitive functionality
#
################################################################################

# Note that our graph holds some data at each edge.
# If edges[a][b] = c, then there is an edge from a -> b with some data c.
mutable struct Graph{T, M}
  edges::Dict{T, Dict{T, M}}
  degrees::Dict{T, Int}

  function Graph{T, M}() where {T, M}
    z = new{T, M}()
    z.edges = Dict{T, Dict{T, M}}()
    z.degrees = Dict{T, Int}()
    return z
  end
end

# append vertex to a graph
function Base.append!(G::Graph{T, M}, a::T) where {T, M}
  G.edges[a] = Dict{T, M}()
  G.degrees[a] = 0
end

# append edge to a graph
function Base.append!(G::Graph{T, M}, e::Tuple{T, T}, data::M) where {T, M}
  if !haskey(G.edges, e[1])
    append!(G, e[1])
  end
  if !haskey(G.edges, e[2])
    append!(G, e[2])
  end
  G.edges[e[1]][e[2]] = data
  G.degrees[e[1]] += 1
  G.degrees[e[2]] += 1
end

# delete vertex from graph
function Base.delete!(G::Graph{T, M}, a::T) where {T, M}
  if haskey(G.edges, a)
    for e in keys(G.edges[a])
      G.degrees[e] -= 1
    end
    Base.delete!(G.edges, a)
  end
  for e in keys(G.edges)
    if haskey(G.edges[e], a)
      G.degrees[e] -= 1
      Base.delete!(G.edges[e], a)
    end
  end
  Base.delete!(G.degrees, a)
end

# Use breadth-first algorithm to find shortes path from root to target.
# The return value is the array of visited vertices in reversed order.
function find_shortest(G::Graph{T, M}, root::T, target::T) where {T, M}
  S = T[]
  Swithparent = Tuple{T, T}[]
  Q = T[]
  push!(S, root)
  push!(Q, root)

  found = false

  while length(Q) > 0 && found == false
    current = pop!(Q)
    if current == target
      found = true
      break
    end
    for n in keys(G.edges[current])
      if !(n in S)
        push!(S, n)
        push!(Swithparent, (current, n))
        if n == target
          found = true
          break
        end
        prepend!(Q, n)
      end
    end
  end

  if !found
    return false, T[]
  end

  current = target
  path = [ target ]

  while current != root
    for s in Swithparent
      if s[2] == current
        current = s[1]
        push!(path, current)
        break
      end
    end
  end

  return true, path
end

# Use breadth-first algorithm to check if root1 and root2 have a path to a
# common vertex. Return values are arrays of visited vertices in reversed
# order.
function find_common(G::Graph{T, M}, root1::T, root2::T) where {T, M}
  S1 = T[]
  S1withparent = Tuple{T, T}[]
  Q1 = T[]
  push!(S1, root1)
  push!(Q1, root1)

  S2 = T[]
  S2withparent = Tuple{T, T}[]
  Q2 = T[]
  push!(S2, root2)
  push!(Q2, root2)

  found = false

  current = root1

  while length(Q1) > 0 && found == false
    current = pop!(Q1)

    if current in S2
      found = true
    end

    for n in keys(G.edges[current])
      if !(n in S1)
        push!(S1, n)
        push!(S1withparent, (current, n))
        prepend!(Q1, n)
      end
    end
  end
  
  while length(Q2) > 0  && found == false
    current = pop!(Q2)
    if current in S1
      found = true
      break
    end
    for n in keys(G.edges[current])
      if !(n in S2)
        push!(S2, n)
        push!(S2withparent, (current, n))
        prepend!(Q2, n)
      end
    end
  end

  if found
    target = current
    
    path1 = [ target ]

    while current != root1
      for s in S1withparent
        if s[2] == current
          current = s[1]
          push!(path1, current)
          break
        end
      end
    end
    
    path2 = [ target ]

    current = target

    while current != root2
      for s in S2withparent
        if s[2] == current
          current = s[1]
          push!(path2, current)
          break
        end
      end
    end
    return true, path1, path2
  end

  false, T[], T[]
end


################################################################################
#
#  Lattice of groups
#
################################################################################

mutable struct GrpAbLattice
  nvertices::Int                    # Number of vertices.
  weak_vertices::WeakKeyDict        # Store the vertices as weak references.
  graph::Graph{UInt, fmpz_mat}      # Abstract graph with vertices the
                                    # object_id's.
  block_gc::Dict{GrpAbFinGen, Bool} # Groups go here if they must not be gc'ed.

  function GrpAbLattice()
    z = new()
    z.weak_vertices = WeakKeyDict()
    z.graph = Graph{UInt, fmpz_mat}()
    z.nvertices = 0
    z.block_gc = Dict{GrpAbFinGen, Bool}()
    return z
  end
end

# the global group lattice
const GroupLattice = GrpAbLattice()

# append group to group lattice
function Base.append!(G::GrpAbLattice, a::GrpAbFinGen)
  @show "appending $a with"
  @show object_id(a)
  G.nvertices += 1
  G.weak_vertices[a] = true
  append!(G.graph, object_id(a))
end

# append map to group lattice
function Base.append!(G::GrpAbLattice, f::Map)
  if haskey(G.weak_vertices, domain(f)) && haskey(G.weak_vertices, codomain(f))
    error("super bad")
  end
  _is_consistent(G)
  if !haskey(G.weak_vertices, domain(f))
    append!(G, domain(f))
  end
  if !haskey(G.weak_vertices, codomain(f))
    append!(G, codomain(f))
  end
  append!(G.graph, (object_id(domain(f)), object_id(codomain(f))), deepcopy(f.map))
  if G.graph.degrees[object_id(domain(f))] > 1
    G.block_gc[domain(f)] = true
  end
  if G.graph.degrees[object_id(codomain(f))] > 1
    G.block_gc[codomain(f)] = true
  end
end

# update the group lattice after a group was gc'ed.
function update!(G::GrpAbLattice, A::GrpAbFinGen)
  #@show length(G.graph.edges)
  #@show length(G.weak_vertices)
  # find out which group was gc'ed
  delete!(G.weak_vertices, A)
  deleted = setdiff(keys(G.graph.edges), (object_id(k) for k in keys(G.weak_vertices)))
  #@show deleted
  for k in deleted
    Base.delete!(G.graph, k)
    G.nvertices -= 1
  end

  # This is only run after something has been deleted
  # So we only need to delete stuff from block_gc
  for a in keys(G.block_gc)
    @assert haskey(G.graph.degrees, object_id(a))
    if G.graph.degrees[object_id(a)] < 2
      Base.delete!(G.block_gc, a)
    end
  end
end

# Test if in the given group lattice L, there is a path from G to H.
# If so, the matrix describing this map returned.
# TODO: Reduce the entries of the final matrix
function can_map_into(L::GrpAbLattice, G::GrpAbFinGen, H::GrpAbFinGen)
  b, p = find_shortest(L.graph, object_id(G), object_id(H))
  if b
    l = length(p)
    m = L.graph.edges[p[l]][p[l-1]]
    for i in l-1:-1:2
      m = m*(L.graph.edges[p[i]][p[i-1]])
    end
    return true, m
  else
    return false, fmpz_mat(0, 0)
  end
end

# Test if in the given group lattice L, the groups G and H can be mapped into
# a common group. If so, return this group and the matrices describing the maps.
function can_map_into_overstructure(L::GrpAbLattice, G::GrpAbFinGen,
                                                     H::GrpAbFinGen)
  b, pG, pH = find_common(L.graph, object_id(G), object_id(H))
  if b
    @assert pG[1] == pH[1]
    M = G
    found = false
    for g in keys(L.weak_vertices)
      if object_id(g) == pG[1]
        M = g
        found = true
        break
      end
    end

    @assert found

    lG = length(pG)
    mG = L.graph.edges[pG[lG]][pG[lG - 1]]
    for i in lG-1:-1:2
      mG = mG*(L.graph.edges[pG[i]][pG[i - 1]])
    end

    lH = length(pH)
    mH = L.graph.edges[pH[lH]][pH[lH - 1]]
    for i in lH-1:-1:2
      mH = mH*(L.graph.edges[pH[i]][pH[i-1]])
    end

    return true, M, mG, mH
  else
    return false, G, fmpz_mat(0, 0), fmpz_mat(0, 0)
  end
end

# This registered at the finalizer of each group in the global group lattice
function _delete_and_clean_up!(G::GrpAbFinGen)
  update!(GroupLattice, G)
  #_is_consistent(GroupLattice)
  nothing
end

function _is_consistent(L::GrpAbLattice)
  @assert length(L.weak_vertices) == length(L.graph.edges)
  G = L.graph
  @show L
  for k in keys(L.weak_vertices)
    @assert haskey(G.edges, object_id(k))
    @assert haskey(G.degrees, object_id(k))
  end
  for k in keys(L.block_gc)
    @assert haskey(G.edges, object_id(k))
    @assert haskey(G.degrees, object_id(k))
  end
  #println("consistent")
end
