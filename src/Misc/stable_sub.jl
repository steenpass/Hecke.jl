
add_verbose_scope(:StabSub)

###############################################################################
#
#  Tools for ZpnGModules
#
###############################################################################

#
#  Lifts a matrix from F_p to Z/p^nZ
#

function lift(M::fq_nmod_mat, R::Nemo.NmodRing)

  x=factor(fmpz(R.n))
  @assert length(x.fac)==1
  @assert order(parent(M[1,1]))==first(keys(x.fac))
  N=zero_matrix(R,rows(M),cols(M))
  for i=1:rows(M)
    for j=1:cols(M)
      N[i,j]=FlintZZ(coeff(M[i,j],0))
    end
  end
  return N
  
end


#
#  Action of a matrix on an element of the group
#

function *(x::GrpAbFinGenElem, M::nmod_mat)

  G=parent(x)
  @assert ngens(G)==rows(M)
  R=parent(M[1,1]) 
  coeff=zero_matrix(R,1,cols(M))
  for i=1:cols(M)
    coeff[1,i]=x.coeff[1,i]
  end
  y=coeff*M
  return G(lift(y))
  
end


#
#  Smith Normal Form for a ZpnGModule
#


function Nemo.snf(M::ZpnGModule)

  A=M.V
  G=M.G
  if issnf(A)
    return M, GrpAbFinGenMap(A,A, identity_matrix(FlintZZ, ngens(A)))
  end
  S,mS=snf(A)
  W=[mS(s) for s in gens(S)]
  R=M.R
  H=Array{nmod_mat,1}(undef, length(G))
  for i=1:length(G)
    N=zero_matrix(R, ngens(S),ngens(S))
    for j=1:length(W)
      y=mS\(W[j]*G[i])
      for k=1:ngens(S)
        N[j,k]=y.coeff[1,k]
      end
    end
    H[i]=N
  end
  return ZpnGModule(S,H), mS
    
end

function is_stable(act::Array{T, 1}, mS::GrpAbFinGenMap) where T <: Map{GrpAbFinGen, GrpAbFinGen} 

  S=mS.header.domain
  for s in gens(S)
    x=mS(s)
    for g in act
      if !haspreimage(mS,g(x))[1]
        return false
      end
    end
  end
  return true

end

function issubmodule(M::ZpnGModule, S::nmod_mat)
  
  @assert issnf(M.V)
  s, ms=subm_to_subg(M,S)
  for x in gens(s)
    el=ms(x)
    for g in M.G
      if !haspreimage(ms,el*g)[1]
        return false
      end
    end
  end
  return true
  
end

#
#  Given a group $G$ and a group of automorphisms of G, this function returns the corresponding ZpnGModule
#

function action(V::GrpAbFinGen, act::Array{T,1}) where T<: Map{GrpAbFinGen, GrpAbFinGen} 

  expon=Int(exponent(V))
  @assert length(factor(order(V)).fac)==1
  RR=ResidueRing(FlintZZ, expon, cached=false)
  act_mat=Array{nmod_mat,1}(undef, length(act))
  for z=1:length(act)
    A=zero_matrix(RR,ngens(V), ngens(V))
    for i=1:ngens(V)
      y=act[z](V[i])
      for j=1:ngens(V)
        A[i,j]=y.coeff[1,j]
      end
    end
    act_mat[z]=A
  end
  return ZpnGModule(V,act_mat)

end


#################################################################################
#
#  Duality
#
#################################################################################

#=
  To get the transpose map of a homomorphism, and so the action on the dual module,
  you need to transpose the matrix and multiply or the elements for a power of p, assuming that the group is in snf.
  Precisely, let p^e be the exponent of the group and let p^t_i be the powers of p appearing in the snf. Then 
  ( p^(e-t_1)     )                  ( p^(e-t_1)     )
  (          ...  )  x ^t A  = A  x  (          ...  )
  (              1)                  (              1)
=#


function dual_module(M::ZpnGModule)

  @assert issnf(M.V)
  G1=deepcopy(M.G)
  for i=1:length(G1)
    for j=1:ngens(M.V)-1
      for k=j+1:ngens(M.V)
        x=divexact(M.V.snf[k], M.V.snf[j])
        G1[i][j,k]=divexact(G1[i][j,k].data, x)
        G1[i][k,j]*=x
      end
    end
    transpose!(G1[i])
  end 
  return ZpnGModule(M.V, G1)
  
end

function _dualize(M::nmod_mat, V::GrpAbFinGen, v::Array{fmpz,1})    
  #
  #  First, compute the kernel of the corresponding homomorphisms
  # 
  K=DiagonalGroup([V.snf[end] for j=1:rows(M)])
  A=lift(transpose(M))
  for j=1:rows(A)
    for k=1:cols(A)
      A[j,k]*=v[j]
    end
  end 
  mH=Hecke.GrpAbFinGenMap(V,K,A)
  newel=kernel_as_submodule(mH)
  return MatrixSpace(parent(M[1,1]),rows(newel), cols(newel), false)(newel)
  
end

function _dualize_1(M::nmod_mat, snf_struct::Array{fmpz,1})

  A=nullspace(transpose(M))
  B=vcat(transpose(A),zero_matrix(M[1,1].parent, cols(A),cols(A)))
  for j=1:cols(A)
    B[rows(A)+j,j]=snf_struct[j]
  end
  S=nullspace(B)
  C=vcat(transpose(A),zero_matrix(M[1,1].parent, cols(A),cols(A)))
  return S*C
 
end

function kernel_as_submodule(h::GrpAbFinGenMap)
  G = domain(h)
  H = codomain(h)
  hn, t = hnf_with_transform(vcat(h.map, rels(H))) 
  for i=1:rows(hn)
    if iszero_row(hn, i)
      return view(t, i:rows(t), 1:ngens(G))
    end
  end
  error("JH")
end



########################################################################################################
#
#  Quotients and subgroups of ZpnGModules
#
########################################################################################################


function quo(M::ZpnGModule, S:: nmod_mat)

  Q,mQ=quo(M.V,lift(S))
  return ZpnGModule(Q,M.G), mQ
  
end

function sub(M::ZpnGModule, S::nmod_mat)

  sg,msg=sub(M.V,lift(S))
  G=Array{nmod_mat,1}(undef, length(M.G))
  for k=1:length(M.G)
    A=zero_matrix(M.R, ngens(sg), ngens(sg))
    for i=1:ngens(sg)
      x=msg(sg[i])*M.G[k]
      x=haspreimage(msg, x)[2].coeff
      for j=1:ngens(sg)
        A[i,j]=x[1,j]
      end
    end
    G[k]=A
  end
  return ZpnGModule(sg,G), msg

end


function sub(M::ZpnGModule, n::Int)
  
  sg,msg=sub(M.V,n)
  G=Array{nmod_mat,1}(undef, length(M.G))
  for k=1:length(M.G)
    A=zero_matrix(M.R, ngens(sg), ngens(sg))
    for i=1:ngens(sg)
      x=msg(sg[i])*M.G[k]
      x=haspreimage(msg, x)[2].coeff
      for j=1:ngens(sg)
        A[i,j]=x[1,j]
      end
    end
    G[k]=A
  end
  return ZpnGModule(sg,G), msg
  
end


function _exponent_p_sub(M::ZpnGModule)

  @assert issnf(M.V)
  G=M.G
  V=M.V
  p=M.p
  F, a = Nemo.FiniteField(p, 1, "a", cached=false)
  v=fmpz[divexact(V.snf[i],p) for i=1:ngens(V)]
  G1=Array{fq_nmod_mat,1}(undef, length(G))
  MS=MatrixSpace(F,ngens(V),ngens(V), false)
  for s=1:length(G1)
    G1[s]=MS(0)
    for i=1:ngens(V)
      for j=1:ngens(V)
        if G[s][i,j]!=0 && v[i] <= v[j]
          G1[s][i,j]=divexact((G[s][i,j].data)*v[i],v[j])
        end
      end
    end
  end
  return FqGModule(G1)
  
end

function subm_to_subg(M::ZpnGModule, S::nmod_mat; op=sub)
  
  return op(M.V,lift(S))
  
end


##########################################################################
#
#  Minimal Submodules function
#
##########################################################################

function minimal_submodules(M::ZpnGModule, ord::Int=-1)
  
  
  R=M.R
  S,mS=snf(M)
  N=_exponent_p_sub(S)
  if ord==-1
    list_sub=minimal_submodules(N)
  else
    list_sub=minimal_submodules(N, ngens(S.V)-ord)
  end
  list=Array{nmod_mat,1}(undef, length(list_sub))
  v=[M.p^(valuation(S.V.snf[i], M.p)-1) for i=1:ngens(S.V)]
  W=MatrixSpace(R,1, ngens(M.V), false)
  for z=1:length(list)
    list[z]= vcat([W((mS(S.V([FlintZZ(coeff(list_sub[z][k,i],0))*v[i] for i=1:ngens(S.V)]))).coeff) for k=1:rows(list_sub[z])])
  end
  return list

end

###########################################################################################
#
#  Maximal Submodules function
#
###########################################################################################


function maximal_submodules(M::ZpnGModule, ind::Int=-1)
  
  R=M.R
  S,mS=snf(M)
  N=dual_module(S)
  if ind==-1
    minlist=minimal_submodules(N)
  else 
    minlist=minimal_submodules(N,ind)
  end
  list=Array{nmod_mat,1}(undef, length(minlist))
  W=MatrixSpace(R,1,ngens(S.V), false)
  v=[divexact(fmpz(R.n),S.V.snf[j]) for j=1:ngens(S.V) ]
  for x in minlist
    K=DiagonalGroup([fmpz(R.n) for j=1:rows(x)])
    A=lift(transpose(x))
    for j=1:rows(A)
      for k=1:cols(A)
        A[j,k]*=v[j]
      end
    end 
    mH=Hecke.GrpAbFinGenMap(S.V,K,A)
    sg,msg=kernel(mH)
    push!(list, vcat([ (mS(msg(y))).coeff for y in gens(sg)]))
  end
  return list

end


##########################################################################################
#
#  Composition factors and series
#
##########################################################################################

#
#  Given a list of square matrices G, it returns a list of matrices given by the minors 
#  (n-s) x (n-s) of the matrices G[i] mod p 
#

function _change_ring(G::Array{nmod_mat,1}, F::Nemo.FqNmodFiniteField, s::Int)
  
  G1=Array{fq_nmod_mat,1}(undef, length(G))
  n=rows(G[1])
  for i=1:length(G)
    M=zero_matrix(F,n-s+1,n-s+1)
    for j=s:n
      for k=s:n
        M[j-s+1,k-s+1]=(G[i][j,k]).data
      end
    end
    G1[i]=M
  end
  return G1

end

#
#  Cut the module in submodules with exponent p, returning the quotients p^i M /p^(i+1) M
#

function _mult_by_p(M::ZpnGModule)

  G=M.G
  V=M.V
  p=M.p
  @assert issnf(V)
  #
  #  First, the quotient M/pM
  #
  F,a=Nemo.FiniteField(p,1,"a", cached=false)
  n=ngens(V)
  Gq=_change_ring(G,F,1)
  spaces=FqGModule[FqGModule(Gq)]
  #
  #  Now, the others
  #
  s=valuation(fmpz(M.R.n),p)
  j=1
  for i=2:s
    while !divisible(V.snf[j],p^i)
      j+=1
    end
    GNew=_change_ring(G,F,j)
    push!(spaces, FqGModule(GNew))
  end
  return spaces  

end

function composition_factors(M::ZpnGModule)

  N,_=snf(M)
  list_sub=_mult_by_p(N)
  list=[]
  for x in list_sub
    append!(list, composition_factors(x))
  end
  for i=1:length(list)
    j=i+1
    while j<=length(list)
      if Hecke.isisomorphic(list[i][1], list[j][1])
        list[i][2]+=list[j][2]
        deleteat!(list,j)
      else 
        j+=1
      end    
    end
  end
  return list

end

#######################################################################################
#
#  Submodules function
#
#######################################################################################

Markdown.doc"""
    submodules(M::ZpnGModule; typequo, typesub, order)
> Given a ZpnGModule M, the function returns all the submodules of M. 
> 
"""

function submodules(M::ZpnGModule; typequo=Int[-1], typesub=Int[-1], ord=-1)

  if typequo!=[-1]
    return submodules_with_quo_struct(M,typequo)
  elseif typesub!=[-1]
    return submodules_with_struct(M,typesub)
  elseif ord!=-1
    return submodules_order(M,ord)
  else 
    return submodules_all(M)
  end
  
end


function submodules_all(M::ZpnGModule)
  
  R=M.R
  S,mS=snf(M)  
  minlist=minimal_submodules(S)
  list=nmod_mat[MatrixSpace(R,length(S.V.snf),length(S.V.snf), false)(1), MatrixSpace(R,1,length(S.V.snf), false)(0)]
  
  #
  #  Find minimal submodules, factor out and recursion on the quotients
  #
  
  for x in minlist
    N,_=quo(S,x)
    newlist=submodules_all(N)
    for y in newlist
      push!(list,vcat(y,x))
    end
  end
  
  append!(list,minlist)
  #
  #  Eliminate redundance via Howell form
  #
  w=fmpz[divexact(fmpz(R.n), S.V.snf[j]) for j=1:ngens(S.V)]
  list=_no_redundancy(list,w)
  
  #
  #  Writing the submodules in terms of the given generators and returning an iterator
  #
  W=MatrixSpace(R,rows(mS.map), cols(mS.map), false)
  MatSnf=W(mS.map)
  return (x*MatSnf for x in list)
  
end



function submodules_with_struct_cyclic(M::ZpnGModule, ord::Int)

  R=M.R
  #
  #  We search for an element in p^(ord-1)*G
  #
  s,ms=sub(M, M.p^(ord-1))
  S,mS=snf(s)
  N=_exponent_p_sub(S)
  submod=fq_nmod_mat[]
  if N.dim==1
    push!(submod, identity_matrix(N.K, 1))
  else
    @vtime :StabSub 1 submod=minimal_submodules(N,1,composition_factors(N))
  end
  list1=Array{nmod_mat,1}(undef, length(submod))
  v=fmpz[(M.p)^(valuation(S.V.snf[i], M.p)-1) for i=1:ngens(S.V)]
  for i=1:length(submod)
    @views list1[i]=lift(submod[i],R)
    @assert rows(list1[i])==1
    for k=1:cols(list1[i])
      @inbounds list1[i][1,k]*=v[k]
    end
  end  
  W=MatrixSpace(R,rows(mS.map), cols(ms.map), false)
  MatSnf=W(mS.map*ms.map)  
  for j=1:length(list1)
    @inbounds @views list1[j]=list1[j]*MatSnf
  end
  if ord==1
    return list1
  end
  list=nmod_mat[]
  for x in list1  
    L,_=quo(M,x)
    newlist=submodules_with_struct_cyclic(L,ord-1)
    i=1
    while i<=length(newlist)
      @views if iszero(M.V(lift(M.p^(ord-1)*newlist[i])))
        deleteat!(newlist,i)
      else 
        i+=1
      end
    end
    append!(list, newlist)
  end
  return list

end

function submodule_with_struct_exp_p(M::ZpnGModule, l::Int)
    
    R=M.R
    S,mS=snf(M)
    N=Hecke._exponent_p_sub(S)
    lf=composition_factors(N)
    list=nmod_mat[]
    v=fmpz[divexact(S.V.snf[i], M.p) for i=1:ngens(S.V)]
    submod=submodules(N,ngens(S.V)-l,comp_factors=lf)
    list1=Array{nmod_mat,1}(undef, length(submod))
    for i=1:length(submod) 
      list1[i]=lift(submod[i],R) 
      for s=1:rows(list1[i])
        for t=1:cols(list1[i])
          @inbounds list1[i][s,t]*=v[t]
        end
      end
    end
    W=MatrixSpace(R,rows(mS.map), cols(mS.map), false)
    MatSnf=W(mS.map)
    return ( x *MatSnf for x in list1)
    
end



function submodules_with_struct(M::ZpnGModule, typesub::Array{Int,1})
  
  # If the group is cyclic, it is easier 
  if length(typesub)==1
    return submodules_with_struct_cyclic(M,typesub[1])
  end
  sort!(typesub)
  # If the subgroups we are searching for have exponent p, it is easier
  if typesub[end]==1
    return submodule_with_struct_exp_p(M,length(typesub))
  end
  R=M.R
  
  a=1
  while a<length(typesub) && typesub[end-a]==typesub[end]
    a+=1
  end
  
  S1,mS1=snf(M)
  s,ms=sub(S1, M.p^(typesub[end]-1))
  S,mS=snf(s)
  N=_exponent_p_sub(S)
  @vtime :StabSub 1 submod=submodules(N,(N.dim)-a)
  #
  #  Write the modules as elements of S
  #
  list1=Array{nmod_mat,1}(undef, length(submod))
  v=[divexact(S.V.snf[i], M.p) for i=1:ngens(S.V)]
  for i=1:length(submod)
    @inbounds list1[i]=lift(submod[i],R)
    for j=1:rows(list1[i])
      for k=1:cols(list1[i])
        @inbounds list1[i][j,k]*=v[k]
      end
    end 
  end 
  
  
  auxmat=mS.imap*ms.map
  W=MatrixSpace(R,rows(auxmat), cols(auxmat), false)
  auxmat=W(auxmat)
  for j=1:length(list1)
    @inbounds list1[j]=list1[j]*auxmat
  end
  #
  #  I create the group to check if the candidates are isomorphic to it
  #

  l=[M.p^(typesub[j]) for j=1:length(typesub)]
  G=DiagonalGroup(l)
  G.issnf=true
  G.snf=l
  
  #
  #  I create the type of the group I am searching for in the quotient
  #
  new_subtype=deepcopy(typesub)
  for j=length(typesub)-a+1:length(typesub)
    @inbounds new_subtype[j]-=1
  end
  
  #
  #  Recursion on the quotient
  #
  
  list=nmod_mat[]
  for x in list1  
    L, _=quo(S1,x)
    newlist=collect(submodules_with_struct(L,new_subtype))
    for j=1:length(newlist)
      @inbounds newlist[j]=vcat(newlist[j],x)
    end
    i=1
    while i<=length(newlist)
      t,mt=subm_to_subg(S1,newlist[i])
      if !isisomorphic(t,G)
        deleteat!(newlist,i)
      else 
        i+=1
      end
    end
    if !isempty(newlist)
      append!(list, newlist)
    end
  end
  #
  #  Check for redundancy
  #
  w=fmpz[divexact(fmpz(R.n), S1.V.snf[j]) for j=1:ngens(S1.V)]
  list=_no_redundancy(list,w)

  #
  #  Write the submodules in terms of the set of given generators
  #  and return an iterator over the list
  #

  W=MatrixSpace(R,ngens(S1.V), ngens(M.V), false)
  MatSnf=W(mS1.map)
  
  return (el*MatSnf for el in list)
  
end


function _no_redundancy(list::Array{nmod_mat,1},w::Array{fmpz,1})

  
  #
  #  Howell form of every candidate, embedding them in a free module
  #
  for i=1:length(list)
    if cols(list[i])>=rows(list[i])
      @inbounds list[i]=vcat(list[i],zero_matrix(parent(list[i][1,1]),cols(list[i])-rows(list[i]), cols(list[i])))
    end
    for j=1:cols(list[i])
      for k=1:rows(list[i])
        @inbounds list[i][k,j]*=w[j]
      end
    end
    howell_form!(list[i])
    list[i]=sub(list[i],1:cols(list[i]),1:cols(list[i]))
  end
  #
  #  Now, check if they are equal
  #
  i=1
  while i<length(list)
    k=i+1
    while k<=length(list)
      if list[i]==list[k]
        deleteat!(list,k)
      else 
        k+=1
      end
    end
    i+=1
  end
  
  #
  #  Write them again not embedded
  #
  for i=1:length(list)
    for j=1:cols(list[i])
      for k=1:j
        list[i][k,j]=divexact(list[i][k,j].data,w[j])
      end
    end
  end
  return list

end


function submodules_order(M::ZpnGModule, ord::Int)
  
  
  R=M.R
  S,mS=snf(M)
  @assert exponent(S.V)==fmpz(R.n)
  N=Hecke._exponent_p_sub(S)
  lf=composition_factors(N)
  list=nmod_mat[]
  v=fmpz[divexact(S.V.snf[i], M.p) for i=1:ngens(S.V)]
  for i=1:ord-1
    minlist=minimal_submodules(N,i,lf)
    for x in minlist  
      A=lift(x,R) 
      for s=1:rows(A)
        for t=1:cols(A)
          A[s,t]*=v[t]
        end
      end
      L, _=quo(S,A)
      newlist=Hecke.submodules_order(L,ord-i)
      for z=1:length(newlist)
        push!(list,vcat(newlist[z],A))
      end
    end
  end
  
  #
  #  Check for redundancy
  #
  w=fmpz[divexact(fmpz(R.n), S.V.snf[j]) for j=1:ngens(S.V)]
  list=_no_redundancy(list,w)
  
  #
  #  Write the submodules in terms of the set of given generators
  #
  
  W=MatrixSpace(R,rows(mS.map), cols(mS.map), false)
  MatSnf=W(mS.map)
  for j=1:length(list)
    list[j]=list[j]*MatSnf #vcat([W(( mS( S.V([list[j][k,i].data for i=1:ngens(S.V)]))).coeff)  for k=1:rows(list[j])])
  end
  
  #
  #  Add the minimal submodules
  # 
  
  minlist=minimal_submodules(N,ord, lf)
  for x in minlist
    push!(list, vcat([W((mS( S.V([FlintZZ(coeff(x[k,i],0))*((M.p)^(v[i]-1)) for i=1:ngens(S.V)]))).coeff) for k=1:rows(x) ]))
  end
  return (x for x in list)
  
end



function submodules_with_quo_struct(M::ZpnGModule, typequo::Array{Int,1})
  
  R=M.R 
  S,mS=snf(M)
  sort!(typequo)
  l=[(M.p)^typequo[i] for i=1:length(typequo)]
  wish=DiagonalGroup(l)
  wish.issnf=true
  wish.snf=l
  if isisomorphic(wish,S.V)
    return nmod_mat[zero_matrix(R, 1, ngens(M.V))]
  end
  if length(l)>length(S.V.snf)
    return nmod_mat[]
  end
  for i=1:length(typequo)
    if !divisible(S.V.snf[ngens(S.V)+1-i],fmpz((M.p)^typequo[length(typequo)+1-i]))
      return nmod_mat[]
    end
  end

  #
  #  Dual Module and candidate submodules
  #
  M_dual=dual_module(S)
  candidates=submodules_with_struct(M_dual, typequo)

  #
  #  Dualize the modules
  #
  v=[divexact(S.V.snf[end],S.V.snf[j]) for j=1:ngens(S.V) ]
  list=(_dualize(x, S.V, v) for x in candidates) 
   
  #
  #  Write the submodules in terms of the given generators
  #
  W=MatrixSpace(R, rows(mS.map), cols(mS.map), false)
  MatSnf=W(mS.map)
  return (final_check_and_ans(x, MatSnf, M) for x in list)
  
end

function final_check_and_ans(x::nmod_mat, MatSnf::nmod_mat, M::ZpnGModule)
  
  y=x*MatSnf
  @hassert :RayFacElem 1 issubmodule(M, y)
  return y 

end
