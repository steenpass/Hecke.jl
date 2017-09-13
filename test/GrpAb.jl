@testset "Abelian groups (GrpAb)" begin
  include("GrpAb/GrpAbFinGen.jl")
  include("GrpAb/SubgroupEnum.jl")
  include("GrpAb/Lattice.jl")
end

#Base.Test.print_test_results(r, 3) 
