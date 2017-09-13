@testset "Lattice" begin
  G = DiagonalGroup([2, 3, 4])
  H1, mH1 = sub(G, [G[1]])
  H2, mH2 = sub(G, [G[2]])
  H3, mH3 = sub(G, [G[3]])

  @test H1[1] + G[1] == 2*G[1]
  @test H1[1] + H2[1] == G[1] + G[2]
  @test H1[1] + H2[1] + H3[1] == G[1] + G[2] + G[3]

  H1 = 1
  mH1 = 1
  gc()

  @test H2[1] + G[2] == 2*G[2]

  H4, mH4 = sub(H3, [2*H3[1]])

  @test H4[1] + G[1] == G[1] + 2*G[3]

  Q, mQ = quo(G, [2*G[3]])

  G = 1
  gc()
  gc()
  @test H2[1] + Q[3] == 1*Q[2] + 1*Q[3]
end
