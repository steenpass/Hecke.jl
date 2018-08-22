@testset "Elements" begin
  @testset "Constructors" begin
    M = FlintZZ[1 2 3; 4 5 6]
    G = @iinfered AbelianGroup(M)
    N = FlintZZ[1 2 3]
    a = @iinfered GrpAbFinGenElem(G, N)
    @test parent(a) == G
    @test a.coeff == N

    G = @iinfered DiagonalGroup([3, 0])
    N = FlintZZ[1 1]
    a = @iinfered GrpAbFinGenElem(G, N)
    @test @iinfered parent(a) == G
    @test a.coeff == N
  end

  @testset "Generators" begin
    M = FlintZZ[1 2 3; 4 5 6]
    G = AbelianGroup(M)
    ge = @iinfered gens(G)
    @test length(ge) == 3
    @test ge[1] == G[1]
    @test ge[2] == G[2]
    @test ge[3] == G[3]
  end

  @testset "Parent" begin
    G = @iinfered DiagonalGroup([3, 0])
    N = FlintZZ[1 1]
    a = @iinfered GrpAbFinGenElem(G, N)
    @test @iinfered parent(a) == G
  end

  @testset "String I/O" begin
    G = DiagonalGroup([3, 0])
    N = FlintZZ[1 1]
    a = GrpAbFinGenElem(G, N)
    @test isa(string(a), String)
  end

  @testset "Hashing" begin
    G = DiagonalGroup([3, 0])
    N = FlintZZ[1 1]
    a = GrpAbFinGenElem(G, N)
    @test isa(hash(a), UInt)
  end

  @testset "Indexing" begin
    G = DiagonalGroup([3, 0])
    N = FlintZZ[1 2]
    a = GrpAbFinGenElem(G, N)
    @test @iinfered a[1] == 1
    @test @iinfered a[2] == 2
  end

  @testset "Comparison" begin
    G = DiagonalGroup([3, 0])
    N = FlintZZ[1 2]
    a = GrpAbFinGenElem(G, N)
    b = GrpAbFinGenElem(G, deepcopy(N))
    @test @iinfered a == b

    H = DiagonalGroup([3, 0])
    c = GrpAbFinGenElem(H, N)
    @test_throws ErrorException a == c
  end

  @testset "Arithmetic" begin
    G = DiagonalGroup([3, 3, 3])
    a = G[1]
    b = G[2]
    c = G([1, 1, 0])
    @test a + b == c
    @test -a == G([2, 0, 0])

    aa = @iinfered(2 * a)
    @test aa == G([2, 0, 0])
    
    aa = @iinfered(a * 2)
    @test aa == G([2, 0, 0])

    
    aa = @iinfered(fmpz(2) * a)
    @test aa == G([2, 0, 0])
  end

  @testset "Neutral element" begin
    G = DiagonalGroup([3, 3, 3])
    a = G[1]
   
    aa = @iinfered(a * fmpz(2))
    @test aa == G([2, 0, 0])

    @test !iszero(a)

    c = G([0, 0, 0])
    @test iszero(c)
  end

  @testset "Parent object overloading" begin
    G = DiagonalGroup([3, 3, 3])
    a = @iinfered G(fmpz[1, 1, 1])
    @test parent(a) == G
    @test a.coeff == FlintZZ[1 1 1]

    a = @iinfered G([1, 1, 1])
    @test parent(a) == G
    @test a.coeff == FlintZZ[1 1 1]

    M = FlintZZ[1 1 1]
    a = @iinfered G(M)
    M[1, 1] = 3
    @test parent(a) == G
    @test a.coeff == FlintZZ[1 1 1]

    a = @iinfered G[1]
    @test parent(a) == G
    @test a.coeff == FlintZZ[1 0 0]
  end

  @testset "Order" begin
    G = DiagonalGroup([3, 3, 0])
    a = G[1]
    @test @iinfered order(a) == 3

    G = DiagonalGroup([3, 5, 0])
    a = G[1]
    @test @iinfered order(a) == 3

    a = G[3]
    @test_throws ErrorException order(a)
  end

  @testset "Random elements" begin
    G = DiagonalGroup([3, 5])
    a = @iinfered rand(G)
    @test parent(a) == G

    G = DiagonalGroup([3, 15])
    a = @iinfered rand(G)
    @test parent(a) == G

    G = DiagonalGroup([3, 0])
    @test_throws ErrorException rand(G)

    a = @iinfered rand(G, 10)
    @test parent(a) == G
    @test -10 <= a[2] <= 10

    a = @iinfered rand(G, fmpz(10))
    @test parent(a) == G
    @test -10 <= a[2] <= 10

    G = DiagonalGroup([3, 0, 5, 0])
    @test_throws ErrorException rand(G)

    a = @iinfered rand(G, 10)
    @test parent(a) == G
    @test -10 <= a[2] <= 10
    @test -10 <= a[4] <= 10

    a = @iinfered rand(G, fmpz(10))
    @test parent(a) == G
    @test -10 <= a[2] <= 10
    @test -10 <= a[4] <= 10
  end

  @testset "Iterator" begin
    G = DiagonalGroup([2, 0])
    @test_throws ErrorException begin for a in G end end
    G = DiagonalGroup([fmpz(2)^100])
    @test_throws ErrorException begin for a in G end end

    G = DiagonalGroup([3, 5, 10])
    @test length(G) == 3*5*10
    @test length(collect(G)) == 3*5*10

    G = DiagonalGroup([3, 9, 27])
    @test length(G) == 3*9*27
    @test length(collect(G)) == 3*9*27
  end

 
  @testset "Helper" begin
    @testset "Reduce mod Hermite normal form" begin
      a = FlintZZ[21 32 43]
      H = FlintZZ[2 0 0 ; 0 3 0 ; 0 0 5]
      Hecke.reduce_mod_hnf!(a, H)
      @test a == FlintZZ[1 2 3]

      a = FlintZZ[1 3 42]
      H = FlintZZ[1 1 14 ; 0 2 11 ; 0 0 17]
      Hecke.reduce_mod_hnf!(a, H)
      @test a == FlintZZ[0 0 0]

      a = FlintZZ[0 0 1]
      H = FlintZZ[1 32 62 ; 0 45 90 ; 0 0 0]
      Hecke.reduce_mod_hnf!(a, H)
      @test a == FlintZZ[0 0 1]
    end
    
    @testset "Smith normal form with transform" begin
      M = MatrixSpace(FlintZZ,1,1)([0])
      S = MatrixSpace(FlintZZ,1,1)([0])
      T,L,R = snf_with_transform(M, true, true)
      @test S == T
      @test L*M*R == T

      M = MatrixSpace(FlintZZ,1,1)([1])
      S = MatrixSpace(FlintZZ,1,1)([1])
      T,L,R = snf_with_transform(M, true, true)
      @test S == T
      @test L*M*R == T

      M = FlintZZ[834 599 214 915 ; 784 551 13 628 ; 986 5 649 100 ; 504 119 64 310 ]
      S = FlintZZ[1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 36533330310]
      T,L,R = snf_with_transform(M, true, true)
      @test S == T
      @test L*M*R == T
      T,L,R = snf_with_transform(M, false, true)
      T,L,R = snf_with_transform(M, true, false)
      T,L,R = snf_with_transform(M, false, false)
    end
  end
end
