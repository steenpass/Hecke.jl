@testset "NfRel" begin
  @testset "issubfield" begin
    global Qx, x = QQ["x"]
    f = x^2 + 12x - 92
    global K, a = NumberField(f, "a")
    Ky, y = K["y"]

    L, b = NumberField(y^2 + y + 1, "b")
    M, c = NumberField(y^6 + y^3 + 1, "c")

    d, LtoM = Hecke.issubfield(L, M)

    @test d == true
    @test parent(LtoM(b)) == M
  end

  @testset "isisomorphic" begin
    global Qx, x = QQ["x"]
    f = x^2 + 12x - 92
    global K, a = NumberField(f, "a")
    Ky, y = K["y"]

    g = y^5 - 5
    L, b = NumberField(g, "b")
    bb = 5*b - 2
    h = minpoly(bb)
    L2, b2 = NumberField(h, "b2")

    c, LtoL2 = Hecke.isisomorphic(L, L2)
    @test c == true
    @test parent(LtoL2(b)) == L2

    i = g - 1
    L3, b3 = NumberField(i, "b3")
    d, LtoL3 = Hecke.isisomorphic(L, L3)
    @test d == false
  end
end
