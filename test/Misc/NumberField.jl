@testset "Misc/NumberField" begin
  @testset "issubfield" begin
    global Qx, x = QQ["x"]
    K, a = NumberField(x^2 + 1, "a")
    L, b = NumberField(x^4 + 1, "b")

    c, KtoL = Hecke.issubfield(K, L)
    @test c == true
    @test parent(KtoL(a)) == L

    OK = MaximalOrder(K)
    OL = MaximalOrder(L)
    c, KtoL = Hecke.issubfield(K, L)
    @test c == true
    @test parent(KtoL(a)) == L
  end

  @testset "isisomorphic" begin
    global Qx, x = QQ["x"]
    f = x^5 + 12x - 92
    K, a = NumberField(f, "a")

    g = x^5 - 172x^4 + 7024x^3 + 8656448x^2 + 55735552x + 45796197888
    global K2, a2 = NumberField(g, "a2")

    c, KtoK2 = Hecke.isisomorphic(K, K2)
    @test c == true
    @test parent(KtoK2(a)) == K2

    OK = MaximalOrder(K)
    OK2 = MaximalOrder(K2)
    c, KtoK2 = Hecke.isisomorphic(K, K2)
    @test c == true
    @test parent(KtoK2(a)) == K2

    h = f - 1
    global K3, a3 = NumberField(h, "a3")
    d, KtoK3 = Hecke.isisomorphic(K, K3)
    @test d == false
  end
end
