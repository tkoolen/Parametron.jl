module SimpleQPTest

using Compat
using Compat.Test
using SimpleQP

@testset "LinearFunction" begin
    x = [Variable(1), Variable(2)]
    A1 = [1.0 2.0; 3.0 4.0]
    f1 = A1 * x
    f2 = 1.0 * A1 * x
    vals = Dict(x => [2.0, 5.0])
    @test f1(vals) == f2(vals)
    @test f1(vals) == A1 * vals[x]

    f3 = 2 * f2
    @test f3(vals) == 2 * f2(vals)

    f4 = f2 * 2.0
    @test f4(vals) == f3(vals)

    A2 = [4.0 5.0; 6.0 7.0]
    f5 = 1.5 * A1 * x + 2.5 * A2 * x
    @test f5(vals) == 1.5 * A1 * vals[x] + 2.5 * A2 * vals[x]

    f6 = -f5
    @test f6(vals) == -f5(vals)

    f7 = A1 * x - 0.5 * A2 * x
    @test f7(vals) == A1 * vals[x] - 0.5 * A2 * vals[x]
end

end
