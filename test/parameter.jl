module ParameterTest

using Compat
using Compat.Test
using SimpleQP

import SimpleQP: setdirty!
import ..MockModel

@testset "Parameter" begin
    # out-of-place parameter
    val = Ref(true)
    model = MockModel()
    p1 = Parameter{Bool}(() -> val[], model)
    @test p1()

    # test caching
    val[] = false
    @test p1()
    setdirty!(p1)
    @test !p1()

    # in-place parameter
    A = zeros(3, 4)
    Aval = Ref(1.0)
    f = let Aval = Aval
        m -> m .= Aval[]
    end
    p2 = Parameter(f, A, model)
    @test p2() === A
    @test all(p2() .== 1)

    # caching
    Aval[] = 2.0
    @test p2() === A
    @test all(p2() .== 1.0)
    setdirty!(p2)
    @test all(p2() .== 2.0)
    setdirty!(p2)
    allocs = @allocated p2()
    @test allocs == 0
end

end
