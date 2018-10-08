module UtilTest

using Parametron
using Test

using Parametron: sort_and_combine!
using Base.Sort: QuickSort

combinepair(x, y) = first(x) => last(x) + last(y)

@testset "sort_and_combine!" begin
    n = 100
    @test sort_and_combine!([2 => 1.0, 2 => 2.0]; combine=combinepair, by=first) == [2 => 3.0]
    @test sort_and_combine!([3 => 1.0, 2 => 4.0, 3 => 2.0, 1 => 2.0], combine=combinepair, by=first) == [1 => 2.0, 2 => 4.0, 3 => 3.0]
    @test length(sort_and_combine!(Pair.(1 : n, rand(n)), combine=combinepair, by=first)) == n
    @test length(sort_and_combine!(Pair.(fill(3, n), rand(n)), combine=combinepair, by=first)) == 1

    v = [rand(1 : round(Int, n / 2)) => rand() for i = 1 : n]
    global allocs
    for i = 1 : 2
        vcopy = copy(v)
        allocs = @allocated sort_and_combine!(vcopy; combine=combinepair, by=first, alg=QuickSort)
    end
    @test allocs == 0
end

end # module
