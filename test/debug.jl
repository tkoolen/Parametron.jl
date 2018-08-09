module DebugTest

using Compat
using Compat.Test
using Parametron

import Parametron: MockModel, setdirty!, findallocs

@testset "findallocs" begin
    model = MockModel()
    x = Variable(model)
    p = Parameter{Int}(() -> 3, model)
    expr = @expression p * 4 + x
    setdirty!(model); findallocs(devnull, expr)
end

end
