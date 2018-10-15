module DebugTest

using Test
using Parametron

import Parametron: mock_model, setdirty!, findallocs

@testset "findallocs" begin
    model = mock_model()
    x = Variable(model)
    p = Parameter{Int}(() -> 3, model)
    expr = @expression p * 4 + x
    setdirty!(model); findallocs(devnull, expr)
end

end
