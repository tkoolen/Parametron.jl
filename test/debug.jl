module DebugTest

using Compat
using Compat.Test
using SimpleQP

import SimpleQP: MockModel, setdirty!, findallocs

@testset "findallocs" begin
    model = MockModel()
    p = Parameter(model) do
        3
    end
    expr = @expression p * 4
    setdirty!(model); findallocs(devnull, expr)
end

end
