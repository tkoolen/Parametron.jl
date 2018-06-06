module SimpleQPTest

using SimpleQP

import SimpleQP: setdirty!

struct MockModel
    params::Vector{Parameter}
end
MockModel() = MockModel(Parameter[])
SimpleQP.setdirty!(model::MockModel) = foreach(setdirty!, model.params)
SimpleQP.addparameter!(model::MockModel, param::Parameter) = (push!(model.params, param); param)

include("parameter.jl")
include("functions.jl")
include("lazyexpression.jl")
include("model.jl")

end # module
