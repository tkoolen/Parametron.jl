import SimpleQP

struct MockModel
    params::Vector{SimpleQP.Parameter}
end
MockModel() = MockModel(SimpleQP.Parameter[])
SimpleQP.setdirty!(model::MockModel) = foreach(SimpleQP.setdirty!, model.params)
SimpleQP.addparameter!(model::MockModel, param::SimpleQP.Parameter) = (push!(model.params, param); param)
