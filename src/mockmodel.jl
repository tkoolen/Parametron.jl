import SimpleQP

"""
A 'mock model' used only for demonstrations and tests.
"""
struct MockModel
    params::Vector{SimpleQP.Parameter}
    nextvarindex::Base.RefValue{Int}
end
MockModel() = MockModel(SimpleQP.Parameter[], Ref(1))
SimpleQP.setdirty!(model::MockModel) = foreach(SimpleQP.setdirty!, model.params)
SimpleQP.addparameter!(model::MockModel, param::SimpleQP.Parameter) = (push!(model.params, param); param)
Variable(model::MockModel) = (ret = Variable(model.nextvarindex[]); model.nextvarindex[] += 1; ret);
