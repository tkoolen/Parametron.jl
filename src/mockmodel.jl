struct MockModel
    params::Vector{Parametron.Parameter}
    nextvarindex::Base.RefValue{Int}
end
MockModel() = MockModel(Parametron.Parameter[], Ref(1))
Parametron.setdirty!(model::MockModel) = foreach(Parametron.setdirty!, model.params)
Parametron.addparameter!(model::MockModel, param::Parametron.Parameter) = (push!(model.params, param); param)
Variable(model::MockModel) = (ret = Variable(model.nextvarindex[]); model.nextvarindex[] += 1; ret);

function mock_model()
    optimizer = MOIU.MockOptimizer(Parametron.ParametronMOIModel{Float64}())
    Parametron.Model(optimizer)
end
