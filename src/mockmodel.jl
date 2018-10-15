Base.@deprecate MockModel() mock_model()

function mock_model()
    optimizer = MOIU.MockOptimizer(Parametron.ParametronMOIModel{Float64}())
    Parametron.Model(optimizer)
end
