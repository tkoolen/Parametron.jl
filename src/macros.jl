# mul_extract_linear_terms(x, y, z...) = mul_extract_linear_terms(mul_extract_linear_terms(x, y), z...)
# mul_extract_linear_terms(x, y) = x * y
# mul_extract_linear_terms(x::Matrix{Float64}, y::Vector{Variable}) = LinearTerm(x, y)

# maybe_constant(x) = x
# maybe_constant(x::Vector{Float64}) = Constant(x)

# function process_expression(expr)
#     expr = postwalk(expr) do x # turn Vector{Float64} into Constant
#         if @capture(x, f_(args__))
#             args = map(arg -> arg isa Symbol ? :(SimpleQP.maybe_constant($arg)) : arg, args)
#             :($f($(args...)))
#         else
#             x
#         end
#     end
#     expr
# end

macro constraint(m, expr)
    addcon = if @capture(expr, >=(lhs_, rhs_))
        :(SimpleQP.add_nonnegative_constraint!)
    elseif @capture(expr, ==(lhs_, rhs_))
        :(SimpleQP.add_zero_constraint!)
    elseif @capture(expr, <=(lhs_, rhs_))
        :(SimpleQP.add_nonpositive_constraint!)
    else
        throw(ArgumentError("Relation not recognized"))
    end
    quote
        f = convert(SimpleQP.AffineFunction, $(esc(lhs))) - convert(SimpleQP.AffineFunction, $(esc(rhs)))
        $addcon($(esc(m)), f)
    end
end
