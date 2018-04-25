wrapconstant(x) = x
wrapconstant(x::Vector{Float64}) = Constant(x)

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
    lhs = postwalk(x -> x isa Symbol ? :(SimpleQP.wrapconstant($x)) : x, lhs)
    rhs = postwalk(x -> x isa Symbol ? :(SimpleQP.wrapconstant($x)) : x, rhs)
    quote
        f = convert(SimpleQP.AffineFunction, $(esc(lhs)) - $(esc(rhs)))
        $addcon($(esc(m)), f)
    end
end
