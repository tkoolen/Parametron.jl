macro constraint(m, expr)
    @assert expr.head == :call
    @assert length(expr.args) == 3

    relation = expr.args[1]
    addcon = if relation == :>=
        :(SimpleQP.add_nonnegative_constraint!)
    elseif relation == :(==)
        :(SimpleQP.add_zero_constraint!)
    elseif relation == :<=
        :(SimpleQP.add_nonpositive_constraint!)
    else
        error("relation $(relation) not recognized.")
    end

    lhs = expr.args[2]
    rhs = expr.args[3]

    quote
        local f = SimpleQP.AffineFunction($(esc(lhs)) - $(esc(rhs)))
        $addcon($(esc(m)), f)
    end
end