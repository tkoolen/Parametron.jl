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
