"""
$(SIGNATURES)

Sort the vector `v` in-place using `Base.sort!`, then combine all entries for which `by(v[i]) == by(v[j])`
using the function `v[i] = combine(v[i], v[j])`. This may result in `v` being resized to a shorter length.

See documentation for [`Base.sort!`](@ref) regarding the keyword arguments.
"""
function sort_and_combine!(v::AbstractVector;
        alg::Base.Sort.Algorithm=Base.Sort.defalg(v), combine, by=identity, lt=isless, rev::Bool=false, order::Base.Ordering=Base.Order.Forward)
    # For context, see: https://github.com/JuliaOpt/MathOptInterface.jl/issues/429#issuecomment-406232629.
    isempty(v) && return v
    sort!(v, alg=alg, lt=lt, by=by, rev=rev, order=order)
    j = 1
    @inbounds for i in eachindex(v)[2 : end]
        x = v[i]
        if lt(by(v[j]), by(x))
            j += 1
            v[j] = x
        else
            v[j] = combine(v[j], x)
        end
    end
    resize!(v, j)
    v
end
