"""
$(SIGNATURES)

Sort the vector `v` in-place using `Base.sort!`, then combine all entries which compare equal
(as determined using the `lt` and `by` functions), using the function `combined = combine(x, y)`.

See documentation for [`Base.sort!`](@ref) regarding the `alg`, `by` and `lt` keyword arguments.
"""
function sort_and_combine!(v::AbstractVector;
        alg::Base.Sort.Algorithm=Base.Sort.defalg(v), combine, by=identity, lt=isless)
    # For context, see: https://github.com/JuliaOpt/MathOptInterface.jl/issues/429#issuecomment-406232629.
    isempty(v) && return v
    sort!(v, alg=alg, lt=lt, by=by)
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
