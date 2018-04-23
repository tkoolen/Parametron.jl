function quad(A::Symmetric{S, <:SparseMatrixCSC}, x::Vector{T}) where {S, T}
    @boundscheck size(A, 1) == length(x) || error()
    ret = zero(promote_type(S, T))
    I, J, coeffs = findnz(A.data)
    uplo = A.uplo
    upper = uplo == 'U'
    @inbounds for k in eachindex(I)
        i = I[k]
        j = J[k]
        if (i <= j) == upper
            coeff = coeffs[k]
            val = coeff * x[i] * x[j]
            ret += i == j ? val : 2 * val
        end
    end
    ret
end
