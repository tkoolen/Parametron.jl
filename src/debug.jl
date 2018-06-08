function findallocs(io::IO, x, depth = 0, argnum = nothing)
    depth > 0 && print(io, "  "^depth)
    argnum != nothing && print(io, "[$argnum]: ")
    if x isa LazyExpression || x isa Parameter
        x isa Parameter && setdirty!(x)
        allocs = @allocated x()
        print(io, "$x: ")
        color = allocs > 0 ? :light_red : :green
        print_with_color(color, io, allocs)
        print(io, " bytes")
        println(io)
    else
        println(io, typeof(x))
    end
    if x isa LazyExpression
        for (argnum, arg) in enumerate(x.args)
            findallocs(io, arg, depth + 1, argnum)
        end
    end
end

findallocs(x) = findallocs(stdout, x)
