var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#Parametron.jl-1",
    "page": "Home",
    "title": "Parametron.jl",
    "category": "section",
    "text": "Parametron makes it easy to set up and efficiently (ideally, with zero allocation) solve instances of a parameterized family of optimization problems.Parametron interfaces with various solvers through MathOptInterface.jl, similar to JuMP. However, unlike JuMP, the focus of Parametron is on efficient and user-friendly problem modification."
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": ""
},

{
    "location": "#Installing-Julia-1",
    "page": "Home",
    "title": "Installing Julia",
    "category": "section",
    "text": "Download links and more detailed instructions are available on the Julia website. The latest version of Parametron.jl requires Julia 0.7, but we recommend downloading 1.0 (the latest stable Julia release at the time of writing).warning: Warning\nDo not use apt-get or brew to install Julia, as the versions provided by these package managers tend to be out of date."
},

{
    "location": "#Installing-Parametron-1",
    "page": "Home",
    "title": "Installing Parametron",
    "category": "section",
    "text": "To install the latest tagged release of Parametron, start Julia and enter Pkg mode by pressing ]. Then simply runadd ParametronTo use the latest master version and work on the bleeding edge (generally, not recommended), instead runadd Parametron#masterA third option is to clone the repository (to the directory printed by julia -e \'import Pkg; println(Pkg.devdir())\'):dev Parametron"
},

{
    "location": "#Usage-1",
    "page": "Home",
    "title": "Usage",
    "category": "section",
    "text": "See the README for usage examples."
},

{
    "location": "#API-1",
    "page": "Home",
    "title": "API",
    "category": "section",
    "text": "See the API section for detailed documentation of exported functions."
},

{
    "location": "api/functions/#",
    "page": "Functions",
    "title": "Functions",
    "category": "page",
    "text": ""
},

{
    "location": "api/functions/#Parametron.Functions",
    "page": "Functions",
    "title": "Parametron.Functions",
    "category": "module",
    "text": "The Functions module provides types that represent decision variables and functions of these decision variables that possess a certain structure, such as being affine or quadratic in the decision variables.\n\nVarious operators are overloaded to make it easy to construct such functions.\n\nExample\n\njulia> x, y = Variable.(1 : 2)\n2-element Array{Parametron.Functions.Variable,1}:\n Parametron.Functions.Variable(1)\n Parametron.Functions.Variable(2)\n\njulia> fun = 3 * (x + y) + 4\n3 * x1 + 3 * x2 + 4\n\njulia> typeof(fun)\nParametron.Functions.AffineFunction{Int64}\n\njulia> fun.linear\n2-element Array{Parametron.Functions.LinearTerm{Int64},1}:\n\n 3 * x1\n 3 * x2\n\njulia> fun.constant\nBase.RefValue{Int64}(4)\n\nThe Functions module also provides in-place versions of certain common operations, e.g., add!, subtract!, and vecdot!, which may be used to evaluate operations performed on the functions into a preallocated destination without any heap allocations.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Functions-1",
    "page": "Functions",
    "title": "Functions",
    "category": "section",
    "text": "CurrentModule = Parametron.FunctionsFunctions"
},

{
    "location": "api/functions/#Parametron.Functions.Variable",
    "page": "Functions",
    "title": "Parametron.Functions.Variable",
    "category": "type",
    "text": "struct Variable\n\nRepresents a single decision variable.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.LinearTerm",
    "page": "Functions",
    "title": "Parametron.Functions.LinearTerm",
    "category": "type",
    "text": "struct LinearTerm{T}\n\nRepresents a scalar linear term, i.e. a decision variable scaled by a coefficient.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.QuadraticTerm",
    "page": "Functions",
    "title": "Parametron.Functions.QuadraticTerm",
    "category": "type",
    "text": "struct QuadraticTerm{T}\n\nRepresents a scalar quadratic term, i.e. the product of two decision variables scaled by a coefficient.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.AffineFunction",
    "page": "Functions",
    "title": "Parametron.Functions.AffineFunction",
    "category": "type",
    "text": "struct AffineFunction{T}\n\nA scalar affine function represented by a sum of LinearTerms and a constant.\n\nAffineFunction overloads the call operator, which can be used to evalute the function given values for the decision variables. The call operator takes an AbstractDict{Variable, T} collection, which associates values with variables.\n\nExample\n\njulia> x, y = Variable.(1 : 2);\n\njulia> affinefun = 2 * x + 3 * y + 4\n2 * x1 + 3 * x2 + 4\n\njulia> vals = Dict(x => 4, y => -3);\n\njulia> affinefun(vals)\n3\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.QuadraticFunction",
    "page": "Functions",
    "title": "Parametron.Functions.QuadraticFunction",
    "category": "type",
    "text": "struct QuadraticFunction{T}\n\nA scalar quadratic function represented by a sum of QuadraticTerms and an AffineFunction.\n\nQuadraticFunction overloads the call operator, which can be used to evalute the function given values for the decision variables. The call operator takes an AbstractDict{Variable, T} collection, which associates values with variables.\n\nExamples\n\njulia> x, y = Variable.(1 : 2);\n\njulia> quadraticfun = 2 * x^2 + 3 * x * y - 2 * y + 4\n2 * x1 * x1 + 3 * x1 * x2 + -2 * x2 + 4\n\njulia> vals = Dict(x => 4, y => -3);\n\njulia> quadraticfun(vals)\n6\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Types-1",
    "page": "Functions",
    "title": "Types",
    "category": "section",
    "text": "Variable\nLinearTerm\nQuadraticTerm\nAffineFunction\nQuadraticFunction"
},

{
    "location": "api/functions/#Parametron.Functions.canonicalize",
    "page": "Functions",
    "title": "Parametron.Functions.canonicalize",
    "category": "function",
    "text": "Re-express a term or function in a canonical form.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.canonicalize-Tuple{QuadraticTerm}",
    "page": "Functions",
    "title": "Parametron.Functions.canonicalize",
    "category": "method",
    "text": "canonicalize(term)\n\n\nRe-express the QuadraticTerm term so that the index of term.rowvar is less than or equal to that of term.colvar.\n\nExample\n\njulia> x, y = Variable.(1 : 2);\n\njulia> term = 3 * y * x\n3 * x2 * x1\n\njulia> canonicalize(term)\n3 * x1 * x2\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.canonicalize-Tuple{AffineFunction}",
    "page": "Functions",
    "title": "Parametron.Functions.canonicalize",
    "category": "method",
    "text": "canonicalize(f)\n\n\nReturn a canonicalized version of f::AffineFunction, namely with linear terms sorted by variable index and with terms corresponding to the same variable combined.\n\nExample\n\njulia> x, y = Variable.(1 : 2);\n\njulia> f = y + x - 2 * y + 3\n1 * x2 + 1 * x1 + -2 * x2 + 3\n\njulia> canonicalize(f)\n1 * x1 + -1 * x2 + 3\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.canonicalize-Tuple{QuadraticFunction}",
    "page": "Functions",
    "title": "Parametron.Functions.canonicalize",
    "category": "method",
    "text": "canonicalize(f)\n\n\nReturn a canonicalized version of f::QuadraticFunction. See canonicalize(f::QuadraticTerm) and canonicalize(f::AffineFunction) for more details. Quadratic terms are ordered lexicographically by (term.rowvar, term.colvar).\n\nExample\n\njulia> x, y = Variable.(1 : 2);\n\njulia> f = x * y + y * x + y + y + x - y\n1 * x1 * x2 + 1 * x2 * x1 + 1 * x2 + 1 * x2 + 1 * x1 + -1 * x2 + 0\n\njulia> canonicalize(f)\n2 * x1 * x2 + 1 * x1 + 1 * x2 + 0\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.canonicalize!",
    "page": "Functions",
    "title": "Parametron.Functions.canonicalize!",
    "category": "function",
    "text": "In-place version of canonicalize.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.prune_zero",
    "page": "Functions",
    "title": "Parametron.Functions.prune_zero",
    "category": "function",
    "text": "prune_zero(f; kwargs...)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:299.\n\nprune_zero(f; kwargs...)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:415.\n\nPrune terms with (approximately) zero coefficients.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.prune_zero!",
    "page": "Functions",
    "title": "Parametron.Functions.prune_zero!",
    "category": "function",
    "text": "prune_zero!(f; atol)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:295.\n\nprune_zero!(f; atol)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:410.\n\nIn-place version of prune_zero.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Exported-functions-1",
    "page": "Functions",
    "title": "Exported functions",
    "category": "section",
    "text": "canonicalize\ncanonicalize(::QuadraticTerm)\ncanonicalize(::AffineFunction)\ncanonicalize(::QuadraticFunction)\ncanonicalize!prune_zero\nprune_zero!"
},

{
    "location": "api/functions/#Parametron.Functions.add!",
    "page": "Functions",
    "title": "Parametron.Functions.add!",
    "category": "function",
    "text": "Add x to f, modifying f.\n\nadd!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:452.\n\nadd!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:453.\n\nadd!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:454.\n\nadd!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:455.\n\nadd!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:457.\n\nadd!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:458.\n\nadd!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:459.\n\nadd!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:461.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.subtract!",
    "page": "Functions",
    "title": "Parametron.Functions.subtract!",
    "category": "function",
    "text": "Subtract x from f, modifying f.\n\nsubtract!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:474.\n\nsubtract!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:475.\n\nsubtract!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:476.\n\nsubtract!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:478.\n\nsubtract!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:487.\n\nsubtract!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:488.\n\nsubtract!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:489.\n\nsubtract!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:490.\n\nsubtract!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:491.\n\nsubtract!(f, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:493.\n\nsubtract!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:502.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.muladd!",
    "page": "Functions",
    "title": "Parametron.Functions.muladd!",
    "category": "function",
    "text": "Multiply x by y and add the result to dest.\n\nmuladd!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:516.\n\nmuladd!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:524.\n\nmuladd!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:527.\n\nmuladd!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:535.\n\nmuladd!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:538.\n\nmuladd!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:546.\n\nmuladd!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:549.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.vecdot!",
    "page": "Functions",
    "title": "Parametron.Functions.vecdot!",
    "category": "function",
    "text": "Take the dot product of vectors x and y, storing the result in dest.\n\nvecdot!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:713.\n\nvecdot!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:716.\n\nvecdot!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:718.\n\nvecdot!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:720.\n\nvecdot!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:722.\n\nvecdot!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:724.\n\nvecdot!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:726.\n\nvecdot!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:728.\n\nvecdot!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:730.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.vecadd!",
    "page": "Functions",
    "title": "Parametron.Functions.vecadd!",
    "category": "function",
    "text": "Add vector x to vector y, storing the result in dest.\n\nvecadd!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:754.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.vecsubtract!",
    "page": "Functions",
    "title": "Parametron.Functions.vecsubtract!",
    "category": "function",
    "text": "Subtract vector y from vector x, storing the result in dest.\n\nvecsubtract!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:754.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.matvecmul!",
    "page": "Functions",
    "title": "Parametron.Functions.matvecmul!",
    "category": "function",
    "text": "Compute the matrix-vector product A * x, storing the result in y.\n\nmatvecmul!(y, A, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:779.\n\nmatvecmul!(y, A, x)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:804.\n\nmatvecmul!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:828.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.bilinearmul!",
    "page": "Functions",
    "title": "Parametron.Functions.bilinearmul!",
    "category": "function",
    "text": "Compute the bilinear form transpose(x) * A * y, storing the result in dest.\n\nbilinearmul!(dest, Q, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:845.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.scale!",
    "page": "Functions",
    "title": "Parametron.Functions.scale!",
    "category": "function",
    "text": "Scale a vector by a number and store the result in dest.\n\nscale!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:877.\n\nscale!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:888.\n\nscale!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:899.\n\nscale!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:910.\n\nscale!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:918.\n\nscale!(dest, x, y)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:923.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#Parametron.Functions.vcat!",
    "page": "Functions",
    "title": "Parametron.Functions.vcat!",
    "category": "function",
    "text": "Vertically concatenate a number of vectors, storing the result in y.\n\nvcat!(y, args)\n\ndefined at /home/travis/build/tkoolen/Parametron.jl/src/functions.jl:989.\n\n\n\n\n\n"
},

{
    "location": "api/functions/#In-place-math-functions-(unexported)-1",
    "page": "Functions",
    "title": "In-place math functions (unexported)",
    "category": "section",
    "text": "add!\nsubtract!\nmuladd!\nvecdot!\nvecadd!\nvecsubtract!\nmatvecmul!\nbilinearmul!\nscale!\nvcat!"
},

{
    "location": "api/parameter/#",
    "page": "Parameter",
    "title": "Parameter",
    "category": "page",
    "text": ""
},

{
    "location": "api/parameter/#Parametron.Parameter",
    "page": "Parameter",
    "title": "Parametron.Parameter",
    "category": "type",
    "text": "struct Parameter{T, F, InPlace}\n\nRepresents a \'placeholder\' for data; a value that may be filled in later.\n\nParameters can be evaluated by simply calling them with no arguments.\n\nParameters keep track of whether they\'ve already been evaluated using dirty flag. To reevalute a parameter, the dirty flag must first be set using setdirty!(parameter). The update function will then be called when the parameter itself is called.\n\nExamples\n\njulia> value = Ref(1)\nBase.RefValue{Int64}(1)\n\njulia> p = Parameter{Int}(() -> value[], model)\nParameter{Int64, …}(…)\n\njulia> p()\n1\n\njulia> value[] = 2\n2\n\njulia> p()\n1\n\njulia> Parametron.setdirty!(p); p()\n2\n\n\n\n\n\n"
},

{
    "location": "api/parameter/#Parametron.Parameter-Tuple{Any,Any}",
    "page": "Parameter",
    "title": "Parametron.Parameter",
    "category": "method",
    "text": "Parameter(f, model)\n\n\nCreate a new \'out-of-place\' Parameter with an update function f that takes no arguments. The type of the output is determined upon construction by calling the update function.\n\nwarning: Warning\nExplicitly specifying the return value type using the Parameter{T}(f, model) constructor is preferred, as using this constructor can lead to type inference issues.\n\n\n\n\n\n"
},

{
    "location": "api/parameter/#Parametron.Parameter-Union{Tuple{Any}, Tuple{T}} where T",
    "page": "Parameter",
    "title": "Parametron.Parameter",
    "category": "method",
    "text": "Parameter(model; val)\n\n\nCreate a new \'in-place\' Parameter that always returns val. This constructor may be used to create Parameters that use val as a work buffer that is manually/externally updated.\n\nwarning: Warning\nBy using this constructor, the automated mechanism for lazily updating a Parameter\'s value when necessary is essentially circumvented, and the user is responsible for ensuring that the Parameter\'s value is updated at the appropriate time.\n\n\n\n\n\n"
},

{
    "location": "api/parameter/#Parametron.Parameter-Union{Tuple{F}, Tuple{T}, Tuple{F,Any}} where F where T",
    "page": "Parameter",
    "title": "Parametron.Parameter",
    "category": "method",
    "text": "Parameter(f, model)\nParameter(model; val)\n\n\nCreate a new \'out-of-place\' Parameter with an update function f that takes no arguments and returns a value of type T.\n\n\n\n\n\n"
},

{
    "location": "api/parameter/#Parametron.Parameter-Union{Tuple{F}, Tuple{T}, Tuple{F,T,Any}} where F where T",
    "page": "Parameter",
    "title": "Parametron.Parameter",
    "category": "method",
    "text": "Parameter(f, val, model)\nParameter(model; val)\n\n\nCreate a new \'in-place\' Parameter with an update function f that takes val as its argument and updates it in place.\n\n\n\n\n\n"
},

{
    "location": "api/parameter/#Parameter-1",
    "page": "Parameter",
    "title": "Parameter",
    "category": "section",
    "text": "Modules = [Parametron]\nPages   = [\"parameter.jl\"]"
},

{
    "location": "api/lazyexpression/#",
    "page": "LazyExpression",
    "title": "LazyExpression",
    "category": "page",
    "text": ""
},

{
    "location": "api/lazyexpression/#LazyExpression-1",
    "page": "LazyExpression",
    "title": "LazyExpression",
    "category": "section",
    "text": ""
},

{
    "location": "api/lazyexpression/#Parametron.@expression",
    "page": "LazyExpression",
    "title": "Parametron.@expression",
    "category": "macro",
    "text": "Create a new LazyExpression and apply optimizations to it to reduce allocations and improve performance.\n\nExpressions that do not depend on Parameters or other [LazyExpression]s are simply evaluated straight away.\n\nExamples\n\nCreating an expression that represents p * x1, where p is a parameter that always evaluates to 2:\n\njulia> x1 = Variable(model);\n\njulia> p = Parameter{Int}(() -> 2, model)\nParameter{Int64, …}(…)\n\njulia> expr = @expression p * x1\nLazyExpression{FunctionWrapper{…}(LazyExpression{typeof(*), …}(…))}(…)\n\njulia> expr()\n2 * x1\n\nCreating an expression that represents p ⋅ x, where p is a parameter that evaluates to [1, 2] and x is a vector of two variables:\n\njulia> model = Parametron.mock_model();\n\njulia> x = Variable.(1 : 2);\n\njulia> p = Parameter(model, val=[1, 2])\nParameter{Array{Int64,1}, …}(…)\n\njulia> expr = @expression p ⋅ x\nLazyExpression{FunctionWrapper{…}(LazyExpression{typeof(Parametron.Functions.vecdot!), …}(…))}(…)\n\njulia> expr()\n1 * x1 + 2 * x2 + 0\n\njulia> @allocated expr()\n0\n\nNote that evaluating the expression does not allocate, because the ⋅ operation is optimized and transformed into a call to the in-place Functions.vecdot! function.\n\n\n\n\n\n"
},

{
    "location": "api/lazyexpression/#The-@expression-macro-1",
    "page": "LazyExpression",
    "title": "The @expression macro",
    "category": "section",
    "text": "@expression"
},

{
    "location": "api/lazyexpression/#Parametron.LazyExpression",
    "page": "LazyExpression",
    "title": "Parametron.LazyExpression",
    "category": "type",
    "text": "struct LazyExpression{F, A}\n\nRepresents an expression that may be evaluated at a later time, by storing both a function, f, and a tuple of function arguments, args.\n\nLazyExpressions are typically not manually constructed by a user, and hence are not exported. Instead, LazyExpressions should be created using the [@expression] macro.\n\nA LazyExpression may be evaluated by simply calling it with no arguments.\n\nExample\n\njulia> a = ones(2); b = ones(2);\n\njulia> expr = Parametron.LazyExpression(+, a, b)\nLazyExpression{typeof(+), …}(…)\n\njulia> expr()\n2-element Array{Float64,1}:\n 2.0\n 2.0\n\njulia> b .= 2\n2-element Array{Float64,1}:\n 2.0\n 2.0\n\njulia> expr()\n2-element Array{Float64,1}:\n 3.0\n 3.0\n\n\n\n\n\n"
},

{
    "location": "api/lazyexpression/#The-LazyExpression-type-1",
    "page": "LazyExpression",
    "title": "The LazyExpression type",
    "category": "section",
    "text": "Parametron.LazyExpression"
},

{
    "location": "api/lazyexpression/#Parametron.wrap",
    "page": "LazyExpression",
    "title": "Parametron.wrap",
    "category": "function",
    "text": "wrap(expr)\n\n\nWrap a LazyExpression in a FunctionWrappers.FunctionWrapper and return a new LazyExpression with the FunctionWrapper as the function f and an empty tuple as the arguments arg.\n\nThe type parameters of the returned LazyFunction depend only on the type of the value returned by expr. This is useful when a common interface is needed for different LazyExpressions that share the same return value type.\n\n\n\n\n\n"
},

{
    "location": "api/lazyexpression/#Wrapping-LazyExpressions-1",
    "page": "LazyExpression",
    "title": "Wrapping LazyExpressions",
    "category": "section",
    "text": "Parametron.wrap"
},

{
    "location": "api/model/#",
    "page": "Model",
    "title": "Model",
    "category": "page",
    "text": ""
},

{
    "location": "api/model/#Model-1",
    "page": "Model",
    "title": "Model",
    "category": "section",
    "text": ""
},

{
    "location": "api/model/#Parametron.Model",
    "page": "Model",
    "title": "Parametron.Model",
    "category": "type",
    "text": "Model(optimizer)\n\n\nCreate a new Model, representing an optimization problem to be solved by the optimizer optimizer (a MathOptInterface.AbstractOptimizer).\n\n\n\n\n\n"
},

{
    "location": "api/model/#Creating-a-Model-1",
    "page": "Model",
    "title": "Creating a Model",
    "category": "section",
    "text": "Model"
},

{
    "location": "api/model/#Parametron.Functions.Variable-Tuple{Model}",
    "page": "Model",
    "title": "Parametron.Functions.Variable",
    "category": "method",
    "text": "Variable(m)\n\n\nCreate a new decision variable (Variable) associated with the model.\n\n\n\n\n\n"
},

{
    "location": "api/model/#Creating-Decision-variables-1",
    "page": "Model",
    "title": "Creating Decision variables",
    "category": "section",
    "text": "Variable(::Model)"
},

{
    "location": "api/model/#Parametron.@constraint",
    "page": "Model",
    "title": "Parametron.@constraint",
    "category": "macro",
    "text": "Add a constraint to the model using operators ==, <=, >=, or in/∈.\n\nin/∈ may only be used for single variables with a right hand side that is one of:\n\nℤ or Integers\n{0, 1} or ZeroOne\n\nExamples\n\nThe constraint x >= zeros(2) can be added to a model as follows:\n\njulia> x = [Variable(model) for i = 1 : 2];\n\njulia> @constraint(model, x >= zeros(2))\n\nThe constraint that variable x[1] should be an integer can be expressed using:\n\njulia> @constraint(model, x[1] ∈ ℤ)\n\n\n\n\n\n"
},

{
    "location": "api/model/#Parametron.@objective",
    "page": "Model",
    "title": "Parametron.@objective",
    "category": "macro",
    "text": "Set the objective function of the model.\n\nExamples\n\nLet model be a Model instance. The objective \'minimize x ⋅ x\' can be added as follows:\n\njulia> x = [Variable(model) for i = 1 : 2];\n\njulia> @objective model Minimize x ⋅ x;\n\n\n\n\n\n"
},

{
    "location": "api/model/#Parametron.setobjective!",
    "page": "Model",
    "title": "Parametron.setobjective!",
    "category": "function",
    "text": "setobjective!(m, sense, expr)\n\n\nSet the objective function and optimization sense (Minimize or Maximize).\n\n\n\n\n\n"
},

{
    "location": "api/model/#Adding-constraints-and-an-objective-function-1",
    "page": "Model",
    "title": "Adding constraints and an objective function",
    "category": "section",
    "text": "@constraint\n@objective\nsetobjective!"
},

{
    "location": "api/model/#Parametron.solve!",
    "page": "Model",
    "title": "Parametron.solve!",
    "category": "function",
    "text": "solve!(m)\n\n\nSolve the model m. (Re-)evaluate constraint and objective expressions, update the optimizer\'s internal representation of the problem, and start the optimization procedure.\n\n\n\n\n\n"
},

{
    "location": "api/model/#Parametron.setdirty!",
    "page": "Model",
    "title": "Parametron.setdirty!",
    "category": "function",
    "text": "setdirty!(model)\n\n\nMark all parameters associated with the model as \'dirty\' (out of date), meaning they must be updated upon their next evaluation.\n\n\n\n\n\n"
},

{
    "location": "api/model/#Parametron.initialize!",
    "page": "Model",
    "title": "Parametron.initialize!",
    "category": "function",
    "text": "initialize!(m)\n\n\nCopy the problem to be solved to the optimizer.\n\nUsers should generally not need to call this function directly, as it is automatically called the first time solve! is called on a Model.\n\n\n\n\n\n"
},

{
    "location": "api/model/#Parametron.update!",
    "page": "Model",
    "title": "Parametron.update!",
    "category": "function",
    "text": "Re-evaluate the expressions used to build the constraints and objective function of Model m.\n\nUsers should generally not need to call this function directly, as it is automatically called in solve!.\n\n\n\n\n\n"
},

{
    "location": "api/model/#Solving-1",
    "page": "Model",
    "title": "Solving",
    "category": "section",
    "text": "solve!\nParametron.setdirty!\nParametron.initialize!\nParametron.update!"
},

{
    "location": "api/model/#Parametron.value",
    "page": "Model",
    "title": "Parametron.value",
    "category": "function",
    "text": "value(m, x)\n\n\nReturn the value of variable x as determined by the optimizer.\n\n\n\n\n\n"
},

{
    "location": "api/model/#Parametron.objectivevalue",
    "page": "Model",
    "title": "Parametron.objectivevalue",
    "category": "function",
    "text": "objectivevalue(m)\n\n\nReturn the value of the objective function at the solution found by the optimizer.\n\n\n\n\n\n"
},

{
    "location": "api/model/#Parametron.terminationstatus",
    "page": "Model",
    "title": "Parametron.terminationstatus",
    "category": "function",
    "text": "terminationstatus(m)\n\n\nReturn the termination status of the solver.\n\n\n\n\n\n"
},

{
    "location": "api/model/#Parametron.primalstatus",
    "page": "Model",
    "title": "Parametron.primalstatus",
    "category": "function",
    "text": "primalstatus(m)\n\n\nReturn information regarding the primal of the problem.\n\n\n\n\n\n"
},

{
    "location": "api/model/#Parametron.dualstatus",
    "page": "Model",
    "title": "Parametron.dualstatus",
    "category": "function",
    "text": "dualstatus(m)\n\n\nReturn information regarding the dual of the problem.\n\n\n\n\n\n"
},

{
    "location": "api/model/#Accessing-solver-results-1",
    "page": "Model",
    "title": "Accessing solver results",
    "category": "section",
    "text": "value\nobjectivevalue\nterminationstatus\nprimalstatus\ndualstatus"
},

{
    "location": "api/debug/#",
    "page": "Debugging Utilities",
    "title": "Debugging Utilities",
    "category": "page",
    "text": ""
},

{
    "location": "api/debug/#Parametron.findallocs",
    "page": "Debugging Utilities",
    "title": "Parametron.findallocs",
    "category": "function",
    "text": "findallocs(x)\n\n\nUtility function that can be used to track down allocations in LazyExpressions.\n\nExamples\n\nThe following session shows the output of findallocs if the expression doesn\'t allocate:\n\njulia> x = [Variable(model) for i in 1 : 2];\n\njulia> param = Parameter{Int}(() -> 2, model)\nParameter{Int64, …}(…)\n\njulia> expr = @expression param * x\nLazyExpression{FunctionWrapper{…}(LazyExpression{typeof(Parametron.Functions.scale!), …}(…))}(…)\n\njulia> Parametron.findallocs(expr)\nLazyExpression{FunctionWrapper{…}(LazyExpression{typeof(Parametron.Functions.scale!), …}(…))}(…): 0 bytes\n  [1]: Array{LinearTerm{Int64},1}\n  [2]: Parameter{Int64, …}(…): 0 bytes\n  [3]: Array{Variable,1}\n\nIn this session, param allocates, and findallocs reports the allocation:\n\njulia> x = [Variable(model) for i in 1 : 2];\n\njulia> param = Parameter(() -> zeros(2), model)\nParameter{Array{Float64,1}, …}(…)\n\njulia> expr = @expression param ⋅ x\nLazyExpression{FunctionWrapper{…}(LazyExpression{typeof(Parametron.Functions.vecdot!), …}(…))}(…)\n\njulia> Parametron.findallocs(expr)\nLazyExpression{FunctionWrapper{…}(LazyExpression{typeof(Parametron.Functions.vecdot!), …}(…))}(…): 0 bytes\n  [1]: AffineFunction{Float64}\n  [2]: Parameter{Array{Float64,1}, …}(…): 96 bytes\n  [3]: Array{Variable,1}\n\n\n\n\n\n"
},

{
    "location": "api/debug/#Debugging-Utilities-1",
    "page": "Debugging Utilities",
    "title": "Debugging Utilities",
    "category": "section",
    "text": "findallocs"
},

]}
