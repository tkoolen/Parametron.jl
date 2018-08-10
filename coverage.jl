import Pkg
Pkg.add("Coverage")
using Coverage
cd(@__DIR__) do
    Codecov.submit(Codecov.process_folder())
end
