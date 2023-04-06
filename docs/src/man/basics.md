
# Getting started with Epidemics.jl

!!! note "Experimental package"

    Epidemics.jl is not available from the package registry, but the latest development version can be installed from Github.

```julia
julia> using Pkg

julia> Pkg.add("https://github.com/pratikunterwegs/Epidemics.jl")
```

The version of Epidemics.jl that is installed can be checked with the `status` command.

```julia
julia> ]

(@v1.8) pkg> status Epidemics
```

## Basic usage

```@example
using Epidemics
output = epidemic()

first(output, 5)
```
