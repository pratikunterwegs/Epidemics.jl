
# Getting started with Epidemics.jl

## Setting up the Environment

If want to use the Epidemics.jl package you need to install it first.
You can do it using the following commands:

```julia
julia> using Pkg

julia> Pkg.add("Epidemcs")
```

or

```julia
julia> ] # press ']' to enter the pkg command line

(@v1.8) pkg> add Epidemics
```

If you want to make sure everything works as expected you can run the tests
bundled with Epidemics.jl, but be warned that it will take more than 30
minutes:

```julia
julia> using Pkg

julia> Pkg.test("Epidemics") # Warning! This could take some time.
```

Additionally, it is recommended to check the version of Epidemics.jl that
you have installed with the `status` command.

```julia
julia> ]

(@v1.9) pkg> status Epidemics
```

```jldoctest epidemics
julia> using Epidemics
```

!!! note "Advanced installation configuration"

    Epidemics.jl related notice: some notice here

```@repl
a = 1
b = 2
a + b
```

```@example
using Epidemics
output = epidemic()

first(output, 5)
```
