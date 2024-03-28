
"""
    make_combinations(v...)

Make combinations of input vectors.

# Arguments
- `v...`: Any number of input `Vector`s.

# Returns
A `Vector` of `tuple`s with combinations of `Vector`s passed in `v...`.
If all vectors in `v...` have the same length, returns a vector of tuples of
length 1.
Combinations are returned in ascending order of the elements of each `i`-th
element of the input vectors.

# Example
```julia
vec1 = [1, 2, 3]
vec2 = ["a", "b"]

result = make_combinations(vec1, vec2)

# note order of results
println(result)
"""
function make_combinations(v...)
    # no input checking on this internal function
    # also no input checking as used to make combinations of
    # infection parameters and composable elements

    # recycle vectors and return a tuple
    vec(collect(Iterators.product(v...)))
end

export recycle_vectors
