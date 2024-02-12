
"""
    recycle_vectors(v...)

Recycles input vectors to match the length of the longest vector following
Tidyverse rules.

# Arguments
- `v...`: Any number of input vectors. All vectors must be of the same length,
or all but one vector must have length 1.

# Returns
An array of vectors passed in `v...`. If all vectors in `v...` have the same
length, returns a vector of these vectors unmodified. If vectors must be
recycled, returns all vectors repeated for the same number of times as the
longest vector.

# Example
```julia
vec1 = [1, 2, 3]
vec2 = [4]

result = recycle_vectors(vec1, vec2)
println(result)
"""
function recycle_vectors(v...)
    # convert from tuple to array
    v = [i for i in v]
    vec_lengths = length.(v)

    # return unmodified if all have same length
    if length(unique(vec_lengths)) == 1
        return v
    end

    longest = maximum(vec_lengths)
    # prefer not to use modification in place with `map!()`
    v[vec_lengths .!= longest] = map(x -> repeat([x], longest), v[vec_lengths .!= longest])
    v
end

export recycle_vectors
