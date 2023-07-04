module Hello 

using Artifacts 

function __init__()
    path = joinpath(artifact"hello_artifact", "hello.txt")
    
    println(read(path,String))
end

end #module
