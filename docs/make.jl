import Literate
import Documenter
import DynamicSDDiP
import JuMP
import SDDP


const INPUT_DIR = joinpath(@__DIR__, "src")
const EXAMPLE_DIR = joinpath(@__DIR__, "src/example")

_sorted_files(dir, ext) = sort(filter(f -> endswith(f, ext), readdir(dir)))

function list_of_sorted_files(prefix, dir, ext = ".md")
    return Any["$(prefix)/$(file)" for file in _sorted_files(dir, ext)]
end

for filename_jl in list_of_sorted_files(INPUT_DIR, INPUT_DIR, ".jl")
    filename = replace(filename_jl, dirname(filename_jl) * "/" => "")

    Literate.markdown(
        filename_jl,
        INPUT_DIR;
        documenter = true,
    )

    # Literate.notebook(filename_jl, INPUT_DIR; execute = false, credit = false)
end

# for filename_jl in list_of_sorted_files(EXAMPLE_DIR, EXAMPLE_DIR, ".jl")
#     filename = replace(filename_jl, dirname(filename_jl) * "/" => "")

#     Literate.markdown(
#         filename_jl,
#         EXAMPLE_DIR;
#         documenter = true,
#     )

#     # Literate.notebook(filename_jl, EXAMPLE_DIR; execute = false, credit = false)
# end

# Documenter.makedocs(;
#     sitename = "DynamicSDDiP.jl",
#     authors = "Christian Fuellner",
#     clean = true,
#     pages = [
#         "Home" => "index.md",
#         "Code" => [
#             "model_assumptions.md",
#         ],
#         "Hydrothermal Example" => [
#             "example/experiment_description.md",
#         ],
#     ],
# )

# Documenter.deploydocs(;
#     repo = "github.com/ChrisFuelOR/DynamicSDDiP.jl.git",
#     push_preview = true,
# )