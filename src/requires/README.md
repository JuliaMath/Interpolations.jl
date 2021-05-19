This directory contains files that facilitate the composition of Interpolations.jl with other julia packages.

To add extra code that is loaded using Requries.jl create a new file named after the package of interest and then
modify `__init__` in Interpolations.jl with a `@require` statement to include the file.