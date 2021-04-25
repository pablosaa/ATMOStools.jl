using Documenter
using ATMOStools

makedocs(
    sitename = "ATMOStools",
    format = Documenter.HTML(),
    modules = [ATMOStools]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
