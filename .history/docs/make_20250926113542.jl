using Documenter
using MeanFieldToolbox

makedocs(
    build       =   "build" ,
    sitename    =   "MeanFieldToolbox.jl"    ,
    modules     =   [MeanFieldToolbox.MFTDecompose, MeanFieldToolbox.MFTBonds, MeanFieldToolbox.TBMFT, MeanFieldToolbox.BDGMFT, MeanFieldToolbox.Build, MeanFieldToolbox.MFTIter, MeanFieldToolbox.MFTRun, MeanFieldToolbox.MFTResume, MeanFieldToolbox.MFTPlot, MeanFieldToolbox.InteractionConvert]   ,
    pages = [
        "Introduction"              =>  "index.md",
        "MFTDecompose"              =>  "MFTDecompose.md",
        "MFTBonds"                  =>  "MFTBonds.md",
        "TightBindingMFT"           =>  "TightBindingMFT.md",
        "BdGMFT"                    =>  "BdGMFT.md",
        "Build"                     =>  "Build.md",
        "MFTIterator"               =>  "MFTIterator.md",
        "MFTRun"                    =>  "MFTRun.md",
        "MFTResume"                 =>  "MFTResume.md",
        "MFTPlot"                   =>  "MFTPlot.md",
        "InteractionConvert"        =>  "InteractionConvert.md"
    ]
)

deploydocs(
    repo = "github.com/andrewkhardy/MeanFieldToolbox.jl.git",
    devbranch = "main"
)