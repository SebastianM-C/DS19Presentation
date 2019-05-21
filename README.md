# DS19Presentation
Materials for the SIAM DS19 Presentation
![conference_logo](assets/DS19_logo.png)

# Julia instructions

In order to reproduce the results you need to reproduce the environment in which
they were created. For this task, julia projects use `Manifest.toml` files. Thus,
the package manager will install the exact same package versions as the ones
that were used by the author.
Moreover, the main dataset was uploaded to figshare and with the help of the
DataDeps package it will be automatically downloaded on the first use.

After cloning the repository and changing to the `DS19Presentation` directory
open a julia REPL in that directory and follow the instructions bellow:

1. Enter `Pkg` mode and activate the project. In order to enter the `Pkg` mode press `]`
at the julia prompt. The prompt will change from `julia>` to `(v1.1) pkg>`.

```julia
(v1.1) pkg> activate .
```

2. Download the required packages

```julia
(DS19Presentation) pkg> instantiate

julia> using DS19Presentation
[ Info: Recompiling stale cache file C:\Users\sebastian\.julia\compiled\v1.1\DS19Presentation\aQbRc.ji for DS19Presentation [c52eaf16-6d25-11e9-34da-d5c7cd6dfd33]WARNING: using GLAbstraction.update! in module GLMakie conflicts with an existing identifier.

Update message: DynamicalSystems v1.3

A method that estimates the predictability properties
of a
dynamical system has been implemented, following the work of:

Wernecke, H., Sándor, B. & Gros, C.
*How to test for partially predictable chaos*.

See the function `predictability`.

WARNING: could not import Base.quit into AtomShell

julia> WARNING: both AbstractPlotting and DynamicalSystems export "step!"; uses of it in module InteractiveChaos must be qualified
WARNING: both CSSUtil and Base export "empty"; uses of it in module Interact must be qualified
```

3. Load the main dataset. The first time you will do this, you will be asked if
you want to download it. To exit the `Pkg` mode press `backspace`.

```julia
julia> g=load()
This program has requested access to the data dependency DS19 presentation materials test.
which is not currently installed. It can be installed
automatically, and you will not see this message again.

Dataset: DS19 presentation materials test
Website: https://figshare.com/articles/DS19_test_database/8078699
Author: Sebastian Micluța-Câmpeanu
Date of Publication: 2019-05-14T10:46:58Z
License: MIT (https://opensource.org/licenses/MIT)

This is the database containing the data to reproduce
the plots in the DS19 presentation of the author. This is just a test version.

Please cite this dataset: Micluța-Câmpeanu, Sebastian
(2019): DS19 presentation materials test. figshare. Conference contribution.



Do you want to download the dataset from Any["https://ndownloader.figshare.com/files/15142844", "https://ndownloader.figshare.com/files/15142850"] to "C:\Users\sebastian\.julia\datadeps\DS19 presentation materials test"?
[y/n]
y
{105620, 209345} Int64 storage graph
```

3. Set the `Plots.jl` backend to PGFPlots. This is just a technical detail required
because some of the plots are written as .tex files.

```julia
julia> pgfplots()
Plots.PGFPlotsBackend()
```
4. Decide what slide to reproduce.
You can now call the functions that were used for the slides.
To view available functions you can use tab completion. Type `slide` and then
press `tab` twice and you will see all the available options. After deciding what
slide you want to reproduce you can use the help mode to access the documentation
for that slide. In order to enter the help mode press `?` at the REPL. After that,
type the name of the function you want to know more about. For example:

```julia

help?> slide7
search: slide7 slide17_18 slide8 slide6 slide5 slide4

  No documentation found.

  DS19Presentation.slide7 is a Function.

  # 1 method for generic function "slide7":  [1] slide7(g; saveimage, savevideo) in DS19Presentation at C:\Users\sebastian\Documents\Physics\Projects\DS19Presentation\src\DS19Presentation.jl:132
```

5. Call the appropriate function.
```julia
julia> slide8(g)

```
