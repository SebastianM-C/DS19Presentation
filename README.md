# DS19Presentation

![conference_logo](assets/DS19_logo.png)

This repository contains Materials for the SIAM DS19 Presentation
entitled "Large-Scale Numerical Investigations into the Dynamics of Nonlinear Classical Systems".

In order to ensure the reproducibility of the results all the
computations were performed inside a [singularity container](https://www.sylabs.io/singularity/). The benefits of this approach are
twofold, as having containers simplifies the deployment in HPC scenarios. Singularity containers are created from definition files.
The definition file used in this case is `julia.def`. After building
the container(`sudo singularity build name.sif julia.def`) and opening a shell inside (`singularity shell name.sif`) one can start
reproducing the Julia environment, as presented bellow.

# Notes

- It is not necessary to use singularity to reproduce the results,
but by using it you can be sure that the environment is the same.
- The dataset required to reproduce the results has about 1GB, so
depending on the Internet connection it may take a longer time to
download.
- If you encounter errors or problems when trying to reproduce the results, please open an issue with the detailed steps that you made
and the errors or problems encountered.

# Julia instructions

In order to reproduce the results you need to reproduce the environment in which
they were created. For this task, [Julia projects](https://julialang.github.io/Pkg.jl/v1/) use `Manifest.toml` files. Thus,
the package manager will install the exact same package versions as the ones
that were used by the author.
Moreover, the main dataset was uploaded to [figshare](https://figshare.com/articles/DS19_presentation_materials/8146106) and with the help of the
[DataDeps package](https://github.com/oxinabox/DataDeps.jl) it will be automatically downloaded on the first use.

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
```
After the instalations finish, exit `Pkg` mode by pressing `backspace` (you should now have the `julia>` prompt again).
```julia
julia> using DS19Presentation
[ Info: Precompiling DS19Presentation [c52eaf16-6d25-11e9-34da-d5c7cd6dfd33]
WARNING: could not import Base.quit into AtomShell
```

3. Load the main dataset. The first time you will do this, you will be asked if
you want to download it. To exit the `Pkg` mode press `backspace`.

```julia
julia> g=load()
This program has requested access to the data dependency DS19 presentation materials.
which is not currently installed. It can be installed automatically, and you will not see this message again.

Dataset: DS19 presentation materials
Website: https://figshare.com/articles/DS19_presentation_materials/8146106
Author: Sebastian Micluța-Câmpeanu
Date of Publication: 2019-05-17T22:19:05Z
License: MIT (https://opensource.org/licenses/MIT)

Here are the materials for the SIAM DS19 presentation entitled "Large-Scale Numerical Investigations into the Dynamics of Nonlinear Classical Systems". The dataset used for all the computations is included as the binary file graph.jls and it requires julia 1.1.0 and the two julia packages linked below to be correctly read. See the DS19Presentation GitHub repository linked below for more details on how to reproduce the presented figures.

Please cite this dataset: Micluța-Câmpeanu, Sebastian (2019): DS19 presentation materials. figshare. Conference contribution.



Do you want to download the dataset from Any["https://ndownloader.figshare.com/files/15181136"] to "/home/user/.julia/datadeps/DS19 presentation materials"?
[y/n]
y
{105627, 209359} Int64 storage graph
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
help?> slide5
search: slide5 slide15 slide_a5_6 slide24_25 slide6 slide4 slide28 slide27 slide13 slide_a8 slide_a7 slide_a2 slide_a1 slide8_9 slide_a3_4 slide30_31 slide22_23

  slide5(g; saveimage=false, savevideo=false)

  This function returns the time node and the scene for the simulation of the nuclear surface and its sections. Call as t, sc = slide5(g); and then display the
  scene using display(sc) and use animate(t, (0, 40)) to animate it.
```

5. Call the appropriate function.
```julia
julia> slide15(g)

```
