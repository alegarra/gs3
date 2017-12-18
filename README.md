<!-- pandoc README.md -f commonmark -t html -s -o README.html -->

# gs3

Software for genomic prediction and whole-genome data analysis, which name stands for *Genomic Selection --- Gibbs Sampling --- Gauss Seidel*.

The copyright of this program is owned by the INRA.

This program has been written by:
* Andrés Legarra (INRA, UMR 1388 GenPhySE)
* Anne Ricard (INRA - IFCE, UMR 1388 GenPhySE)
* Olivier Filangi (INRA - INRIA)

Other contributors:
* Timothée Flutre (INRA)

The license of this program is the [GPL](https://www.gnu.org/licenses/gpl.txt), version 3 or later.

This program has been partially financed by FEDER European funds through the [POCTEFA](http://www.poctefa.eu/) program, and by the ANR through the [Rules & Tools](https://forge-dga.jouy.inra.fr/projects/rules-tools) project.


## Compilation and installation

Run `./compile.sh` to get help.

Example: `./compile.sh gfortran`

This will generate a binary named `gs3` in the current directory.
Before using it on Windows, you should rename it as `gs3.exe`.

Older binaries can also be downloaded from [here](http://genoweb.toulouse.inra.fr/~alegarra/gs3_folder/).


## Documentation

Read the PDF file in `manual/`.
