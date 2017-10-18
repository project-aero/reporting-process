# Directory Structure 

Created directory ```eric-figs/``` contains code and output for figures 1, 2, and 3.

```eric-figs/ReportingProcess-ericfigs.Rmd``` is the top-level document. (see below)

```eric-figs/data/``` contains source data files, in this case consisting of saved snapshots of the complete description of the random number generator state, required to reproduce simulations and figures exactly. R rscripts 
 
```eric-figs/design/``` will contain Illustrator files under Git LFS version control (see below).

```eric-figs/images/``` will contain print-ready images exported from Illustrator

```eric-figs/install/``` contains any packages not installable from CRAN or GitHub.  These figures rely on a custom SPAERO version which Aemon wants to keep separate from the public version.

```eric-figs/output/plots``` contains plots output from R

```eric-figs/output/data``` will contain any intermediate data output from R

```eric-figs/R/``` contains help R script(s)

```eric-figs/sketches/``` contains unfinished R experiments or notes, not to be included in final released materials.

Top level code is in ```eric-figs/ReportingProcess-ericfigs.Rmd```, which calls helper scripts, loads all required packages including custom SPAERO version, and outputs figures.  *Since the working directory for R code chunks inside an .Rmd is forced to the directory of the .Rmd document itself, the working directory does not need to be set before previewing the .Rmd or running chunks, so long as the directory structure relative to the .Rmd document is preserved.* 

The ```.Rmd``` document verbosely loads dependencies, documents session info, displays folder hierarchy, then shows code used to generate figure(s); figures are output as PDFs, saved to ```eric-figs/output/plots/xxx.pdf```, then read back in, converted to PNGs and re-output to ```eric-figs/output/plots/xxx.png```.  The PNG file is then displayed in the markdown.  This way we have separate print quality output files (pdf and png) for the paper, and an image of the plot included in the R markdown for completeness and reproducibility.

For figures that will undergo cleaning up in external design software (Illustrator), the design files will be stored in ```eric-figs/design```, which will be version controlled with Git Large File Storage <https://help.github.com/articles/collaboration-with-git-large-file-storage/>.  With Git LFS enabled, only pointers to Illustrator files are pushed to github, not the actual files.  Users without LFS installed will only have access to the pointers.  Output from the Illustrator files (PDF and PNG) will be stored in ```eric-figs/images``` and displayed in the top level ```.Rmd``` file next to the direct output from R.  We can test Git LFS for shared repositories after this initial test.

