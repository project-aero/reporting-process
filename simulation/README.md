The ``simulation`` subdirectory contains:
1) simulation results used to generate the figures in the paper contained in the zip archive.
2) All scripts used to generate the simulation results in the directory ``src``.
3) Docker files (see below)
4) A copy of the R package spaero
5) Some additional non-essential but related scripts in ``src``.

To facilitate reproducibility, the scripts were run within a Docker container based on an image that is available on Docker Hub. The file ``build-image`` was used to create the image based on
the ``Dockerfile``. To reproduce the results on systems with Docker installed, run the bash script ``run-scripts`` the ``simulation`` subdirectory. This script will:

1) Pull the Docker container down from Docker Hub and create a time-stamped child directory 
2) Run the Makefile in ``src`` to reproduce the full simulation results in the time-stamped child directory

Reproducing the full output is computationally intensive and not recommended on a desktop computer. The full output comes to about 30GB and includes many intermediate results which are not needed to recreate the figures. The script ``src/make-light.sh`` is included to generate a lighter compressed zip version of the results and only includes the data necessary for figures.

It is also possible to run the scripts without Docker, although one might have to consult the Docker image to determine the correct software versions to use. One important dependency, the spaero
R package for calculating the moving window statistics, is included in this subdirectory. The right sequence in which the scripts must be run to produce the results can be seen in ``src/Makefile``.
