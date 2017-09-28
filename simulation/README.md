This subdirectory contains files used to produce simulation results as
well as the simulation results that were used to generate the figures
in the paper. The latter is contained in the zip archive
``work-2017-08-29_17-29-29-light.zip``. Scripts used to generate the
results as well as some other non-essential but related scripts are in
``src``.

To enhance reproduciblity, the scripts were run within a
Docker container based on an image that is available on Docker
Hub. The file ``build-image`` was used to create the image based on
the ``Dockerfile``. On systems with Docker installed, when the bash
script ``run-scripts`` is run, it should pull the Docker container
down from Docker Hub and run the Makefile in ``src`` within the
container to reproduce the full simulation results. Note that the
script is meant to be run from the parent directory to ``src`` and it
will create a time-stamped child directory containing the simulation
output. This output includes many intermediate results that were
removed from the previously-mentioned zip archive to make its size
suitable for distribution. The script ``src/make-light.sh`` was used
to generate this lighter version of the results.

It is also possible to run the scripts without Docker, although one
might have to consult the Docker image to determine the correct
versions of the software to use. One important dependency, the spaero
R package for calculating the moving window statistics, is included in
this subdirectory. The right sequence in which the scripts must be run
to produce the results can be seen in ``src/Makefile``.
