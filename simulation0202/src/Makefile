RS = '/usr/bin/Rscript'
WD = /home/docker/work
SD = /home/docker/src
DD = /home/docker/data
VPATH = $(SD):$(WD)

all: checkpoint-03.rda

checkpoint-03.rda: analyze_observations.R checkpoint-02.rda
	$< >$(<F).so 2>$(<F).se

checkpoint-02.rda: simulate_observation.R checkpoint-01.rda
	$< >$(<F).so 2>$(<F).se

checkpoint-01.rda: simulate_cases.R
	$< >$(<F).so 2>$(<F).se
