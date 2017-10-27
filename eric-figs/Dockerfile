FROM r-base
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN apt-get update && apt-get install -y -q --no-install-recommends \
libcurl4-openssl-dev \
libxml2-dev

RUN install2.r --error packrat
RUN mkdir /root/work
COPY ./packrat /root/work/packrat
COPY ./.Rprofile /root/work
RUN /usr/bin/Rscript -e "setwd(\"/root/work\"); source(\".Rprofile\"); packrat::restore()"
