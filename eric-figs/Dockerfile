FROM r-base:3.4.3
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN apt-get update && apt-get install -y -q --no-install-recommends \
libcurl4-openssl-dev \
libgdal-dev \
libnlopt-dev \
libproj-dev \
libpng-dev \
libpoppler-cpp-dev \
librsvg2-dev \
libssl-dev \
libwebp-dev \
libxml2-dev

RUN install2.r --error packrat
RUN mkdir /root/work
COPY ./packrat /root/work/packrat
COPY ./.Rprofile /root/work/
RUN /usr/bin/Rscript -e "setwd(\"/root/work\"); source(\".Rprofile\"); packrat::restore()"
