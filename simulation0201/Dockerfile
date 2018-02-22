FROM r-base
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN apt-get update && apt-get install -y -q --no-install-recommends \
  openssh-server \
  r-cran-coda \
  r-cran-curl \
  r-cran-desolve \
  r-cran-doparallel \
  r-cran-evaluate \
  r-cran-foreach \
  r-cran-formatr \
  r-cran-ggplot2 \
  r-cran-httr \
  r-cran-igraph \
  r-cran-iterators \
  r-cran-jsonlite \
  r-cran-lazyeval \
  r-cran-memoise \
  r-cran-mime \
  r-cran-mvtnorm \
  r-cran-nloptr \
  r-cran-openssl \
  r-cran-praise \
  r-cran-rcpp \
  r-cran-r6 \
  r-cran-testthat \
  r-cran-tgp \
  r-cran-xtable \
  r-cran-yaml

RUN install2.r --error \
  earlywarnings \
  knitr \
  lintr \
  np \
  pomp \
  rmarkdown

COPY ./spaero_0.2.0.1.tar.gz /home/docker/spaero_0.2.0.1.tar.gz
RUN /usr/bin/R CMD check --no-build-vignettes --no-manual /home/docker/spaero_*.tar.gz > /home/docker/check.so \
&& /usr/bin/R CMD INSTALL --build /home/docker/spaero_*.tar.gz

RUN mkdir /var/run/sshd && echo 'docker:screencast' | chpasswd \
  && sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd
EXPOSE 22

CMD ["/usr/sbin/sshd", "-D"]
