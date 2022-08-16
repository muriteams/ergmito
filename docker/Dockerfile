FROM rocker/r-base:latest

RUN apt-get update && \
	apt-get install -y --no-install-recommends libssl-dev libcurl4-openssl-dev
RUN apt-get install -y --no-install-recommends \
	r-cran-mass r-cran-rmarkdown r-cran-knitr \
	r-cran-lmtest r-cran-rcpparmadillo && \
	Rscript -e 'utils::install.packages(c("fmcmc","ergm", "texreg"))'

RUN Rscript -e 'utils::install.packages(c("sna", "covr","tinytest"))'

CMD ["bash"]

