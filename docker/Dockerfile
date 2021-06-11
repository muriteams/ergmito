FROM rocker/r-base:latest

RUN apt-get update && \
	apt-get install -y --no-install-recommends libssl-dev
RUN apt-get install -y --no-install-recommends \
	r-cran-mass r-cran-rmarkdown r-cran-knitr \
	r-cran-lmtest && \
	Rscript -e 'utils::install.packages(c("fmcmc","ergm", "texreg"))'

CMD ["bash"]

