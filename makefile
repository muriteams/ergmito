all:
	cd ../ && R CMD build ergmito/ && R CMD check --as-cran --use-valgrind ergmito_*.tar.gz
