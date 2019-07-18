all:
	cd ../ && R CMD build ergmito/ && R CMD check --use-valgrind ergmito_*.tar.gz
