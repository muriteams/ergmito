check:
	$(MAKE) build && R CMD check ../ergmito_*.tar.gz && $(MAKE) clean

checkv:
	$(MAKE) build && R CMD check --use-valgrind ../ergmito_*.tar.gz && $(MAKE) clean

build:
	cd ../ && R CMD build ergmito/ && cd ergmito

clean:
	rm  ../ergmito_*.tar.gz

install:
	Rscript -e "Rcpp::compileAttributes();devtools::document()" && R CMD INSTALL --preclean .

.PHONY: check checkv build clean 
