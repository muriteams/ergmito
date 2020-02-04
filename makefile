build:
	$(MAKE) ergmito.tar.gz

check: ergmito.tar.gz
	R CMD check ergmito.tar.gz

checkv: ergmito.tar.gz
	R CMD check --use-valgrind ergmito.tar.gz

ergmito.tar.gz: R/*.R
	rm ergmito.tar.gz;\
		R CMD build . && \
		mv ergmito*.tar.gz ergmito.tar.gz
clean: 
	rm ergmito.tar.gz; rm -rf ergmito.Rcheck

.PHONY: check checkv build clean 
