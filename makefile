check: ergmito.tar.gz
	R CMD check ergmito.tar.gz

checkv: ergmito.tar.gz
	R CMD check --use-valgrind ergmito.tar.gz

ergmito.tar.gz: R/*.R
	rm ergmito.tar.gz;\
		R CMD build . && \
		mv ergmito*.tar.gz ergmito.tar.gz
clean: ergmito.tar.gz
	rm ergmito.tar.gz

.PHONY: check checkv build clean 
