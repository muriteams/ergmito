REPO:=muriteams/ergmito

build:
	docker build -t $(REPO):release -f Dockerfile . && \
		docker build -t $(REPO):devel -f Dockerfile.dev .
push: build
	docker push $(REPO):release && \
		docker push $(REPO):devel
