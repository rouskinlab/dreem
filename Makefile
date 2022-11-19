DOCKER_IMAGE := ydmt/dreem
VERSION := $(shell git describe --always --dirty --long)

default: 
	python setup.py install

init:
	pip install -r requirements.txt

build-image:
	docker build .
		-f ./Dockerfile
		-t $(DOCKER_IMAGE):$(VERSION)

push-image:
	docker push $(DOCKER_IMAGE):$(VERSION)

upgrade-dependencies:
	pip uninstall -y dreem
	rm -f requirements.txt
	pip freeze > requirements.txt
	python setup.py install

push_to_pypi:
	rm -fr dist
	python3 -m build
	twine upload -r pypi dist/*
	rm -fr dist


