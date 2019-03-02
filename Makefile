.PHONY: clean clean-test clean-pyc clean-build docs help
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

try:
	from urllib import pathname2url
except:
	from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT


define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/

init: ## install pi
	pip install pipenv
	pipenv install --dev --skip-lock

lint: ## check style with black
	black --check --py36 pistis/ tests/

test: lint ## run tests quickly with the default Python
	pipenv run pytest

test-all: lint ## run tests on every Python version with tox
	pipenv run tox

coverage: ## check code coverage quickly with the default Python
	coverage run --source pistis -m pytest
	coverage report -m
	coverage html
	$(BROWSER) htmlcov/index.html

docs: ## generate Sphinx HTML documentation, including API docs
	rm -f docs/pistis.rst
	rm -f docs/modules.rst
	sphinx-apidoc -o docs/ pistis
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.html

servedocs: docs ## compile the docs watching for changes
	watchmedo shell-command -p '*.rst' -c '$(MAKE) -C docs html' -R -D .

release: clean lint ## package and upload a release
	python setup.py sdist upload
	python setup.py bdist_wheel upload

dist: clean lint ## builds source and wheel package
	python setup.py sdist
	twine upload -r pypi dist/`ls -t dist | head -1`
	echo "Make sure you have synchronised README files with pandoc --from=markdown --to=rst --output=README.rst README.md"


install: clean ## install the package to the active Python's site-packages
	pipenv install -e .	--skip-lock
