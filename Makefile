install:
	pip install -r requirements.txt

install-dev: install
	pip install -e .[dev]

lint:
	flake8 pymerk tests --max-line-length=120 --ignore=N802,W503

test:
	pytest tests