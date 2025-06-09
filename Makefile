PYTHON := python3

virtual:
	$(PYTHON) -m venv virtual

init: virtual
	pip install -r requirements.txt

prioritize: virtual
	$(PYTHON) py/prioritize.py

scalability: virtual
	$(PYTHON) py/scalability.py

generate-input: virtual
	$(PYTHON) tools/generate-scalability-input.py

test: virtual
