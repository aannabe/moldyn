init:
	pip install -r requirements.txt

test:
	pytest

# Target for cleaning up temporary files
clean:
	find . -name "*.pyc" -exec rm -f {} +
	find . -name "__pycache__" -exec rm -rf {} +
	rm -rf build dist *.egg-info

.PHONY: init test clean
