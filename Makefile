.PHONY: install install-dev test test-cov lint format docs docs-serve clean help

# Default target
help:
	@echo "vise development commands"
	@echo ""
	@echo "  install      Install package in editable mode"
	@echo "  install-dev  Install with dev dependencies"
	@echo "  test         Run tests"
	@echo "  test-cov     Run tests with coverage"
	@echo "  test-api     Run API tests only"
	@echo "  lint         Run linter (ruff)"
	@echo "  format       Format code (ruff)"
	@echo "  docs         Build documentation"
	@echo "  docs-serve   Serve documentation locally"
	@echo "  clean        Clean build artifacts"
	@echo ""

# Installation
install:
	pip install -e .

install-dev:
	pip install -e ".[dev]"

install-all:
	pip install -e ".[dev,bolztrap]"

# Testing
test:
	pytest tests/ -v --tb=short

test-cov:
	pytest tests/ -v --cov=vise --cov-report=html --cov-report=term

test-api:
	pytest tests/api/ -v --tb=short

test-fast:
	pytest tests/ -v --tb=short -x --ignore=tests/analyzer/vasp/test_make_effective_mass.py --ignore=tests/util/phonopy/

# Linting and formatting
lint:
	ruff check vise/

format:
	ruff check --fix vise/
	ruff format vise/

# Documentation
docs:
	cd docs && make html


docs-serve:
	cd docs/build/html && python -m http.server 8000


# Cleaning
clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf vise.egg-info/
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -rf htmlcov/
	rm -rf .ruff_cache/
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete

# Build
build:
	python -m build

# Type checking
typecheck:
	mypy vise/

# Run CLI help
cli-help:
	vise --help
