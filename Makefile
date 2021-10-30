black: ## Black format every python and notebook file to line length 100
	find . -type f -name "*.py" | xargs black --line-length=100;
	find . -type f -name "*.ipynb" | xargs jblack --line-length=100;
	find . -type f -name "*.py" | xargs absolufy-imports;
	make clean;

clean: ## Remove pycache, .ipynb_checkpoints/, and .empty/
	find . -type d -name "__pycache__" | xargs rm -r;
	find . -type d -name ".ipynb_checkpoints" | xargs rm -r;
	find . -type d -name ".empty" | xargs rm -r;
