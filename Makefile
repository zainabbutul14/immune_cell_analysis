PYTHON=python3

.PHONY: setup pipeline dashboard clean

setup:
	$(PYTHON) -m pip install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) analysis.py

dashboard:
	$(PYTHON) -m streamlit run app.py

clean:
	rm -f cell_counts.db
	rm -rf outputs
