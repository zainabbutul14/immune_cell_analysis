# Immune Cell Analysis
Bioinformatics pipeline for immune cell population analysis using Python, SQLite and Streamlit.

# Project Overview
This project was developed as part of a technical assessment to analyze immune-cell population data from a clinical trial. The objective was to build a complete and reproducible bioinformatics workflow that transforms raw immune-cell count data into meaningful biological insights through database design, statistical analysis, and an interactive Streamlit dashboard.

Instead of treating the dataset as one large spreadsheet, I organized the data into a relational SQLite database and developed a modular workflow using Python. This approach keeps the analysis reproducible, easier to maintain, and scalable for future studies.

---

# Workflow
```
cell-count.csv
       │
       ▼
Step 1 Database Design & Data Loading
       │
       ▼
Step 2 Relative Frequency Analysis
       │
       ▼
Step 3 Statistical Analysis
       │
       ▼
Step 4 Baseline Subset Analysis
       │
       ▼
Interactive Streamlit Dashboard
```
The complete workflow can be reproduced using
```bash
make pipeline
```
This command executes the complete workflow from database creation to statistical analysis and automatically generates all required output files.

---

# Step 1 – Database Design & Data Loading
To organize the clinical dataset efficiently, I designed a normalized SQLite database consisting of five related tables:

- Projects
- Subjects
- Samples
- Cell Populations
- Cell Counts

The `load_data.py` script automatically:

- Creates the SQLite database
- Initializes the database schema
- Validates the input dataset
- Loads every record from `cell-count.csv`

Separating the data into related tables reduces duplicated information while maintaining relationships using primary and foreign keys.

---

# Step 2 – Relative Frequency Analysis

For every biological sample, the total immune-cell count is calculated by summing the counts of the five immune-cell populations.

The relative frequency of each population is calculated using:

**Relative Frequency (%) = (Cell Count / Total Cell Count) × 100**

The generated summary table contains:

- Sample ID
- Total Cell Count
- Cell Population
- Cell Count
- Relative Frequency (%)

This provides a normalized representation of immune-cell composition across all samples.

---

# Step 3 – Statistical Analysis

To investigate immune-cell populations associated with treatment response, the analysis focuses only on:

- Melanoma patients
- PBMC samples
- Miraclib-treated patients

The workflow:

- Filters the required subset of samples
- Compares responders and non-responders
- Creates boxplots for each immune-cell population
- Performs a two-sided Mann–Whitney U test
- Applies Benjamini–Hochberg False Discovery Rate (FDR) correction

The objective of this analysis is to investigate whether immune-cell population frequencies differ between responders and non-responders and to identify patterns that may support future predictive models.

After multiple-testing correction, no immune-cell population remained statistically significant, indicating that additional biological data may be required before identifying robust treatment-response biomarkers.

---

# Step 4 – Baseline Subset Analysis

To investigate early treatment effects, the database is queried to identify:

- Melanoma PBMC baseline samples (`time_from_treatment_start = 0`)
- Patients treated with **miraclib**

For this subset, the workflow reports:

- Number of baseline samples from each project
- Number of responder and non-responder subjects
- Number of male and female subjects

Finally, the workflow calculates the **average B-cell count for male melanoma responders at baseline (time = 0) across all sample and treatment types**, as required in the assessment.

The database schema was intentionally designed so that future treatments, including **quintazide**, can be incorporated without modifying the database structure.

These summaries provide a quick overview of the baseline cohort and help explore potential early treatment effects.

---

# Database Design

```
Projects
    │
    ▼
Subjects
    │
    ▼
Samples
    │
    ▼
Cell Counts
      ▲
      │
Cell Populations
```

### Why this design?

Clinical datasets often contain repeated information. A single patient can contribute multiple samples collected at different time points, while each sample contains measurements for several immune-cell populations.

Instead of storing this information repeatedly, I separated the dataset into five related tables connected through primary and foreign keys. This reduces redundancy, improves data consistency, and simplifies future analyses.

### Scalability

The database can easily support:

- Hundreds of research projects
- Thousands of patient samples
- Additional immune-cell populations
- New treatments (such as **quintazide**)
- Future analytical workflows including longitudinal studies and machine learning

Because each biological entity is stored separately and linked through primary and foreign keys, the same schema can easily accommodate larger studies without requiring structural changes.

---

# Project Structure

```
immune_cell_analysis/
│
├── load_data.py
├── analysis.py
├── app.py
├── schema.sql
├── Makefile
├── requirements.txt
├── README.md
├── cell-count.csv
└── outputs/
```

### load_data.py

Creates the SQLite database, initializes the schema, validates the dataset, and loads all records into the database.

### analysis.py

Performs all analytical tasks required in the assessment, including:

- Relative frequency calculations
- Statistical analysis
- Boxplot generation
- Baseline subset analysis
- Output table generation

### app.py

Builds an interactive Streamlit dashboard for exploring the generated results.

Separating the project into independent scripts keeps each stage of the workflow focused on a single responsibility, making the project easier to understand, maintain, and extend.

---

# Running the Project

Clone the repository:

```bash
git clone https://github.com/zainabbutul14/immune_cell_analysis.git
cd immune_cell_analysis
```

Install the required packages:

```bash
make setup
```

Run the complete workflow:

```bash
make pipeline
```

This launches the interactive Streamlit dashboard for exploring the generated analysis results.

Launch the interactive dashboard:

```bash
make dashboard
```

This launches the Streamlit application locally for interactive exploration of the analysis results.

The dashboard will be available locally at:

```
http://localhost:8501
```

---

# Dashboard

The Streamlit dashboard includes four interactive sections:

- Cell Population Frequencies
- Miraclib Response Analysis
- Baseline Subset Analysis
- Data Explorer

## Dashboard Preview

<img width="1468" height="819" alt="Screenshot 2026-07-16 at 9 15 41 PM" src="https://github.com/user-attachments/assets/c61187a9-fa8d-4219-81a6-99f28247a55e" />

---

# GitHub Repository
The complete source code is available at:
https://github.com/zainabbutul14/immune_cell_analysis
