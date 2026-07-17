import os
import sqlite3

import pandas as pd


# ---------------------------------------------------------
# File locations
# ---------------------------------------------------------

# Get the folder where this Python script is saved
project_folder = os.path.dirname(os.path.abspath(__file__))

csv_file = os.path.join(project_folder, "cell-count.csv")
database_file = os.path.join(project_folder, "cell_counts.db")
schema_file = os.path.join(project_folder, "schema.sql")


# These are the five immune-cell populations in the dataset
cell_populations = [
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte"
]


# ---------------------------------------------------------
# Check required files
# ---------------------------------------------------------

def check_files():
    """
    Check whether the CSV file and schema file are available.
    """

    if not os.path.exists(csv_file):
        raise FileNotFoundError(
            "cell-count.csv was not found in the project folder."
        )

    if not os.path.exists(schema_file):
        raise FileNotFoundError(
            "schema.sql was not found in the project folder."
        )

    print("All required files were found.")


# ---------------------------------------------------------
# Read the CSV file
# ---------------------------------------------------------

def read_data():
    """
    Read cell-count.csv using pandas.
    """

    data = pd.read_csv(csv_file)

    print("CSV file was read successfully.")
    print("Number of rows:", len(data))
    print("Number of columns:", len(data.columns))

    print("\nColumn names:")
    print(data.columns.tolist())

    return data


# ---------------------------------------------------------
# Validate the CSV data
# ---------------------------------------------------------

def validate_data(data):
    """
    Check that the CSV contains the required columns
    and that the values are valid.
    """

    required_columns = [
        "project",
        "subject",
        "condition",
        "age",
        "sex",
        "treatment",
        "response",
        "sample",
        "sample_type",
        "time_from_treatment_start",
        "b_cell",
        "cd8_t_cell",
        "cd4_t_cell",
        "nk_cell",
        "monocyte"
    ]

    missing_columns = []

    for column in required_columns:
        if column not in data.columns:
            missing_columns.append(column)

    if missing_columns:
        raise ValueError(
            f"These columns are missing from the CSV: {missing_columns}"
        )

    # Sample IDs should be unique
    if data["sample"].duplicated().any():
        raise ValueError(
            "Duplicate sample IDs were found in the sample column."
        )

    # Cell counts should not be negative
    if (data[cell_populations] < 0).any().any():
        raise ValueError(
            "Negative cell-count values were found."
        )

    print("Data validation completed successfully.")


# ---------------------------------------------------------
# Create the SQLite database
# ---------------------------------------------------------

def create_database():
    """
    Create a new SQLite database.

    If an old database already exists, remove it so that
    every run starts with a clean database.
    """

    if os.path.exists(database_file):
        os.remove(database_file)

    connection = sqlite3.connect(database_file)

    # Make sure SQLite uses foreign-key relationships
    connection.execute("PRAGMA foreign_keys = ON")

    print("SQLite database created successfully.")

    return connection


# ---------------------------------------------------------
# Create database tables
# ---------------------------------------------------------

def create_tables(connection):
    """
    Read schema.sql and create the database tables.
    """

    with open(schema_file, "r") as file:
        schema = file.read()

    connection.executescript(schema)

    print("Database tables created successfully.")


# ---------------------------------------------------------
# Check created tables
# ---------------------------------------------------------

def check_tables(connection):
    """
    Display the tables created in the SQLite database.
    """

    query = """
    SELECT name
    FROM sqlite_master
    WHERE type = 'table'
    ORDER BY name
    """

    tables = connection.execute(query).fetchall()

    print("\nTables in the database:")

    for table in tables:
        print(table[0])


# ---------------------------------------------------------
# Load projects
# ---------------------------------------------------------

def load_projects(connection, data):
    """
    Load unique project names into the projects table
    and create a project lookup dictionary.
    """

    unique_projects = data["project"].drop_duplicates()

    project_records = []

    for project in unique_projects:
        project_records.append(
            (project,)
        )

    connection.executemany(
        """
        INSERT INTO projects (project_name)
        VALUES (?)
        """,
        project_records
    )

    connection.commit()

    project_rows = connection.execute(
        """
        SELECT project_id, project_name
        FROM projects
        """
    ).fetchall()

    project_lookup = {}

    for project_id, project_name in project_rows:
        project_lookup[project_name] = project_id

    print("Projects loaded successfully.")
    print("Number of projects loaded:", len(project_lookup))

    return project_lookup


# ---------------------------------------------------------
# Load subjects
# ---------------------------------------------------------

def load_subjects(connection, data, project_lookup):
    """
    Load unique subjects into the subjects table
    and create a subject lookup dictionary.
    """

    subject_columns = [
        "project",
        "subject",
        "condition",
        "age",
        "sex"
    ]

    unique_subjects = data[
        subject_columns
    ].drop_duplicates()

    subject_records = []

    for _, row in unique_subjects.iterrows():

        project_id = project_lookup[
            row["project"]
        ]

        subject_records.append(
            (
                project_id,
                row["subject"],
                row["condition"],
                int(row["age"]),
                row["sex"]
            )
        )

    connection.executemany(
        """
        INSERT INTO subjects (
            project_id,
            subject_name,
            indication,
            age,
            gender
        )
        VALUES (?, ?, ?, ?, ?)
        """,
        subject_records
    )

    connection.commit()

    subject_rows = connection.execute(
        """
        SELECT
            projects.project_name,
            subjects.subject_name,
            subjects.subject_id

        FROM subjects

        JOIN projects
            ON subjects.project_id = projects.project_id
        """
    ).fetchall()

    subject_lookup = {}

    for project_name, subject_name, subject_id in subject_rows:

        subject_lookup[
            (project_name, subject_name)
        ] = subject_id

    print("Subjects loaded successfully.")
    print("Number of subjects loaded:", len(subject_lookup))

    return subject_lookup


# ---------------------------------------------------------
# Load samples
# ---------------------------------------------------------

def load_samples(connection, data, subject_lookup):
    """
    Load sample-level information into the samples table
    and create a sample lookup dictionary.
    """

    sample_records = []

    for _, row in data.iterrows():

        subject_id = subject_lookup[
            (row["project"], row["subject"])
        ]

        response = row["response"]

        # Store missing response values as NULL in SQLite
        if pd.isna(response):
            response = None

        sample_records.append(
            (
                row["sample"],
                subject_id,
                row["treatment"],
                response,
                row["sample_type"],
                int(row["time_from_treatment_start"])
            )
        )

    connection.executemany(
        """
        INSERT INTO samples (
            sample_name,
            subject_id,
            treatment,
            response,
            sample_type,
            time_from_treatment_start
        )
        VALUES (?, ?, ?, ?, ?, ?)
        """,
        sample_records
    )

    connection.commit()

    sample_rows = connection.execute(
        """
        SELECT sample_name, sample_id
        FROM samples
        """
    ).fetchall()

    sample_lookup = {}

    for sample_name, sample_id in sample_rows:
        sample_lookup[sample_name] = sample_id

    print("Samples loaded successfully.")
    print("Number of samples loaded:", len(sample_lookup))

    return sample_lookup


# ---------------------------------------------------------
# Load immune-cell population names
# ---------------------------------------------------------

def load_cell_populations(connection):
    """
    Load the five immune-cell population names
    and create a population lookup dictionary.
    """

    population_records = []

    for population in cell_populations:
        population_records.append(
            (population,)
        )

    connection.executemany(
        """
        INSERT INTO cell_populations (population_name)
        VALUES (?)
        """,
        population_records
    )

    connection.commit()

    population_rows = connection.execute(
        """
        SELECT population_name, population_id
        FROM cell_populations
        """
    ).fetchall()

    population_lookup = {}

    for population_name, population_id in population_rows:
        population_lookup[population_name] = population_id

    print("Cell populations loaded successfully.")
    print(
        "Number of cell populations loaded:",
        len(population_lookup)
    )

    return population_lookup


# ---------------------------------------------------------
# Load cell-count values
# ---------------------------------------------------------

def load_cell_counts(
    connection,
    data,
    sample_lookup,
    population_lookup
):
    """
    Load one cell-count record for each population
    in each biological sample.
    """

    cell_count_records = []

    for _, row in data.iterrows():

        sample_id = sample_lookup[
            row["sample"]
        ]

        for population in cell_populations:

            population_id = population_lookup[
                population
            ]

            count = int(
                row[population]
            )

            cell_count_records.append(
                (
                    sample_id,
                    population_id,
                    count
                )
            )

    connection.executemany(
        """
        INSERT INTO cell_counts (
            sample_id,
            population_id,
            cell_count
        )
        VALUES (?, ?, ?)
        """,
        cell_count_records
    )

    connection.commit()

    number_of_records = connection.execute(
        """
        SELECT COUNT(*)
        FROM cell_counts
        """
    ).fetchone()[0]

    print("Cell-count values loaded successfully.")
    print(
        "Number of cell-count records loaded:",
        number_of_records
    )


# ---------------------------------------------------------
# Final database check
# ---------------------------------------------------------

def show_database_summary(connection):
    """
    Display a summary of the number of rows
    stored in each main database table.
    """

    table_names = [
        "projects",
        "subjects",
        "samples",
        "cell_populations",
        "cell_counts"
    ]

    print("\nDatabase summary:")

    for table_name in table_names:

        query = f"""
        SELECT COUNT(*)
        FROM {table_name}
        """

        count = connection.execute(
            query
        ).fetchone()[0]

        print(f"{table_name}: {count} rows")


# ---------------------------------------------------------
# Run the complete loading process
# ---------------------------------------------------------

def main():
    """
    Run the complete database-loading pipeline.
    """

    check_files()

    data = read_data()

    validate_data(data)

    connection = create_database()

    try:
        create_tables(connection)

        check_tables(connection)

        project_lookup = load_projects(
            connection,
            data
        )

        subject_lookup = load_subjects(
            connection,
            data,
            project_lookup
        )

        sample_lookup = load_samples(
            connection,
            data,
            subject_lookup
        )

        population_lookup = load_cell_populations(
            connection
        )

        load_cell_counts(
            connection,
            data,
            sample_lookup,
            population_lookup
        )

        show_database_summary(connection)

        print(
            "\nDatabase loading completed successfully."
        )

        print(
            "Database created:",
            database_file
        )

    finally:
        connection.close()

        print(
            "Database connection closed."
        )


if __name__ == "__main__":
    main()
