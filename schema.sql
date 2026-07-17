PRAGMA foreign_keys = ON;

CREATE TABLE projects (
    project_id INTEGER PRIMARY KEY AUTOINCREMENT,
    project_name TEXT UNIQUE NOT NULL
);

CREATE TABLE subjects (
    subject_id INTEGER PRIMARY KEY AUTOINCREMENT,
    project_id INTEGER NOT NULL,
    subject_name TEXT NOT NULL,
    indication TEXT,
    age INTEGER,
    gender TEXT,

    UNIQUE(project_id, subject_name),

    FOREIGN KEY(project_id)
        REFERENCES projects(project_id)
);

CREATE TABLE samples (
    sample_id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_name TEXT UNIQUE NOT NULL,
    subject_id INTEGER NOT NULL,
    treatment TEXT,
    response TEXT,
    sample_type TEXT,
    time_from_treatment_start INTEGER,

    FOREIGN KEY(subject_id)
        REFERENCES subjects(subject_id)
);

CREATE TABLE cell_populations (
    population_id INTEGER PRIMARY KEY AUTOINCREMENT,
    population_name TEXT UNIQUE NOT NULL
);

CREATE TABLE cell_counts (
    sample_id INTEGER NOT NULL,
    population_id INTEGER NOT NULL,
    cell_count INTEGER NOT NULL,

    PRIMARY KEY(sample_id, population_id),

    FOREIGN KEY(sample_id)
        REFERENCES samples(sample_id),

    FOREIGN KEY(population_id)
        REFERENCES cell_populations(population_id)
);