import os
import sqlite3

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


# ---------------------------------------------------------
# File locations
# ---------------------------------------------------------

# Get the folder where this script is saved
project_folder = os.path.dirname(os.path.abspath(__file__))

database_file = os.path.join(
    project_folder,
    "cell_counts.db"
)

output_folder = os.path.join(
    project_folder,
    "outputs"
)

plot_folder = os.path.join(
    output_folder,
    "plots"
)


# ---------------------------------------------------------
# Connect to the database
# ---------------------------------------------------------

def connect_to_database():
    """
    Connect to the SQLite database.
    """

    if not os.path.exists(database_file):
        raise FileNotFoundError(
            "cell_counts.db was not found. "
            "Please run python3 load_data.py first."
        )

    connection = sqlite3.connect(database_file)

    print("Connected to the database successfully.")

    return connection


# ---------------------------------------------------------
# Part 2: Create the frequency table
# ---------------------------------------------------------

def create_frequency_table(connection):
    """
    Calculate the total cell count and relative frequency
    of each immune-cell population for every sample.
    """

    frequency_query = """
    WITH sample_totals AS (
        SELECT
            sample_id,
            SUM(cell_count) AS total_count

        FROM cell_counts

        GROUP BY sample_id
    )

    SELECT
        samples.sample_name AS sample,
        sample_totals.total_count AS total_count,
        cell_populations.population_name AS population,
        cell_counts.cell_count AS count,

        100.0 * cell_counts.cell_count
        / sample_totals.total_count AS percentage,

        projects.project_name AS project,
        subjects.subject_name AS subject,
        subjects.indication AS indication,
        subjects.gender AS gender,
        samples.treatment AS treatment,
        samples.response AS response,
        samples.sample_type AS sample_type,
        samples.time_from_treatment_start
            AS time_from_treatment_start

    FROM cell_counts

    JOIN sample_totals
        ON cell_counts.sample_id = sample_totals.sample_id

    JOIN samples
        ON cell_counts.sample_id = samples.sample_id

    JOIN cell_populations
        ON cell_counts.population_id =
           cell_populations.population_id

    JOIN subjects
        ON samples.subject_id = subjects.subject_id

    JOIN projects
        ON subjects.project_id = projects.project_id

    ORDER BY
        samples.sample_name,
        cell_populations.population_name
    """

    frequency_data = pd.read_sql_query(
        frequency_query,
        connection
    )

    frequency_output_file = os.path.join(
        output_folder,
        "cell_population_frequencies.csv"
    )

    frequency_data.to_csv(
        frequency_output_file,
        index=False
    )

    print("\nPart 2 completed.")
    print("Frequency table created successfully.")
    print("Number of rows:", len(frequency_data))
    print("Saved as:", frequency_output_file)

    print("\nFirst ten rows:")
    print(
        frequency_data[
            [
                "sample",
                "total_count",
                "population",
                "count",
                "percentage"
            ]
        ].head(10)
    )

    return frequency_data


# ---------------------------------------------------------
# create dataframe
# ---------------------------------------------------------

def filter_miraclib_response_data(frequency_data):
    """
    Select melanoma PBMC samples treated with miraclib.

    Only responder and non-responder samples are included.
    """

    analysis_data = frequency_data[
        frequency_data["indication"]
        .str.lower()
        .eq("melanoma")
        &
        frequency_data["treatment"]
        .str.lower()
        .eq("miraclib")
        &
        frequency_data["sample_type"]
        .str.lower()
        .eq("pbmc")
        &
        frequency_data["response"]
        .isin(["yes", "no"])
    ].copy()

    print(
        "\nNumber of rows used for the "
        "miraclib response analysis:",
        len(analysis_data)
    )

    return analysis_data


# ---------------------------------------------------------
# Part 3: Statistical analysis
# ---------------------------------------------------------

def perform_statistical_analysis(frequency_data):
    """
    Compare relative cell frequencies between responders
    and non-responders.

    A Mann-Whitney U test is used for each cell population.
    Benjamini-Hochberg correction is then applied because
    five populations are tested.
    """

    analysis_data = filter_miraclib_response_data(
        frequency_data
    )

    result_list = []

    populations = sorted(
        analysis_data["population"].unique()
    )

    for population in populations:

        population_data = analysis_data[
            analysis_data["population"] == population
        ]

        responder_values = population_data[
            population_data["response"] == "yes"
        ]["percentage"]

        nonresponder_values = population_data[
            population_data["response"] == "no"
        ]["percentage"]

        test_result = mannwhitneyu(
            responder_values,
            nonresponder_values,
            alternative="two-sided"
        )

        responder_mean = responder_values.mean()
        nonresponder_mean = nonresponder_values.mean()

        result_list.append(
            {
                "population": population,
                "responder_samples":
                    len(responder_values),
                "nonresponder_samples":
                    len(nonresponder_values),
                "responder_mean_percentage":
                    responder_mean,
                "nonresponder_mean_percentage":
                    nonresponder_mean,
                "difference_in_means":
                    responder_mean - nonresponder_mean,
                "mann_whitney_u":
                    test_result.statistic,
                "p_value":
                    test_result.pvalue
            }
        )

    statistical_results = pd.DataFrame(
        result_list
    )

    correction_results = multipletests(
        statistical_results["p_value"],
        alpha=0.05,
        method="fdr_bh"
    )

    statistical_results[
        "adjusted_p_value"
    ] = correction_results[1]

    statistical_results[
        "significant"
    ] = correction_results[0]

    statistical_results = statistical_results.sort_values(
        "adjusted_p_value"
    )

    statistics_output_file = os.path.join(
        output_folder,
        "miraclib_response_statistics.csv"
    )

    statistical_results.to_csv(
        statistics_output_file,
        index=False
    )

    print("\nPart 3 statistical analysis completed.")
    print("Saved as:", statistics_output_file)

    print("\nStatistical results:")
    print(
        statistical_results[
            [
                "population",
                "responder_mean_percentage",
                "nonresponder_mean_percentage",
                "p_value",
                "adjusted_p_value",
                "significant"
            ]
        ]
    )

    return statistical_results


# ---------------------------------------------------------
# Part 3: Create the boxplot
# ---------------------------------------------------------

def create_response_boxplot(frequency_data):
    """
    Create a boxplot comparing responders and
    non-responders for each immune-cell population.
    """

    plot_data = filter_miraclib_response_data(
        frequency_data
    )

    populations = sorted(
        plot_data["population"].unique()
    )

    plot_values = []
    plot_positions = []
    plot_labels = []

    position = 1

    for population in populations:

        for response in ["yes", "no"]:

            values = plot_data[
                (
                    plot_data["population"] == population
                )
                &
                (
                    plot_data["response"] == response
                )
            ]["percentage"].to_numpy()

            plot_values.append(values)
            plot_positions.append(position)

            plot_labels.append(
                f"{population}\n{response}"
            )

            position = position + 1

        # Add extra space between populations
        position = position + 0.5

    figure, axis = plt.subplots(
        figsize=(14, 7)
    )

    axis.boxplot(
        plot_values,
        positions=plot_positions,
        widths=0.7
    )

    axis.set_xticks(
        plot_positions
    )

    axis.set_xticklabels(
        plot_labels,
        rotation=30,
        ha="right"
    )

    axis.set_xlabel(
        "Immune-cell population and response"
    )

    axis.set_ylabel(
        "Relative frequency (%)"
    )

    axis.set_title(
        "Melanoma PBMC Samples Treated with Miraclib\n"
        "Responders versus Non-Responders"
    )

    axis.grid(
        axis="y",
        alpha=0.3
    )

    figure.tight_layout()

    plot_output_file = os.path.join(
        plot_folder,
        "miraclib_response_boxplots.png"
    )

    figure.savefig(
        plot_output_file,
        dpi=200
    )

    plt.close(figure)

    print("\nBoxplot created successfully.")
    print("Saved as:", plot_output_file)


# ---------------------------------------------------------
# Part 4: Baseline subset analysis
# ---------------------------------------------------------

def perform_baseline_analysis(connection):
    """
    Identify melanoma PBMC samples at baseline
    from patients treated with miraclib.
    """

    baseline_query = """
    SELECT
        projects.project_name AS project,
        subjects.subject_name AS subject,
        samples.sample_name AS sample,
        samples.response AS response,
        subjects.gender AS gender,
        samples.sample_type AS sample_type,
        samples.treatment AS treatment,
        samples.time_from_treatment_start
            AS time_from_treatment_start

    FROM samples

    JOIN subjects
        ON samples.subject_id = subjects.subject_id

    JOIN projects
        ON subjects.project_id = projects.project_id

    WHERE LOWER(subjects.indication) = 'melanoma'
      AND LOWER(samples.sample_type) = 'pbmc'
      AND LOWER(samples.treatment) = 'miraclib'
      AND samples.time_from_treatment_start = 0

    ORDER BY
        projects.project_name,
        subjects.subject_name,
        samples.sample_name
    """

    baseline_data = pd.read_sql_query(
        baseline_query,
        connection
    )

    baseline_output_file = os.path.join(
        output_folder,
        "baseline_melanoma_pbmc_miraclib.csv"
    )

    baseline_data.to_csv(
        baseline_output_file,
        index=False
    )

    print("\nPart 4 baseline subset created.")
    print("Number of baseline samples:", len(baseline_data))
    print("Saved as:", baseline_output_file)

    # -----------------------------------------------------
    # Number of samples from each project
    # -----------------------------------------------------

    samples_by_project = (
        baseline_data
        .groupby("project")
        .size()
        .reset_index(name="sample_count")
    )

    samples_by_project_file = os.path.join(
        output_folder,
        "baseline_samples_by_project.csv"
    )

    samples_by_project.to_csv(
        samples_by_project_file,
        index=False
    )

    print("\nSamples from each project:")
    print(samples_by_project)

    # -----------------------------------------------------
    # Create one row per unique subject
    # -----------------------------------------------------

    unique_subjects = baseline_data.drop_duplicates(
        [
            "project",
            "subject"
        ]
    )

    # -----------------------------------------------------
    # Number of subjects by response
    # -----------------------------------------------------

    subjects_by_response = (
        unique_subjects
        .groupby("response")
        .size()
        .reset_index(name="subject_count")
    )

    subjects_by_response_file = os.path.join(
        output_folder,
        "baseline_subjects_by_response.csv"
    )

    subjects_by_response.to_csv(
        subjects_by_response_file,
        index=False
    )

    print("\nSubjects by response:")
    print(subjects_by_response)

    # -----------------------------------------------------
    # Number of subjects by gender
    # -----------------------------------------------------

    subjects_by_gender = (
        unique_subjects
        .groupby("gender")
        .size()
        .reset_index(name="subject_count")
    )

    subjects_by_gender_file = os.path.join(
        output_folder,
        "baseline_subjects_by_gender.csv"
    )

    subjects_by_gender.to_csv(
        subjects_by_gender_file,
        index=False
    )

    print("\nSubjects by gender:")
    print(subjects_by_gender)

    return baseline_data


# ---------------------------------------------------------
# Part 4: Average B-cell count
# ---------------------------------------------------------

def calculate_average_b_cells(connection):
    """
    Calculate the average B-cell count for male melanoma
    responders at time 0.

    All sample types and treatment types are included,
    as requested in the assignment.
    """

    b_cell_query = """
    SELECT
        AVG(cell_counts.cell_count)
            AS average_b_cell_count

    FROM cell_counts

    JOIN cell_populations
        ON cell_counts.population_id =
           cell_populations.population_id

    JOIN samples
        ON cell_counts.sample_id = samples.sample_id

    JOIN subjects
        ON samples.subject_id = subjects.subject_id

    WHERE LOWER(cell_populations.population_name) = 'b_cell'
      AND LOWER(subjects.indication) = 'melanoma'
      AND UPPER(subjects.gender) = 'M'
      AND LOWER(samples.response) = 'yes'
      AND samples.time_from_treatment_start = 0
    """

    result = connection.execute(
        b_cell_query
    ).fetchone()

    average_b_cells = result[0]

    if average_b_cells is None:
        raise ValueError(
            "No samples matched the B-cell query."
        )

    average_b_cells = round(
        float(average_b_cells),
        2
    )

    answer_data = pd.DataFrame(
        {
            "average_b_cell_count": [
                average_b_cells
            ]
        }
    )

    answer_output_file = os.path.join(
        output_folder,
        "male_melanoma_responder_b_cell_average.csv"
    )

    answer_data.to_csv(
        answer_output_file,
        index=False
    )

    print(
        "\nAverage B-cell count for male melanoma "
        "responders at time 0:"
    )

    print(
        f"{average_b_cells:.2f}"
    )

    print(
        "Saved as:",
        answer_output_file
    )

    return average_b_cells


# ---------------------------------------------------------
# Show significant populations
# ---------------------------------------------------------

def show_significant_results(statistical_results):
    """
    Display cell populations with a significant difference
    after multiple-testing correction.
    """

    significant_results = statistical_results[
        statistical_results["significant"] == True
    ]

    if len(significant_results) == 0:

        print(
            "\nNo immune-cell populations showed a "
            "significant difference after "
            "multiple-testing correction."
        )

    else:

        print(
            "\nSignificant immune-cell populations:"
        )

        print(
            significant_results[
                [
                    "population",
                    "adjusted_p_value"
                ]
            ]
        )


# ---------------------------------------------------------
# Run the complete analysis pipeline
# ---------------------------------------------------------

def main():
    """
    Run Parts 2, 3, and 4 of the analysis.
    """

    # Create output folders if they do not already exist
    os.makedirs(
        output_folder,
        exist_ok=True
    )

    os.makedirs(
        plot_folder,
        exist_ok=True
    )

    connection = connect_to_database()

    try:

        frequency_data = create_frequency_table(
            connection
        )

        statistical_results = perform_statistical_analysis(
            frequency_data
        )

        create_response_boxplot(
            frequency_data
        )

        perform_baseline_analysis(
            connection
        )

        calculate_average_b_cells(
            connection
        )

        show_significant_results(
            statistical_results
        )

        print(
            "\nComplete analysis pipeline "
            "finished successfully."
        )

    finally:

        connection.close()

        print(
            "Database connection closed."
        )


if __name__ == "__main__":
    main()