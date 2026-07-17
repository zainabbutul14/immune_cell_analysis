import os

import pandas as pd
import plotly.express as px
import streamlit as st


# ---------------------------------------------------------
# File locations
# ---------------------------------------------------------

# Get the folder where this script is saved
project_folder = os.path.dirname(
    os.path.abspath(__file__)
)

output_folder = os.path.join(
    project_folder,
    "outputs"
)

frequency_file = os.path.join(
    output_folder,
    "cell_population_frequencies.csv"
)

statistics_file = os.path.join(
    output_folder,
    "miraclib_response_statistics.csv"
)

baseline_file = os.path.join(
    output_folder,
    "baseline_melanoma_pbmc_miraclib.csv"
)

project_counts_file = os.path.join(
    output_folder,
    "baseline_samples_by_project.csv"
)

response_counts_file = os.path.join(
    output_folder,
    "baseline_subjects_by_response.csv"
)

gender_counts_file = os.path.join(
    output_folder,
    "baseline_subjects_by_gender.csv"
)

b_cell_average_file = os.path.join(
    output_folder,
    "male_melanoma_responder_b_cell_average.csv"
)


# ---------------------------------------------------------
# Dashboard page setup
# ---------------------------------------------------------

st.set_page_config(
    page_title="Immune Cell Analysis",
    layout="wide"
)

st.title("Immune Cell Population Analysis")

st.write(
    """
        I built this dashboard to explore immune-cell composition
        across clinical samples and investigate whether relative
        cell frequencies differ between miraclib responders and
        non-responders.
        """
)


# ---------------------------------------------------------
# Check that analysis output files exist
# ---------------------------------------------------------

required_files = [
    frequency_file,
    statistics_file,
    baseline_file,
    project_counts_file,
    response_counts_file,
    gender_counts_file,
    b_cell_average_file
]

missing_files = []

for file in required_files:
    if not os.path.exists(file):
        missing_files.append(file)


if missing_files:

    st.error(
        "Some analysis output files are missing. "
        "Please run python3 load_data.py and "
        "python3 analysis.py first."
    )

    st.write("Missing files:")

    for file in missing_files:
        st.write(file)

    st.stop()


# ---------------------------------------------------------
# Load analysis output files
# ---------------------------------------------------------

@st.cache_data
def load_dashboard_data():
    """
    Load all CSV files used in the dashboard.
    """

    frequency_data = pd.read_csv(
        frequency_file
    )

    statistics_data = pd.read_csv(
        statistics_file
    )

    baseline_data = pd.read_csv(
        baseline_file
    )

    project_counts = pd.read_csv(
        project_counts_file
    )

    response_counts = pd.read_csv(
        response_counts_file
    )

    gender_counts = pd.read_csv(
        gender_counts_file
    )

    b_cell_average = pd.read_csv(
        b_cell_average_file
    )

    return (
        frequency_data,
        statistics_data,
        baseline_data,
        project_counts,
        response_counts,
        gender_counts,
        b_cell_average
    )


(
    frequency_data,
    statistics_data,
    baseline_data,
    project_counts,
    response_counts,
    gender_counts,
    b_cell_average
) = load_dashboard_data()


# ---------------------------------------------------------
# Create dashboard tabs
# ---------------------------------------------------------

tab1, tab2, tab3, tab4 = st.tabs(
    [
        "Cell Frequencies",
        "Miraclib Response",
        "Baseline Subset",
        "Data Explorer"
    ]
)


# ---------------------------------------------------------
# Tab 1: Cell population frequencies
# ---------------------------------------------------------

with tab1:

    st.header(
        "Relative Frequency of Each Cell Population"
    )

    st.write(
        """
        Select one or more samples to view the cell counts
        and relative frequencies for the five immune-cell
        populations.
        """
    )

    sample_names = sorted(
        frequency_data["sample"].unique()
    )

    default_samples = sample_names[:5]

    selected_samples = st.multiselect(
        "Select sample IDs",
        options=sample_names,
        default=default_samples
    )

    if selected_samples:

        selected_data = frequency_data[
            frequency_data["sample"].isin(
                selected_samples
            )
        ].copy()

        display_table = selected_data[
            [
                "sample",
                "total_count",
                "population",
                "count",
                "percentage"
            ]
        ].copy()

        display_table["percentage"] = (
            display_table["percentage"]
            .round(2)
        )

        st.subheader(
            "Cell Frequency Table"
        )

        st.dataframe(
            display_table,
            use_container_width=True,
            hide_index=True
        )

        bar_plot = px.bar(
            selected_data,
            x="sample",
            y="percentage",
            color="population",
            title=(
                "Relative Frequency of Immune-Cell "
                "Populations by Sample"
            ),
            labels={
                "sample": "Sample",
                "percentage": "Relative frequency (%)",
                "population": "Cell population"
            }
        )

        st.plotly_chart(
            bar_plot,
            use_container_width=True
        )

    else:

        st.info(
            "Please select at least one sample."
        )


# ---------------------------------------------------------
# Tab 2: Miraclib responder analysis
# ---------------------------------------------------------

with tab2:

    st.header(
        "Miraclib Responders versus Non-Responders"
    )

    st.write(
        """
        This analysis includes only melanoma PBMC samples
        from patients treated with miraclib.
        """
    )

    response_data = frequency_data[
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

    number_of_responder_samples = (
        response_data[
            response_data["response"] == "yes"
        ]["sample"]
        .nunique()
    )

    number_of_nonresponder_samples = (
        response_data[
            response_data["response"] == "no"
        ]["sample"]
        .nunique()
    )

    metric1, metric2, metric3 = st.columns(3)

    metric1.metric(
        "Responder samples",
        number_of_responder_samples
    )

    metric2.metric(
        "Non-responder samples",
        number_of_nonresponder_samples
    )

    metric3.metric(
        "Cell populations tested",
        response_data["population"].nunique()
    )

    box_plot = px.box(
        response_data,
        x="population",
        y="percentage",
        color="response",
        points="outliers",
        title=(
            "Relative Cell Frequencies in Responders "
            "and Non-Responders"
        ),
        labels={
            "population": "Immune-cell population",
            "percentage": "Relative frequency (%)",
            "response": "Response"
        }
    )

    st.plotly_chart(
        box_plot,
        use_container_width=True
    )

    st.subheader(
        "Statistical Results"
    )

    statistics_display = statistics_data.copy()

    numeric_columns = [
        "responder_mean_percentage",
        "nonresponder_mean_percentage",
        "difference_in_means",
        "mann_whitney_u",
        "p_value",
        "adjusted_p_value"
    ]

    for column in numeric_columns:
        if column in statistics_display.columns:
            statistics_display[column] = (
                statistics_display[column]
                .round(4)
            )

    st.dataframe(
        statistics_display,
        use_container_width=True,
        hide_index=True
    )

    significant_results = statistics_data[
        statistics_data["significant"] == True
    ]

    if len(significant_results) == 0:

        st.info(
            """
            No immune-cell populations showed a statistically
            significant difference after Benjamini-Hochberg
            multiple-testing correction.
            """
        )

    else:

        st.success(
            "Significant immune-cell populations were found."
        )

        st.dataframe(
            significant_results,
            use_container_width=True,
            hide_index=True
        )

    st.caption(
        """
        A two-sided Mann-Whitney U test was used for each
        population. P-values were adjusted using the
        Benjamini-Hochberg false-discovery-rate method.
        """
    )


# ---------------------------------------------------------
# Tab 3: Baseline subset analysis
# ---------------------------------------------------------

with tab3:

    st.header(
        "Baseline Melanoma PBMC Samples Treated with Miraclib"
    )

    st.write(
        """
        This section includes melanoma PBMC samples collected
        at time 0 from patients treated with miraclib.
        """
    )

    average_b_cells = float(
        b_cell_average[
            "average_b_cell_count"
        ].iloc[0]
    )

    baseline_subjects = baseline_data[
        [
            "project",
            "subject"
        ]
    ].drop_duplicates()

    metric1, metric2, metric3 = st.columns(3)

    metric1.metric(
        "Baseline samples",
        len(baseline_data)
    )

    metric2.metric(
        "Unique subjects",
        len(baseline_subjects)
    )

    metric3.metric(
        "Average B-cell count",
        f"{average_b_cells:.2f}"
    )

    st.caption(
        """
        The B-cell average is for male melanoma responders
        at time 0 across all sample types and treatment types.
        """
    )

    left_column, right_column = st.columns(2)

    with left_column:

        st.subheader(
            "Samples from Each Project"
        )

        st.dataframe(
            project_counts,
            use_container_width=True,
            hide_index=True
        )

        project_plot = px.bar(
            project_counts,
            x="project",
            y="sample_count",
            title="Baseline Samples by Project",
            labels={
                "project": "Project",
                "sample_count": "Number of samples"
            }
        )

        st.plotly_chart(
            project_plot,
            use_container_width=True
        )

    with right_column:

        st.subheader(
            "Subjects by Response"
        )

        st.dataframe(
            response_counts,
            use_container_width=True,
            hide_index=True
        )

        response_plot = px.bar(
            response_counts,
            x="response",
            y="subject_count",
            title="Baseline Subjects by Response",
            labels={
                "response": "Response",
                "subject_count": "Number of subjects"
            }
        )

        st.plotly_chart(
            response_plot,
            use_container_width=True
        )

    st.subheader(
        "Subjects by Gender"
    )

    gender_column1, gender_column2 = st.columns(
        [1, 2]
    )

    with gender_column1:

        st.dataframe(
            gender_counts,
            use_container_width=True,
            hide_index=True
        )

    with gender_column2:

        gender_plot = px.bar(
            gender_counts,
            x="gender",
            y="subject_count",
            title="Baseline Subjects by Gender",
            labels={
                "gender": "Gender",
                "subject_count": "Number of subjects"
            }
        )

        st.plotly_chart(
            gender_plot,
            use_container_width=True
        )

    st.subheader(
        "Baseline Sample Records"
    )

    st.dataframe(
        baseline_data,
        use_container_width=True,
        hide_index=True
    )


# ---------------------------------------------------------
# Tab 4: Filter and explore the full data
# ---------------------------------------------------------

with tab4:

    st.header(
        "Explore the Complete Frequency Dataset"
    )

    st.write(
        """
        Use the filters below to explore specific projects,
        indications, treatments, sample types, responses,
        and cell populations.
        """
    )

    filter_column1, filter_column2 = st.columns(2)

    with filter_column1:

        selected_projects = st.multiselect(
            "Project",
            options=sorted(
                frequency_data["project"].unique()
            )
        )

        selected_indications = st.multiselect(
            "Indication",
            options=sorted(
                frequency_data["indication"].unique()
            )
        )

        selected_treatments = st.multiselect(
            "Treatment",
            options=sorted(
                frequency_data["treatment"].unique()
            )
        )

    with filter_column2:

        selected_sample_types = st.multiselect(
            "Sample type",
            options=sorted(
                frequency_data["sample_type"].unique()
            )
        )

        available_responses = (
            frequency_data["response"]
            .dropna()
            .unique()
        )

        selected_responses = st.multiselect(
            "Response",
            options=sorted(
                available_responses
            )
        )

        selected_populations = st.multiselect(
            "Cell population",
            options=sorted(
                frequency_data["population"].unique()
            )
        )

    filtered_data = frequency_data.copy()

    if selected_projects:

        filtered_data = filtered_data[
            filtered_data["project"].isin(
                selected_projects
            )
        ]

    if selected_indications:

        filtered_data = filtered_data[
            filtered_data["indication"].isin(
                selected_indications
            )
        ]

    if selected_treatments:

        filtered_data = filtered_data[
            filtered_data["treatment"].isin(
                selected_treatments
            )
        ]

    if selected_sample_types:

        filtered_data = filtered_data[
            filtered_data["sample_type"].isin(
                selected_sample_types
            )
        ]

    if selected_responses:

        filtered_data = filtered_data[
            filtered_data["response"].isin(
                selected_responses
            )
        ]

    if selected_populations:

        filtered_data = filtered_data[
            filtered_data["population"].isin(
                selected_populations
            )
        ]

    st.write(
        "Number of matching rows:",
        len(filtered_data)
    )

    st.dataframe(
        filtered_data,
        use_container_width=True,
        hide_index=True
    )

    csv_download = filtered_data.to_csv(
        index=False
    ).encode("utf-8")

    st.download_button(
        label="Download Filtered Data",
        data=csv_download,
        file_name="filtered_cell_population_data.csv",
        mime="text/csv"
    )


# ---------------------------------------------------------
# Sidebar information
# ---------------------------------------------------------

st.sidebar.header(
    "About This Dashboard"
)

st.sidebar.write(
    """
    This dashboard was created using Python,
    pandas, SQLite, Plotly, and Streamlit.
    """
)

st.sidebar.write(
    """
    The database design can support additional treatments
    without changing the schema. For example, a future
    treatment such as quintazide could be added as another
    treatment value and would automatically appear in the
    treatment filters.
    """
)
