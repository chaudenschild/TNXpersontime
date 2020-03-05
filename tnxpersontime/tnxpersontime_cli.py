from TNXpersontime.tnxpersontime import PersonTime, RiskAnalysis, TNXDateCsv

# TODO: make cli

if __name__ == "__main__":

    encounter = TNXDateCsv("data/encounter.csv",
                           ["start_date", "end_date"], "%Y%m%d")
    index = TNXDateCsv(
        "data/index_with_outcome.csv", ["index_d",
                                        "date_death", "date_HF"], "%d%b%y"
    )

    persontime30days = PersonTime(
        encounter,
        index,
        index_file_endpoints=["date_death", "date_HF"],
        index_date_alias="index_d",
        patient_id_alias="patient_id",
    )
    persontime6mo = PersonTime(
        encounter,
        index,
        index_file_endpoints=["date_death", "date_HF"],
        index_date_alias="index_d",
        patient_id_alias="patient_id",
    )

    persontime1yr = PersonTime(
        encounter,
        index,
        index_file_endpoints=["date_death", "date_HF"],
        index_date_alias="index_d",
        patient_id_alias="patient_id",
    )

    persontime30days.generate_person_time_df(
        window_days=30, index_offset=365, save_output=True, output_save_path='data/days_30_df.csv'
    )

    persontime6mo.generate_person_time_df(
        window_days=183, index_offset=365, save_output=True, output_save_path='data/months_6_df.csv'
    )

    persontime1yr.generate_person_time_df(
        window_days=365,
        index_offset=365,
        save_output=True,
        output_save_path="data/months_12_df.csv",
    )
