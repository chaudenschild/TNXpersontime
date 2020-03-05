
tnxpersontime.venv.setup <- function(venv='tnxpersontime') {
  reticulate::virtualenv_create('./env')
  reticulate::virtualenv_install("./env", packages = c("pandas", "numpy", "scipy"))
  ret
  reticulate::use_virtualenv("./env", required = TRUE)
  tnx <<- reticulate::import('tnxpersontime')
}
reticulate::virtualenv_root()

read.tnxdatecsv <- function(filepath, dt_cols=NULL, dt_format=NULL) {
  date_csv_obj = tnx$TNXDateCsv(filepath, dt_cols=dt_cols, dt_format=dt_format)
  return(date_csv_obj)
}

read.tnxlabcsv <- function(filepath, dt_cols=NULL, dt_format=NULL) {
  lab_csv_obj = tnx$TNXLabCsv(filepath, dt_cols=dt_cols, dt_format=dt_format)
  return(date_csv_obj)
}

create.tnxpersontime <- function(encounter_object, index_object, index_file_endpoints=NULL, index_date_alias="index_d", patient_id_alias="patient_id") {
  persontime = tnx$PersonTime(encounter_object, index_object, index_file_endpoints, index_date_alias, patient_id_alias)
  return(persontime)
}

create.riskanalysis <- function(persontime, stratification_col, outcome_alias, empty_value='') {
  riskanalysis = tnx$RiskAnalysis(persontime, stratification_col, outcome_alias, empty_value)
  return(riskanalysis)
}

