import datetime as dt
import functools
from collections import defaultdict
from multiprocessing import cpu_count

import dill
import numpy as np
import pandas as pd
import scipy
import scipy.stats
from dask import dataframe as dd


def time_elapsed(func):
    """Decorator for timing a function"""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start = dt.datetime.now()
        out = func(*args, **kwargs)
        print(f"Time elapsed: {dt.datetime.now() - start}")
        return out

    return wrapper


class TNXCsv:
    """Base class for reading in TNX data and converting specified columns to datetime format
    Args:
        filepath (str): filepath to csv
        dt_cols (list of str, optional): list of columns to be converted (default: None)
        dt_format (str, optional): string format to use for dt_cols, based on strftime. E.g. %y%m%d (default: None)
    """

    def __init__(self, filepath, dt_cols=None, dt_format=None):
        self.raw_file = pd.read_csv(filepath)
        self._process_csv(dt_cols, dt_format)

    @classmethod
    def create(cls, filepath, dt_cols=None, dt_format=None, *args, **kwargs):
        """Class method that creates and returns a dataframe directly
        Args:
            filepath (str): filepath to csv
            dt_cols (list of str, optional): list of columns to be converted (default: None)
            dt_format (str, optional): string format to use for dt_cols, based on strftime. E.g. %y%m%d (default: None)
            *args, **kwargs: additional arguments to be passed to TNXCsv subclass constructors
        """
        obj = cls(filepath, dt_cols=None, dt_format=None, *args, **kwargs)
        return obj.df

    def _process_csv(self, dt_cols, dt_format):
        self.df = self.raw_file.copy()
        if dt_cols:
            if isinstance(dt_cols, str):
                dt_cols = [dt_cols]
            for col in dt_cols:
                self.df[col] = pd.to_datetime(self.df[col], format=dt_format)

    def __call__(self):
        return self.df


class TNXLabCsv(TNXCsv):

    def __init__(self, filepath, dt_cols=None, dt_format=None, included_lab_codes='all', code_alias='code'):

        self.raw_file = pd.read_csv(filepath)

        if included_lab_codes != 'all':
            self.raw_file = self.raw_file[[
                c in included_lab_codes for c in self.raw_file[code_alias]]]

        self._process_csv(dt_cols, dt_format)


class TNXProcedureCsv(TNXLabCsv):
    def __init__(self, filepath, dt_cols=None, dt_format=None, included_procedure_codes='all', code_alias='code'):
        super().__init__(filepath, dt_cols, dt_format,
                         included_lab_codes=included_procedure_codes, code_alias=code_alias)


class PersonTime:
    """Performs person-time calculations with encounter and index objects
    Args:
        encounter_object (TNXCsv, pandas DataFrame): TNX-style encounter file, either loaded into a TNXCsv or as a pandas DataFrame
        index_object (TNXCsv, pandas DataFrame): An index file containing a patient id column and an index date column at minimum, with optional endpoints
        index_file_endpoints (list of str, optional): list of columns in the index file that are to be considered as endpoints (default: None)
        patient_id_alias (str, optional): alternate column name for the patient id column (default: 'patient_id')

    Attributes:
        endpoints: index file endpoints
        patient_id_alias: patient id column name
        index date alias: index date column name
        index_file: pandas DataFrame containing index information
        encounter_file: subsetted encounter pandas DataFrame that includes only encounters of patients present in the index file
    """

    def __init__(
        self,
        encounter_object,
        index_object,
        index_file_endpoints=None,
        index_date_alias="index_d",
        patient_id_alias="patient_id",
    ):

        self.endpoints = index_file_endpoints
        self.patient_id_alias = patient_id_alias
        self.index_date_alias = index_date_alias
        self.encounter_file = encounter_object() if isinstance(
            encounter_object, TNXCsv) else encounter_object
        self.index_file = index_object() if isinstance(
            index_object, TNXCsv) else index_object

        self._make_patient_dict()

        self.encounter_file = self.encounter_file[
            [
                idx in self.patient_dict.keys()
                for idx in self.encounter_file[patient_id_alias]
            ]
        ]

    def _make_patient_dict(self):
        """Creates the patient_dict, automatically called on init"""

        self.patient_dict = defaultdict(lambda: {})

        for idx, index_date in zip(
            self.index_file[self.patient_id_alias], self.index_file[self.index_date_alias]
        ):
            self.patient_dict[idx][self.index_date_alias] = index_date

        if self.endpoints:
            for endpoint in self.endpoints:
                for idx, date in zip(
                    self.index_file[self.patient_id_alias], self.index_file[endpoint]
                ):
                    self.patient_dict[idx][endpoint] = date

    def _person_time(self, pt, window_days, index_offset):
        """Calculate a single patient's person-days
        Args:
            pt (pandas DataFrame): a single patient's encounters
            window_days (int): number of days in before-after window to use in the person-days calculation
            index_offset (int): amount of days to offset the index date when beginning patient time tally (e.g. 0 is equivalent to just using the index_date)
        Returns:
            window_time (int), total_time (int): person days-at-risk taking into account windowing and not taking into account windowing, respectively
        """

        window = pd.Timedelta(value=window_days, unit="days")
        pt_id = pt.name  # list(set(pt[self.patient_id_alias]))[0]
        endpoint = [self.patient_dict[pt_id][e] for e in self.endpoints]
        endpoint = [e for e in endpoint if not pd.isnull(e)]
        if endpoint:
            endpoint = min(endpoint)

        def _windowed_time(row):
            """Calculate person days of a single row in encounter file
            Args:
                row (pandas Series): single row of encounter file
            Returns:
                pandas DatetimeIndex
            """
            start = max(
                row.start_date - window,
                self.patient_dict[row.patient_id][self.index_date_alias]
                + pd.Timedelta(value=index_offset, unit="days"),
            )

            end = (
                min(
                    endpoint,
                    row.end_date + window
                    if not pd.isnull(row.end_date)
                    else row.start_date + window,
                )
                if endpoint
                else row.end_date + window
                if not pd.isnull(row.end_date)
                else row.start_date + window
            )

            return pd.date_range(start, end)

        dates = pt.apply(_windowed_time, axis=1)

        d = dates.iloc[0]
        for date in dates.iloc[1:]:
            d = d.union(date)

        window_time = len(d)

        last_encounter = endpoint if endpoint else max(
            max(pt.start_date), max(pt.end_date))
        total_time = (
            len(pd.date_range(self.patient_dict[pt_id][self.index_date_alias]
                              + pd.Timedelta(value=index_offset, unit="days"), last_encounter))
            if d.values.size > 0
            else 0
        )

        return window_time, total_time

    @time_elapsed
    def generate_person_time_df(
        self, window_days, index_offset=0, save_output=True, output_save_path=None
    ):
        """Calculates patient-by-patient days-at-risk and stores it as person_time_df attribute
            Args:
                window_days (int): number of days before and after encounter start dates and end dates to consider in the days-at-risk calculations
                index_offset (int, optional): number of days offset from the index date to begin patient days-at-risk tally (default: 0)
                save_output (bool, optional): if True, saves the output as a csv (default: True)
                output_save_path (str, optional): target save path for the output csv, if save_output is set to True (default: None)
            Returns:
                None
        """
        if save_output:
            assert output_save_path is not None

        print(
            f"Generating person-time data for {len(self.patient_dict)} patients and {len(self.encounter_file)} encounters with {window_days} day window and {index_offset} day offset from index date."
        )
        # processed = self.encounter_file.groupby(self.patient_id_alias).apply(
        #    self._person_time, window_days=window_days, index_offset=index_offset
        # )
        processed = dd.from_pandas(self.encounter_file.set_index(self.patient_id_alias), npartitions=cpu_count()).groupby(self.patient_id_alias).apply(
            self._person_time, window_days=window_days, index_offset=index_offset, meta=(None, 'O')).compute(scheduler='processes')

        out_df = pd.DataFrame(
            processed.tolist(),
            index=processed.index,
            columns=["window_time", "total_time"],
        )
        self.person_time_df = out_df

        if save_output:
            self.person_time_df.to_csv(output_save_path)
            print(f"Output saved to {output_save_path}")

    def load_person_time_df(self, fpath):
        """load a precomputed person_time_df
        Args:
            fpath (str): file path
        """
        self.person_time_df = pd.read_csv(fpath)

        try:
            unknown_pts = [idx not in self.patient_dict.keys()
                           for idx in self.person_time_df[self.patient_id_alias]]
            assert sum(unknown_pts) == 0
        except AssertionError as err:
            print(err)
            print(f'{sum(unknown_pts)} patient ids from loaded file not found in patient dictionary. Did you load the correct file?')

    def save_pkl(self, fname):
        """Pickles PersonTime object
        Args:
            fname (str): file name to save as
        """
        with open(fname, 'wb') as out:
            dill.dump(self, out)

    @staticmethod
    def read_from_pkl(fname):
        """Reads in a PersonTime object from a pickle file
        Args:
            fname (str): file name to read from
        """
        with open(fname, 'rb') as file:
            pt = dill.load(file)
        return pt

    def __getattr__(self, attr):
        if attr == 'person_time_df':
            raise AttributeError(
                'person_time_df not found - either run generate_person_time_df method or load_person_time_df')
        else:
            raise AttributeError(f'{attr} not found')


class RiskAnalysis():
    """Creates risk metrics and rate metrics using the time-at-risk from a PersonTime object
    Args:
        persontime (PersonTime): PersonTime object with person_time_df computed
        exposure_col (str): Column containing the levels of the exposure variable
        outcome_col (str): Column containing date of cases, or nan values if absent. Must be one of the endpoints specified in the PersonTime object.
        empty_value (str, optional): Name of the level for non-case patients (default: '')

    Raises:
        AssertionError if outcome_col not in PersonTime endpoints
    """

    def __init__(self, persontime, exposure_col, outcome_col, empty_value=''):

        self.persontime = persontime

        assert outcome_col in self.persontime.endpoints, f'outcome_col must be one of {self.persontime.endpoints}'

        self._generate_tallied_df(
            exposure_col, outcome_col, empty_value)

    def _se(self, x, type, referent):
        """Calculates standard error of risks and rate ratios and differences
        Args:
            x (pandas Series): A level of total persons and cases in the contingency table
            type (str): The type of desired standard error metric.
            referent: referent level of exposure
        """
        assert type in ['risk_difference', 'risk_ratio',
                        'rate_ratio', 'rate_difference']

        if x.name in [referent, 'total']:
            return np.nan

        # I - C: intervention - control aka exposed - not exposed
        # E - N: events - non-events aka cases - non-cases
        # IPt: patient time, intervention/exposed group
        # CPt: patient time, control
        IN = x.persons - x.cases
        IE = x.cases
        CN = self.tallied_df.loc[referent, 'persons'] - \
            self.tallied_df.loc[referent, 'cases']
        CE = self.tallied_df.loc[referent, 'cases']

        IPt = x.time
        CPt = self.tallied_df.loc[referent, 'time']

        if type == 'risk_ratio':
            return np.sqrt(IN / (IE * (IE + IN)) + CN / (CE * (CE + CN)))
        elif type == 'risk_difference':
            return np.sqrt(((IE * IN) / (IE + IN)**3) + ((CE * CN) / (CE + CN)**3))
        elif type == 'rate_ratio':
            return np.sqrt((1 / IE) + (1 / CE))
        elif type == 'rate_difference':
            return np.sqrt((IE / IPt**2) + (CE / CPt**2))

    def _generate_tallied_df(self, exposure_col, outcome_col, empty_value):
        """Generates count table, stratified by exposure"""

        self.merged_person_time_df = self.persontime.person_time_df.merge(self.persontime.index_file[[
            self.persontime.patient_id_alias, exposure_col]], how='left')

        self.merged_person_time_df[exposure_col] = self.merged_person_time_df[exposure_col].fillna(
            empty_value)

        stratified_times = self.merged_person_time_df.groupby(
            exposure_col).apply(lambda x: sum(x.window_time))

        stratified_cases = self.merged_person_time_df.groupby(exposure_col).apply(lambda x: sum(
            [not pd.isna(self.persontime.patient_dict[idx][outcome_col]) for idx in x[self.persontime.patient_id_alias]]))

        stratified_persons = self.merged_person_time_df.groupby(
            exposure_col).count().iloc[:, 0]

        self.tallied_df = pd.concat([stratified_persons, stratified_cases, stratified_times], axis=1).rename(
            columns={self.persontime.patient_id_alias: 'persons', 0: 'cases', 1: 'time'})

        self.tallied_df.loc['total', :] = self.tallied_df.sum()
        self.tallied_df = self.tallied_df.astype('int')

    def _calculate_df_metrics(self, confidence_level, crude_colname, metric_type, referent, denom, order):
        """Calculate various risks, rates and confidence intervals for both"""

        c_level_str = str(int(confidence_level * 100))
        alpha = abs(scipy.stats.norm.ppf((1 - confidence_level) / 2))

        df = pd.DataFrame(self.tallied_df.apply(
            lambda x, denom: x.cases / x.persons * denom if metric_type == 'risk' else x.cases / x.time * denom, denom=denom, axis=1)).rename(columns={0: crude_colname})

        df[f'crude_{metric_type}_ratio'] = df[crude_colname] / \
            df.loc[referent].values

        se_ratio = self.tallied_df.apply(
            self._se, type=f'{metric_type}_ratio', referent=referent, axis=1)

        df[f'{metric_type}_ratio_lower_bound_{c_level_str}'] = np.exp(np.log(
            df[f'crude_{metric_type}_ratio']) - alpha * se_ratio)
        df[f'{metric_type}_ratio_upper_bound_{c_level_str}'] = np.exp(np.log(
            df[f'crude_{metric_type}_ratio']) + alpha * se_ratio)

        df[f'crude_{metric_type}_difference'] = df.apply(
            lambda x: x[crude_colname] - df.loc[referent, crude_colname], axis=1)

        se_diff = self.tallied_df.apply(
            self._se, type=f'{metric_type}_difference', referent=referent, axis=1) * denom
        df[f'{metric_type}_difference_lower_bound_{c_level_str}'] = df[f'crude_{metric_type}_difference'] - alpha * se_diff
        df[f'{metric_type}_difference_upper_bound_{c_level_str}'] = df[f'crude_{metric_type}_difference'] + alpha * se_diff

        df.iloc[df.index.get_loc(referent), 1:] = np.nan
        df.iloc[df.index.get_loc('total'), 1:] = np.nan

        if order:
            df = df.loc[order]

        return df

    def generate_rate_risk_dfs(self, referent, risk_denom=100, rate_denom=100000, confidence_level=0.95, order=None):
        """Creates rate_df and risk_df attributes in the RiskAnalysis object
        Args:
            referent (str): referent level of the exposure
            risk_denom (int, optional): denominator for the risk calculations (default: 100)
            rate_denom (int, optional): denominator for rate calculations (default: 100000)
            confidence_level (float, optional): confidence interval level, between 0 and 1 (default: 0.95)
            order (list of str, optional): ordering of exposure variables in the respective dfs (default: None)
        Raises:
            ValueError if confidence level is out of bounds
        """

        if confidence_level >= 1.0 or confidence_level <= 0:
            raise ValueError('confidence_level must be between 0 and 1')

        self.risk_df = self._calculate_df_metrics(
            confidence_level=confidence_level, crude_colname=f'crude_risk_per_{person_denom}', metric_type='risk', denom=risk_denom, referent=referent, order=order)

        self.rate_df = self._calculate_df_metrics(
            confidence_level=confidence_level, crude_colname=f'crude_rate_per_{rate_denom}', metric_type='rate', denom=rate_denom, referent=referent, order=order)

    def __getattr__(self, attr):

        if attr in ['risk_df', 'rate_df']:
            raise AttributeError(
                'risk and rate not calculated, run generate_rate_risk_dfs')
        else:
            raise AttributeError(f'{attr} not found')
