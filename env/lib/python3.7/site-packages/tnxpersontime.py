import datetime as dt
from collections import defaultdict

import numpy as np
import pandas as pd
import scipy
import scipy.stats


def time_elapsed(func):
    def wrapper(*args, **kwargs):
        start = dt.datetime.now()
        out = func(*args, **kwargs)
        print(f"Time elapsed: {dt.datetime.now() - start}")
        return out

    return wrapper


class TNXDateCsv:
    def __init__(self, filepath, dt_cols=None, dt_format=None):
        self.raw_file = pd.read_csv(filepath)
        self._process_csv(dt_cols, dt_format)

    def _process_csv(self, dt_cols, dt_format):
        self.df = self.raw_file.copy()
        if dt_cols:
            if isinstance(dt_cols, str):
                dt_cols = [dt_cols]
            for col in dt_cols:
                self.df[col] = pd.to_datetime(self.df[col], format=dt_format)

    def __call__(self):
        return self.df


class TNXLabCsv(TNXDateCsv):
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
        self.encounter_file = encounter_object()
        self.index_file = index_object()

        self._make_patient_dict()

        self.encounter_file = self.encounter_file[
            [
                idx in self.patient_dict.keys()
                for idx in self.encounter_file[patient_id_alias]
            ]
        ]

    def _make_patient_dict(self):

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

        window = pd.Timedelta(value=window_days, unit="days")
        pt_id = list(set(pt[self.patient_id_alias]))[0]
        endpoint = [self.patient_dict[pt_id][e] for e in self.endpoints]
        endpoint = [e for e in endpoint if not pd.isnull(e)]
        if endpoint:
            endpoint = min(endpoint)

        def _windowed_time(row):

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
        self, window_days, index_offset, save_output=True, output_save_path=None
    ):

        if save_output:
            assert output_save_path is not None

        print(
            f"Generating person-time data for {len(self.patient_dict)} patients and {len(self.encounter_file)} encounters with {window_days} day window and {index_offset} day offset from index date."
        )
        processed = self.encounter_file.groupby(self.patient_id_alias).apply(
            self._person_time, window_days=window_days, index_offset=index_offset
        )
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
        self.person_time_df = pd.read_csv(fpath)

        try:
            unknown_pts = [idx not in self.patient_dict.keys()
                           for idx in self.person_time_df[self.patient_id_alias]]
            assert sum(unknown_pts) == 0
        except AssertionError as err:
            print(err)
            print(f'{sum(unknown_pts)} patient ids from loaded file not found in patient dictionary. Did you load the correct file?')


class RiskAnalysis():
    def __init__(self, persontime, stratification_col, outcome_alias, empty_value=''):

        self.persontime = persontime

        self._generate_tallied_df(
            stratification_col, outcome_alias, empty_value)

    def _se(self, x, type, referent):

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

    def _generate_tallied_df(self, stratification_col, outcome_alias, empty_value):

        self.merged_person_time_df = self.persontime.person_time_df.merge(self.persontime.index_file[[
            self.persontime.patient_id_alias, stratification_col]], how='left')

        self.merged_person_time_df[stratification_col] = self.merged_person_time_df[stratification_col].fillna(
            empty_value)

        stratified_times = self.merged_person_time_df.groupby(
            stratification_col).apply(lambda x: sum(x.window_time))

        stratified_cases = self.merged_person_time_df.groupby(stratification_col).apply(lambda x: sum(
            [not pd.isna(self.persontime.patient_dict[idx][outcome_alias]) for idx in x[self.persontime.patient_id_alias]]))

        stratified_persons = self.merged_person_time_df.groupby(
            stratification_col).count().iloc[:, 0]

        self.tallied_df = pd.concat([stratified_persons, stratified_cases, stratified_times], axis=1).rename(
            columns={self.persontime.patient_id_alias: 'persons', 0: 'cases', 1: 'time'})

        self.tallied_df.loc['total', :] = self.tallied_df.sum()
        self.tallied_df = self.tallied_df.astype('int')

    def _calculate_df_metrics(self, confidence_level, crude_colname, metric_type, referent, denom, order):

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

    def generate_rate_risk_dfs(self, referent, person_denom=100, rate_denom=100000, confidence_level=0.95, order=None):

        self.risk_df = self._calculate_df_metrics(
            confidence_level=confidence_level, crude_colname=f'crude_risk_per_{person_denom}', metric_type='risk', denom=person_denom, referent=referent, order=order)

        self.rate_df = self._calculate_df_metrics(
            confidence_level=confidence_level, crude_colname=f'crude_rate_per_{rate_denom}', metric_type='rate', denom=rate_denom, referent=referent, order=order)
