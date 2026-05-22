#!/usr/bin/env python3
"""Early-2020 Washington State COVID-19 daily cases -> TSVs for the seed-timing
chapter. Run from anywhere: writes next to this script.

Source: NYT covid-19-data (us-states.csv), whose US series begins with the
first US case — in Washington, 2020-01-21. We diff the cumulative count to a
daily incidence series, anchored at t = 0 on 2020-01-21.

Produces three files:
  covid_wa_daily.tsv      full daily incidence, t=0 = 2020-01-21, through Apr 2020
  covid_wa_community.tsv  daily with the t=0 import singleton dropped (the
                          21 Jan travel case is not the community-clade founder;
                          a smooth seed gives it ~0 probability -> loglik -inf)
  covid_wa_growth.tsv     community series truncated to the pre-lockdown growth
                          window, days 0-55 (WA interventions ramped 11-23 Mar)

CAVEAT (stated loudly in the chapter): the *detected* case curve is dominated by
testing ramp-up and reporting delay precisely in the early window that carries
the seed time. The first community-acquired case was detected 2020-02-28 after
an estimated 4-6 weeks of cryptic transmission; the genomic introduction
estimate (outbreak-clade common ancestor / tMRCA) is 2020-01-18 to 2020-02-09
(Bedford et al., Cryptic transmission of SARS-CoV-2 in Washington State,
Science 2020, doi:10.1126/science.abc0523). Those genomic numbers — not the
case counts — are the trustworthy anchor for tau.
"""
import io
import urllib.request
from pathlib import Path

import polars as pl

URL = "https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv"
T0 = "2020-01-21"
HERE = Path(__file__).resolve().parent

raw = urllib.request.urlopen(URL, timeout=60).read().decode()
df = (
    pl.read_csv(io.StringIO(raw))
    .filter((pl.col("state") == "Washington") & (pl.col("date") < "2020-05-01"))
    .with_columns(pl.col("date").str.to_date())
    .sort("date")
)
t0 = df["date"].min()
daily = df.with_columns(
    (pl.col("date") - t0).dt.total_days().alias("time"),
    pl.col("cases").diff().fill_null(pl.col("cases")).clip(0).alias("cases_daily"),
).select(time="time", cases="cases_daily")

daily.write_csv(HERE / "covid_wa_daily.tsv", separator="\t")
community = daily.filter(pl.col("time") >= 1)
community.write_csv(HERE / "covid_wa_community.tsv", separator="\t")
community.filter(pl.col("time") <= 55).write_csv(HERE / "covid_wa_growth.tsv", separator="\t")

print(f"wrote covid_wa_{{daily,community,growth}}.tsv to {HERE} (t=0 = {T0})")
print(f"  daily: {daily.height} d, {int(daily['cases'].sum())} cases")
print("  community drops the 21 Jan import; growth truncates to days 0-55 (pre-lockdown)")
print("  detected community onset 2020-02-28; genomic intro window ~2020-01-18..02-09")
