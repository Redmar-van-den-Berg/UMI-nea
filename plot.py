# /// script
# [tool.marimo.runtime]
# auto_instantiate = false
# ///

import marimo

__generated_with = "0.19.11"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import altair as alt

    return alt, mo, pd


@app.cell
def _(pd):
    df = pd.read_csv("perf", sep=" ")
    return (df,)


@app.cell
def _(df, mo):
    # Create a simple UI to select the number of points for a random scatter plot
    fields = [
        ("num_founder", int),
        ("mean_children_num", int),
        ("umi_len", int),
        ("insertion-deletion-substitution", str),
        ("err_rate", float),
    
    ]
    filter = mo.ui.dictionary(
        { field: mo.ui.dropdown(options=df[field], value=type(df[field][0])) for field, type in fields}
    )
    filter
    return fields, filter


@app.cell
def _(filter):
    filter
    return


@app.cell
def _(df, fields, filter):
    # Determine which rows are all true for each specified filter
    series = [df[field] == filter.value[field] for field, type in fields]
    all_true = list()
    nr_items = len(series[0])
    for i in range(len(series[0])):
        a = all((serie[i] for serie in series))
        all_true.append(a)

    filtered=df.loc[all_true]
    # Remove NaN columns
    filtered = filtered.drop('RPU_cutoff', axis=1)
    filtered = filtered.drop('RPU_cutoff_model', axis=1)
    return (filtered,)


@app.cell
def _(filter):
    filter
    return


@app.cell
def _(alt, filtered):
    def barchart(df, column, log=False):
        #column = "V-measure"
        #column = "runtime_in_sec"
        # The bars
        bars=alt.Chart(df).mark_bar().encode(
            x='tool',
            y=alt.Y(f'mean({column}):Q').title(f'Mean {column}').scale(type="log" if log else ""),
            color=alt.Color("tool")
    
        )
        # The error bars
        error_bars = alt.Chart().mark_errorbar(extent='ci').encode(
            x='tool',
            y=column,
        )
        return alt.layer(bars,error_bars, data=df)

    barchart(filtered, "V-measure", log=True)
    barchart(filtered, "runtime_in_sec")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
