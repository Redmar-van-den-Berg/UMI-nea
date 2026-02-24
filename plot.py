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
def _(mo):
    f_browser=mo.ui.file_browser(multiple=False)
    f_browser
    return (f_browser,)


@app.cell
def _(f_browser, pd):
    fname=f_browser.path(index=0)
    fname=fname if fname else "perf"
    df = pd.read_csv(fname, sep=" ")
    # Add CPU time
    df["cputime_in_sec"] = df["runtime_in_sec"] * df["thread"]
    # Remove NaN columns
    df = df.drop('RPU_cutoff', axis=1)
    df = df.drop('RPU_cutoff_model', axis=1)
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
    return (filtered,)


@app.cell
def _(filter):
    filter
    return


@app.cell
def _(alt, filtered):
    def barchart(df, column, log=False, domain=None):
        #column = "V-measure"
        #column = "runtime_in_sec"
        # The bars
        bars=alt.Chart(df).mark_bar().encode(
            x='tool',
            y=alt.Y(f'mean({column}):Q').title(f'Mean {column}').scale(type="log" if log else "", domain=None),
            color=alt.Color("tool")

        )
        # The error bars
        error_bars = alt.Chart().mark_errorbar(extent='ci').encode(
            x='tool',
            y=column,
        )
        return alt.layer(bars,error_bars, data=df)

    barchart(filtered, "V-measure", domain =[0.9,1])
    return (barchart,)


@app.cell
def _(barchart, filtered):
    barchart(filtered, "runtime_in_sec")
    return


@app.cell
def _(barchart, filtered):
    barchart(filtered, "cputime_in_sec")
    return


@app.cell
def _(mo):
    settings = mo.ui.dictionary({
        "field": mo.ui.dropdown(options=["V-measure", "runtime_in_sec", "cputime_in_sec"], value="cputime_in_sec"),
        "mutation": mo.ui.dropdown(options=["1-1-1", "1-1-40"], value="1-1-1"),
        "founders": mo.ui.dropdown(options=[1000, 10000], value=1000),
    })
    settings
    return (settings,)


@app.cell
def _(alt, df, settings):
    field=settings.value["field"]
    founders=settings.value["founders"]

    t = df.loc[df["insertion-deletion-substitution"] == settings.value["mutation"]]
    t = t.loc[t["num_founder"] == settings.value["founders"]]

    size=200

    alt.Chart(t, width=size, height=size).mark_bar().encode(
        x="tool",
        y=f"mean({field}):Q",
        color=alt.Color("tool"),
        row="err_rate",
        column="umi_len",
        tooltip=[field, "insertion-deletion-substitution"]
    ).properties(title=f"Metrics for simulated data with {founders} original reads and 200 children")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
