import numpy as np

import plotly.graph_objects as go
import plotly.express as px

import seaborn as sns
import plotly.subplots as plotly_sp
import plotly.figure_factory as ff


def EDA_plotly(df):
    x_vars = 'p_sat'
    y_vars = ['Rs', 'Bo_psat', 'visc_o_psat']

    figs = []
    for y_var in y_vars:
        fig_px = px.scatter(df, x=x_vars, y=y_var, color='source', log_x=True, log_y=True)
        fig_px.update_layout(showlegend=False)
        figs.append(fig_px)

    # For as many traces that exist per Express figure, get the traces from each plot and store them in an array.
    # This is essentially breaking down the Express fig into it's traces
    figure_traces = {}
    for i, fig in enumerate(figs):
        for j, trace in enumerate(range(len(fig["data"]))):
            if j == 0:
                figure_traces[i] = [fig["data"][trace]]
            else:
                figure_traces[i].append(fig["data"][trace])

    # Create a 1x2 subplot
    fig = plotly_sp.make_subplots(rows=len(figs), cols=1, shared_xaxes=True)

    # Get the Express fig broken down as traces and add the traces to the proper plot within in the subplot
    for i, figure_trace in enumerate(figure_traces.values()):
        for traces in figure_trace:
            fig.add_trace(traces, row=i + 1, col=1)

        if i == 0:
            xaxis = 'xaxis'
            yaxis = 'yaxis'
        else:
            xaxis = f'xaxis{i + 1}'
            yaxis = f'yaxis{i + 1}'
        fig['layout'][xaxis]['title'] = x_vars
        fig['layout'][yaxis]['title'] = y_vars[i]

    fig.update_xaxes(type="log")
    fig.update_yaxes(type="log")

    # set x an y labels
    fig.show()


def EDA_seaborn(df):
    # plot comparing the viscosities
    g = sns.pairplot(df, vars=['gamma_s', 'gamma_c'], hue='source', )
    g.figure.savefig(rf'figures/pairplot_specific_gravity.png')

    # plot psat vs Rs,  psat vs Bo, psat vs visco
    g = sns.pairplot(df, x_vars='p_sat', y_vars=['Rs', 'Bo_psat', 'visc_o_psat'], hue='source',
                     height=3)
    g.figure.savefig(rf'figures/pairplot_p_sat.png')

    # Check correlations all
    g = sns.pairplot(df, vars=['p_sat', 'temperature', 'gamma_s', 'gamma_c', 'API', 'Rs'], hue='HGOR')
    g.figure.savefig(rf'figures/pairplots_all.png')

    # Check correlations with RS
    g = sns.pairplot(df, x_vars=['p_sat', 'temperature', 'gamma_s', 'gamma_c', 'API', 'Rs'],
                     y_vars='Rs',
                     hue='HGOR')
    g.figure.savefig(rf'figures/pairplots_Rs.png')


def plot_log_log(df, measured, calculated, title=None, metrics_df=None, property='Rs'):
    colorsList = ["red", "blue", "green", "purple", "orange", "black"]

    fig = plotly_sp.make_subplots(
        rows=3, cols=1,
        # vertical_spacing=0.03,
        specs=[[{"type": "scatter", 'rowspan': 2}],
               [None],
               [{"type": "table"}]
               ]
    )

    for i, method in enumerate(calculated):

        for hgor, df_gor in df.groupby('HGOR'):

            if hgor:
                name = method + '_H_gor'
                symbol = 'diamond'
            else:
                name = method + '_l_gor'
                symbol = 'x'

            fig.add_trace(go.Scatter(mode="markers", x=df_gor[measured], y=df_gor[method],
                                     name=name,
                                     marker={'color': colorsList[i], 'symbol': symbol}),
                          row=1, col=1
                          )

    # add 45 line
    min_x = df[measured].min().min()
    max_x = df[measured].max().max()

    x_45 = np.linspace(min_x, max_x)
    fig.add_trace(go.Scatter(x=x_45, y=x_45,
                             name='Perfect Trend',
                             line=dict(color='black', dash='dash')), row=1, col=1)

    # Add table with metrics
    if metrics_df is not None:
        metrics_df = metrics_df.reset_index(names='Method')

        fig.add_trace(
            go.Table(header=dict(values=list(metrics_df.columns),
                                 # fill_color='paleturquoise',
                                 align='left'),
                     cells=dict(values=metrics_df.transpose().values.tolist(),
                                # fill_color='lavender',
                                align='left')
                     ), row=3, col=1
        )
        # fig.add_trace(ff.create_table(metrics_df, index=True), row=2, col=1)

    fig.update_xaxes(type="log",
                     title_text="Measured",
                     minor=dict(ticks="inside", showgrid=True))

    fig.update_yaxes(type="log",
                     # scaleanchor="x",
                     scaleratio=1,
                     title_text="Calculated",
                     )

    fig.update_layout(
        title=dict(text=title, font=dict(size=50), automargin=True)
    )

    fig.show()
    fig.write_html(fr"figures/{property}.html")


def plot_pairplots(df, hue='', origin='xom'):
    g1 = sns.pairplot(df, hue=hue)
    g1.figure.savefig(rf'figures/pairplots_{origin}.png')

    g2 = sns.pairplot(df, y_vars='Rs', hue=hue)
    g2.figure.savefig(rf'figures/pairplots_Rs_{origin}.png')


def metrics(measured, calculated):
    # remove nan
    nan_measured = ~np.isnan(measured)
    measured = measured[nan_measured]
    calculated = calculated[nan_measured]

    ln_measured = np.log(measured)
    ln_calculated = np.log(calculated)

    n_samples = measured.shape[0]

    ADE = np.sum(np.abs(ln_measured - ln_calculated))
    LSE = np.sum(np.power(ln_measured - ln_calculated, 2))
    AARE = np.sum(np.abs((measured - calculated) / calculated)) * 100 / n_samples

    metrics_ = {'ADE': ADE, 'LSE': LSE, 'AARE': AARE}

    return metrics_
