import numpy as np

import plotly.graph_objects as go
import plotly.express as px

import matplotlib.pyplot as plt

import seaborn as sns
import plotly.subplots as plotly_sp
import plotly.figure_factory as ff
from pyDOE2 import lhs
import pandas as pd
from scipy import stats


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


def plot_properties(df, measured, calculated, title=None, metrics_df=None, property='Rs', log_axis=True):
    colorsList = ["red", "blue", "green", "purple", "orange", "black", 'cyan']

    fig = plotly_sp.make_subplots(
        rows=3, cols=2,
        shared_xaxes='rows', shared_yaxes='rows',
        # vertical_spacing=0.03,
        subplot_titles=('Linear Scale', 'Log Scale'),
        specs=[[{"type": "scatter", 'rowspan': 2}, {"type": "scatter", 'rowspan': 2}],
               [None, None],
               [{"type": "table"}, {}]
               ]
    )

    for i_scale, scale in enumerate(['linear', 'log'], 1):
        if i_scale == 1:
            showlegend = True
        else:
            showlegend = False

        legendGroup_counter = 0
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
                                         marker={'color': colorsList[i], 'symbol': symbol},
                                         legendgroup=f'group{legendGroup_counter}', showlegend=showlegend),
                              row=1, col=i_scale
                              )
                legendGroup_counter += 1

        # add 45 line
        min_x = df[measured].min().min()
        max_x = df[measured].max().max()

        x_45 = np.linspace(min_x, max_x)
        fig.add_trace(go.Scatter(x=x_45, y=x_45,
                                 name='Perfect Trend',
                                 line=dict(color='black', dash='dash'),
                                 legendgroup=f'group{legendGroup_counter}', showlegend=showlegend),
                      row=1, col=i_scale)

        fig.update_xaxes(type=scale, row=1, col=i_scale)
        fig.update_yaxes(type=scale, row=1, col=i_scale)

    # Add table with metrics
    if metrics_df is not None:
        metrics_df = metrics_df.reset_index(names='Method')

        fig.add_trace(
            go.Table(header=dict(values=list(metrics_df.columns),
                                 # fill_color='paleturquoise',
                                 align='left'),
                     cells=dict(values=metrics_df.transpose().values.tolist(),
                                # fill_color='lavender',
                                align='left'),
                     columnwidth=[400, 80, 80, 80],
                     ), row=3, col=1
        )
        # fig.add_trace(ff.create_table(metrics_df, index=True), row=2, col=1)

    fig.update_xaxes(title_text="Measured",
                     minor=dict(ticks="inside", showgrid=True))

    fig.update_yaxes(  # scaleanchor="x",
        scaleratio=1,
        title_text="Calculated",
    )

    fig.update_layout(
        title=dict(text=title, font=dict(size=50), automargin=True)
    )

    fig.write_html(fr"figures/{property}.html")
    fig.show()


def plot_comparePVT(inputs, df_old, df_new, df_opt=None, x_axis='p', title='', path='', properties=None):
    if properties is None:
        properties = list(df_new.columns)

    n_properties = len(properties)

    fig = plotly_sp.make_subplots(
        rows=2, cols=n_properties,
        shared_xaxes='rows', shared_yaxes='columns',
        # vertical_spacing=0.03,
        subplot_titles=properties + properties,
    )

    showlegend = True

    legends = ['Lab Value', 'Old PVT', 'New PVT', 'Matched PVT']

    for i_scale, scale in enumerate(['linear', 'log'], 1):
        for i_prop, property_i in enumerate(properties, 1):

            fig.add_trace(go.Scatter(mode="markers", x=inputs[x_axis], y=df_old[property_i],
                                     name=legends[1],
                                     marker={'color': "red", 'symbol': 'square'},
                                     legendgroup='group2', showlegend=showlegend),
                          row=i_scale, col=i_prop
                          )

            fig.add_trace(go.Scatter(mode="markers", x=inputs[x_axis], y=df_new[property_i],
                                     name=legends[2],
                                     marker={'color': "blue", 'symbol': 'x'},
                                     legendgroup='group3', showlegend=showlegend),
                          row=i_scale, col=i_prop
                          )

            if df_opt is not None:
                fig.add_trace(go.Scatter(mode="markers", x=inputs[x_axis], y=df_opt[property_i],
                                         name=legends[3],
                                         marker={'color': "green", 'symbol': 'circle-cross'},
                                         legendgroup='group4', showlegend=showlegend),
                              row=i_scale, col=i_prop
                              )

            if property_i in inputs:
                fig.add_trace(go.Scatter(mode="markers", x=inputs[x_axis], y=inputs[property_i],
                                         name=legends[0],
                                         marker={'color': "black", 'symbol': 'circle-open'},
                                         legendgroup='group1', showlegend=showlegend),
                              row=i_scale, col=i_prop
                              )

            fig.update_xaxes(type=scale, title_text=x_axis, row=i_scale, col=i_prop)
            fig.update_yaxes(type=scale, row=i_scale, col=i_prop)  # , title_text=property_i

            if i_prop == 1:
                fig.update_yaxes(title_text=f'<b>{scale}</b>', row=i_scale, col=i_prop)

            if showlegend:
                showlegend = False

    # fig.update_xaxes(minor=dict(ticks="inside", showgrid=True))

    fig.update_layout(
        title=dict(text=f'{title}:  Linear at the top | Log at the bottom', font=dict(size=50),
                   automargin=True)
    )

    fig.write_html(fr"{path}comparingPVTs_{title}.html")
    fig.show()


def plot_pairplots(df, hue='', origin='xom'):
    g1 = sns.pairplot(df, hue=hue)
    g1.figure.savefig(rf'figures/pairplots_{origin}.png')

    g2 = sns.pairplot(df, y_vars='Rs', hue=hue)
    g2.figure.savefig(rf'figures/pairplots_Rs_{origin}.png')


def metrics(measured, calculated, columns=None):
    if columns is not None:
        measured = measured[columns].to_numpy()
        calculated = calculated[columns].to_numpy()

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


def sampling(sampling_type='lhs', nVariables=2, n_samples=100, criterion=None,
             random_state=123, n_psat=50, bounds=None):
    keys = list(bounds.keys())
    keys_expanded = keys + ['sample']

    p_bounds = bounds['p']

    bounds = np.vstack((bounds['API'], bounds['gamma_s'], bounds['temperature'])).T

    if sampling_type == 'lhs':
        X = lhs(nVariables, samples=n_samples, random_state=random_state)

        for i in range(nVariables):
            min_val = bounds[0, i]
            max_val = bounds[1, i]
            X[:, i] = X[:, i] * (max_val - min_val) + min_val

        keys.remove('p')
        dfX = pd.DataFrame(X, columns=keys)
        # print(f'Sampled values for Synthetic case: \n{dfX}')

    else:
        raise Exception(f'Sampling method {sampling_type} ot defined')

    # include pressure
    psat = np.linspace(start=p_bounds[0], stop=p_bounds[1], num=n_psat)
    sampleNumber = np.reshape(np.arange(n_samples, dtype=int), (-1, 1))

    # combinatorial
    # todo: improve Code Quality
    Y = []
    for i in range(n_psat):
        psat_i = np.full((n_samples, 1), psat[i])
        Y.append(np.hstack((psat_i, X, sampleNumber)))

    df = pd.DataFrame(np.concatenate(Y), columns=keys_expanded)
    df = df.sort_values(by='sample').reset_index(drop=True)

    df['sample'] = df['sample'].astype(int)
    return df, dfX


def plot_synthetic_data(correlations_df, input_df, name='', jumpLog='', hueplot=None):
    palette = None
    nplots = len(correlations_df.columns)

    fig, axes = plt.subplots(nrows=nplots, ncols=2, figsize=(10, 12))

    for j, axes_type in enumerate(['linear', 'log']):
        for i, correlation_name in enumerate(correlations_df):
            data = input_df[['p', 'sample']]
            data[correlation_name] = correlations_df[correlation_name]
            if (j == 1) and (i == nplots - 2):
                showlegend = True
            else:
                showlegend = False

            if hueplot is not None:
                palette = sns.color_palette()

            g = sns.lineplot(data=data,
                             x="p",
                             y=correlation_name,
                             hue=hueplot,
                             palette=palette,
                             ax=axes[i, j],
                             legend=showlegend
                             )

            axes[i, j].set_title(correlation_name)
            axes[i, j].set_ylabel('')

            if (j == 1) and (correlation_name == jumpLog):
                axes[i, j].set_xlim(0, 5000)
            else:
                axes[i, j].set_yscale(axes_type)
                axes[i, j].set_xscale(axes_type)
    fig.suptitle(name, fontsize=24)
    plt.tight_layout()
    fig.savefig(fr'figures\synthetic_{name}_hue_{hueplot}')
    plt.show()


def relativeErrorforMatch(dfold, dfnew, columns=None, type_error='abs'):
    if columns is None:
        columns = dfnew.columns

    n_samples = dfnew.shape[0]

    old_values = dfold[columns].to_numpy()
    new_values = dfnew[columns].to_numpy()

    if type_error == 'square':
        error = ((new_values - old_values) / old_values) ** 2
    elif type_error == 'abs':
        error = np.abs(((new_values - old_values) / old_values))
    else:
        raise Exception(f'Optimizer error method not implemented: {type_error}')

    error_prop = np.sum(error, axis=0) / n_samples
    error_total = np.sum(error_prop)

    df_error = pd.DataFrame(error_prop.reshape(1, -1), columns=columns)

    return error_total, df_error


def printInputValues(old=None, new=None):
    columns = ['API', 'Specific_gravity', 'Temperature']

    if old is None:
        rows = ['Calculated_Values']
        new = np.reshape(new, (1, -1))
    else:
        if isinstance(old, pd.DataFrame):
            old = old[['API', 'gamma_s', 'temperature']].iloc[0].values
        rows = ['Original_Values', 'Calculated_Values']
        new = np.vstack((old, new))

    df = pd.DataFrame(new, columns=columns, index=rows)

    print(f'Result from the match:\n{df}')


def concatDF(df1, df2):
    df3 = pd.concat([df1, df2])
    return df3


def metric2df(dict_, errorObj=None):
    df = pd.DataFrame.from_dict(dict_).reset_index(names='property')
    df_long = pd.melt(df, id_vars=['property'], var_name='metric')
    df_out = df_long.set_index(['metric', 'property']).T

    if errorObj is not None:
        df_out.insert(0, 'optError', errorObj)

    return df_out

