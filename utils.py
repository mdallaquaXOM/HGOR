import numpy as np

from statsmodels.nonparametric.smoothers_lowess import lowess
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.pyplot as plt
import pandas as pd
import os
import seaborn as sns
import plotly.subplots as plotly_sp
import plotly.figure_factory as ff
from tqdm import tqdm


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


def plot_log_log(df, measured, calculated, title, metrics_df=None):
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
    fig.write_html(fr"figures/Rs.html")


def plot_pairplots(df, hue='', origin='xom'):
    g1 = sns.pairplot(df, hue=hue)
    g1.figure.savefig(rf'figures/pairplots_{origin}.png')

    g2 = sns.pairplot(df, y_vars='Rs', hue=hue)
    g2.figure.savefig(rf'figures/pairplots_Rs_{origin}.png')


def metrics(measured, calculated):
    ln_measured = np.log(measured)
    ln_calculated = np.log(calculated)

    n_samples = measured.shape[0]

    ADE = np.sum(np.abs(ln_measured - ln_calculated))
    LSE = np.sum(np.power(ln_measured - ln_calculated, 2))
    AARE = np.sum(np.abs((measured - calculated) / calculated)) * 100 / n_samples

    metrics = {'ADE': ADE, 'LSE': LSE, 'AARE': AARE}

    return metrics


class Grace:
    def __init__(self):
        self.tol = 1e-8
        self._span = None

    @property
    def span(self):
        return self._span

    @span.setter
    def span(self, value):
        self._span = value

    def finite_ace(self, x, y):
        x = (x - np.mean(x)) / np.std(x)
        y = (y - np.mean(y)) / np.std(y)

        theta_y = y / np.linalg.norm(y, 2)

        phi_x = 0.
        delta = 1
        error0 = np.ones((x.shape[0], 1))
        count = 1
        tol = self.tol
        while delta > tol:
            phi_x = self.smoother(x=x, y=theta_y, span=self.span)
            theta_y = self.smoother(x=y, y=phi_x, span=self.span)

            theta_y = (theta_y - np.mean(theta_y)) / np.linalg.norm(theta_y, 2)

            error = (theta_y - phi_x) ** 2
            delta = np.linalg.norm(error - error0, 2)
            error0 = error

            if count > 50:
                tol = tol * 10
                count = 1
            else:
                count += 1

        return phi_x, theta_y

    @staticmethod
    def smoother(x, y, span):
        # s = savgol_filter(x, window_length=span, polyorder=polyorder)
        s_ = lowess(endog=y, exog=x, frac=span, return_sorted=False)
        return s_

    @staticmethod
    def linearRegression(phi_x, theta_y):
        # linear regression
        a = phi_x.T @ phi_x
        b = phi_x.T @ theta_y
        beta = np.linalg.solve(a, b)

        yCalc = phi_x @ beta

        SSE = np.sum((theta_y - yCalc) ** 2)
        Syy = np.sum((theta_y - np.mean(theta_y)) ** 2)
        Rsq = 1 - SSE / Syy

        return beta, Rsq

    @staticmethod
    def plot_smoother():
        # plot x vs y
        # plot x vs s
        a = 0

    @staticmethod
    def plot_ace():
        a = 0

    @staticmethod
    def plot_cv_span():
        a = 0

    def crossValidation(self, x, y, n_spans=50):
        n_samples = x.shape[0]

        k_test = np.linspace(0.01, 1., num=n_spans)

        cv_error = np.zeros((n_spans, 1))
        for i_k, k in enumerate(k_test):
            y_pred = np.zeros(x.shape)
            for i in tqdm(range(n_samples), desc=f"Testing span {i_k + 1}/{n_spans}"):
                # excluding i
                x_test = np.delete(x, i, axis=0)
                y_test = np.delete(y, i, axis=0)

                # smoothing the remaining
                s = self.smoother(x_test, y_test, span=k)

                # predicting x_i by linear interpolation
                y_pred[i] = np.interp(x[i], x_test, s)

            cv_error[i_k] = np.sum((y - y_pred) ** 2) / n_samples

        min_cv = np.argmin(cv_error)
        best_span = k_test[min_cv]

        print(f'Best Span: {best_span}')

        ## plot CV
        fig, ax = plt.subplots(1, 1)
        ax.scatter(k_test, cv_error)
        ax.scatter(best_span, cv_error[min_cv])

        plt.show()

        return cv_error


if __name__ == '__main__':
    # testing grace algorithm
    data = pd.read_csv(os.path.join('..', 'Data', 'syn.csv'))
    x = data['x']
    y = data['y']

    # smothered data
    grace = Grace()

    # perform cross validation to find best span
    # span = grace.crossValidation(x, y)

    grace.span = 0.19183673469387755
    s = grace.smoother(x=x, y=y, span=grace.span)

    # plotting the data and smoother
    fig, ax = plt.subplots(1, 1)
    ax.scatter(x, y, label='Original Data')
    ax.scatter(x, s, label=f'Smothered span:{grace.span}')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.grid()
    ax.legend()
    plt.show(block=True)

    # test ACE
    phi_x, theta_y = grace.finite_ace(data['x'], data['y'])

    fig, ax = plt.subplots(3, 1)
    ax[0].scatter(phi_x, theta_y)
    ax[0].set_xlabel('phi_x')
    ax[0].set_ylabel('theta_y')

    ax[1].scatter(x, phi_x)
    ax[1].set_xlabel('x')
    ax[1].set_ylabel('phi_x')

    ax[2].scatter(y, theta_y)
    ax[2].set_xlabel('y')
    ax[2].set_ylabel('theta_y')

    plt.tight_layout()
    plt.show(block=True)
