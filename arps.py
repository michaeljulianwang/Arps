import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, FormatStrFormatter

def split(a, n):
    ''' split splits list 'a' into 'n' sub-lists as evenly as possible '''
    k, m = divmod(len(a), n)
    return list(a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

def Arps(Qi, D, B, dlim, years):
    '''
    Rate-based Arps

    :param Qi: initial production, bbl/day or scf/day
    :param D: D_esi, initial effective decline rate from secant approximation, 1/year
    :param B: hyperbolic exponent
    :param dlim: limiting effective decline rate, below/after which an exponential decline is used, 1/year
    :param years: number of years to project
    :return: q_daily: list with daily production rates, bbl/day or scf/day
    :return: q_monthly: list with monthly production rates, bbl/month

    References:
    https://secure.spee.org/sites/spee.org/files/wp-files/pdf/ReferencesResources/REP06-DeclineCurves.pdf
    http://www.fekete.com/san/webhelp/feketeharmony/harmony_webhelp/content/html_files/reference_material/analysis_method_theory/Traditional_Decline_Theory.htm
    https://petrowiki.org/Production_forecasting_decline_curve_analysis
    '''

    days = np.arange(0, years * 365 + 1)
    D_nom = ((1 - D) ** (-1 * B) - 1) / B # nominal decline as function of secant effective decline, 1/year
    dlim_nom = m.log(1 - dlim) / -1 # nominal decline rate at dlim, 1/year
    qlim = Qi * (dlim_nom / D_nom) ** (1 / B) # bbl/day
    tlim = m.ceil((((Qi / qlim) ** B) - 1) / (B * D_nom)) * 365 # days

    q_daily = [Qi / ((1 + B * D_nom * day / 365) ** (1 / B)) if day < tlim else qlim * m.exp(- 1 * dlim_nom * (day - tlim) / 365) for day in days]
    group_by_month = split(q_daily, years * 12) # produces nested list, grouped in subgroups of ~30 days
    q_monthly = [sum(month) for month in group_by_month]

    return q_daily, q_monthly


def cumArps(Qi, D, B, dlim, month):
    '''
    Cumulative Arps

    :param Qi: initial production, bbl/day
    :param D: D_esi, initial effective decline rate from secant approximation, 1/year
    :param B: hyperbolic exponent
    :param dlim: limiting effective decline rate, below/after which an exponential decline is used, 1/year
    :param month: integer, months to which cumulative production is calculated
    :return: q_cum: total production from t = 0 to t = months
    '''

    D_nom = ((1 - D) ** (-1 * B) - 1) / B
    dlim_nom = -1 * m.log(1 - dlim)
    qlim = Qi * ((dlim_nom / D_nom) ** (1 / B))
    tlim = m.ceil((((Qi / qlim) ** B) - 1) / (B * D_nom)) * 365 # days
    q_switch = (Qi / ((1 - B) * D_nom)) * 365 * (1 - (1 + B * D_nom * tlim / 365) ** (1 - (1 / B)))

    if month * 365 / 12 < tlim:
        q_cum = (Qi / ((1 - B) * D_nom)) * 365 * (1 - (1 + B * D_nom * month / 12) ** (1 - (1 / B)))
    else:
        q_cum = (qlim / dlim_nom * 365) * (1 - m.exp(-dlim_nom * (month * 365 / 12 - tlim) / 365)) + q_switch

    return q_cum

def visualize(q_monthly):

    fig, ax1 = plt.subplots()

    months = range(len(q_monthly))
    monthly_cum = np.cumsum(q_monthly)

    color1 = 'tab:orange'
    ax1.set_xlabel('Months')
    ax1.set_ylabel('Q, rate', color=color1, labelpad=15)
    ax1.plot(months, q_monthly, color=color1)
    ax1.tick_params(axis='y', labelcolor=color1, rotation=45)

    color2 = 'tab:gray'
    ax2 = ax1.twinx()
    ax2.set_ylabel('Q, cumulative', color=color2, labelpad=15)
    ax2.plot(months, monthly_cum, '--', color=color2)
    ax2.tick_params(axis='y', labelcolor=color2, rotation=45)

    ax1.get_yaxis().set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%1.0e'))

    plt.title('Arps Forecast')
    ax1.grid(True, linestyle='dotted')
    plt.tight_layout()
    plt.show()

    return

if __name__ == '__main__':

    ### INPUTS ###
    Qi = 10000
    D = .7
    B = 1.4
    dlim = .1
    years = 40
    months = years * 12
    ### INPUTS ###

    q_daily, q_monthly = Arps(Qi, D, B, dlim, years)
    q_cum = cumArps(Qi, D, B, dlim, months)

    visualize(q_monthly)
