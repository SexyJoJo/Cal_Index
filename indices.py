import math
from itertools import groupby

import numpy as np
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
from sympy import *
import pandas as pd
import metpy.calc as mpcalc
from metpy.units import units

microwave_height = (0, 10, 25, 50, 75, 100, 130, 160, 190, 220, 250, 280, 310, 340, 370, 400, 430, 460, 490, 520,
                    560, 600, 640, 680, 720, 760, 800, 840, 880, 920, 960, 1000, 1040, 1080, 1120, 1160, 1200, 1260,
                    1320, 1380, 1440, 1500, 1560, 1620, 1680, 1740, 1800, 1890, 1980, 2170, 2260, 2350, 2430, 2500,
                    2600,
                    2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3650, 3800, 3950, 4100, 4250, 4400, 4550,
                    4600,
                    4800, 5000, 5200, 5400, 5600, 5800, 6000, 6300, 6600, 6900, 7200, 7500, 7800, 8100, 8400, 8700,
                    9000,
                    9300, 9600, 9800, 10000)
K = 0.286
C_pd = 1004.675


def cal_dewp(data_temp, data_rh):
    """
    给定温度和相对湿度，求出露点温度（列表或单个值）
    :param data_temp: 温度数据集
    :param data_rh: 相对湿度数据集
    :return: 露点温度
    """
    if isinstance(data_temp, list) and isinstance(data_rh, list):
        data_dewp = []
        list_zip = zip(data_temp, data_rh)
        for t, f in list_zip:
            x = 1 - 0.01 * f
            dpd = (14.55 + 0.114 * t) * x + ((2.5 + 0.007 * t) * x) ** 3 + (15.9 + 0.117 * t) * (x ** 14)
            Td = t - dpd
            data_dewp.append(round(Td, 3))
        return data_dewp
    else:
        x = 1 - 0.01 * data_rh
        dpd = (14.55 + 0.114 * data_temp) * x + ((2.5 + 0.007 * data_temp) * x) ** 3 + (15.9 + 0.117 * data_temp) * (
                x ** 14)
        Td = data_temp - dpd
        return round(Td, 3)


def CAPE_index(prs, tem, dewp):
    try:
        tem = get_tv(tem, prs)
        p = pd.Series(prs).values * units.hPa
        T = pd.Series(tem).values * units.degC
        Td = pd.Series(dewp).values * units.degC
        prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')
        cape, cin = mpcalc.cape_cin(p, T, Td, prof, which_el='most_cape')
        return float(str(cape).split()[0])
    except Exception as e:
        print(e)
        return None


def LFC_index(prs, tem, dewp):
    """计算自由对流高度（CAPE的最低点）"""
    try:
        tem = get_tv(tem, prs)
        p = pd.Series(prs).values * units.hPa
        T = pd.Series(tem).values * units.degC
        Td = pd.Series(dewp).values * units.degC
        LFC = mpcalc.lfc(p, T, Td)
        LFC = round(float(str(LFC[0]).split()[0]), 3)
        return LFC
    except Exception:
        return None


def get_tv(tems, prses):
    """计算虚温"""
    tvs = []
    for i in range(len(tems)):
        es = get_es(tems[i])
        q = get_r(es, prses[i])  # 比湿
        tv = tems[i] + 0.61 * q * tems[i]
        tvs.append(tv)
    return tvs


def get_tc(t0, td0):
    """
    给定初始的温度和压强求出抬升凝结温度tc
    :param t0: 初始温度
    :param td0: 初始露点温度
    :return:抬升凝结温度tc
    """
    tc = t0 - ((t0 - td0) / (0.976 - 0.000833 * ((237.3 + td0) ** 2) / (273.15 + td0))) * 0.976
    return round(tc, 3)


def get_pc(t0, p0, tc):
    """
    给定初始的温度，压强，抬升凝结温度求出抬升凝结高度压强pc
    :param t0: 初始温度
    :param p0: 初始压强
    :param tc: 抬升凝结温度
    :return: 抬升凝结高度压强pc
    """
    pc = p0 * ((273.15 + tc) / (273.15 + t0)) ** (1 / K)
    return round(pc, 3)


def T(t):
    """
    将摄氏温度转换为开氏温度
    """
    return t + 273.15


def P_d(P, t):
    """
    获得系统干空气分压
    :param P: 气压
    :param t: 气温
    """
    return P - get_es(t)


def r_s(P, t):
    """
    获得饱和混合比(g/g)
    :param P:气压
    :param t:气温(℃)
    :return:饱和混合比
    """
    # 622 * e_s(t) / P    单位为(g/kg)
    return 0.622 * get_es(t) / P


def get_es(td):
    """
    求饱和水汽压，分冰面和水面,根据温度(°C)得到饱和水汽压(mb) 如果用露点代替温度则可得到水汽压 （天气分析 P5）
    param td: 露点温度 单位℃
    return es: 水汽压 单位hPa
    """
    # if t >= 0:     # 水面饱和水汽压的计算公式
    #     return 6.11 * 10 ** ((7.5 * t) / (237.3 + t))
    # else:        # 冰面饱和水汽压的计算公式
    #     return 6.11 * 10 ** ((9.5 * t) / (265.5 + t))
    # es = 6.112 * math.exp((17.67 * t) / (t + 243.5))  # 修正的公式

    # 焦老师提供的方法
    # 作为水面处理
    if td >= -15:
        a = 17.2693882
        b = 35.86
    # 作为冰面处理
    elif td <= -40:
        a = 21.8745584
        b = 7.66
    # 作为水面共存处理
    else:
        fa = interp1d([-15, -40], [17.2693882, 21.8745584])
        fb = interp1d([-15, -40], [35.86, 7.66])
        a = fa(td)
        b = fb(td)
    es = 6.1078 * math.exp((a * td) / (273.16 + td - b))
    return es


def L_v(t):
    return 2.5 * 10 ** 6 - 2323 * t  # 大气物理学P129-130


def T_se(P, t):
    """
    :param t: 温度
    :param T: 绝对温度
    :param P_d: 系统干空气分压
    :param r_s: 饱和混合比
    :param L_v: 汽化潜热
    :return: 假相当位温
    """
    return T(t) * ((1000 / P_d(P, t)) ** K) * math.exp((r_s(P, t) * L_v(t)) / (C_pd * T(t)))


def T_1(P1, t0, T_se_start):
    """二分法公式，给定3参数可迭代出下一个T"""
    T0 = T(t0)
    return T0 - ((T_se(P1, t0) - T_se_start) / (T_se(P1, t0) - T_se(P1, t0 - 0.1))) * 0.1 - 273.15


def getthe(p, t, td, q):
    if (td - t) >= -0.1:
        tlcl = t
    else:
        tlcl = 56.0 + ((td - 56.0) ** (-1) + 0.00125 * log(t / td)) ** (-1)
    gettheresult = t * ((100000.0 / p) ** (0.2854 * (1.0 - 0.28 * q))) * math.exp(
        ((3376.0 / tlcl) - 2.54) * q * (1.0 + 0.81 * q))

    # csnil = dbstack(1);
    # csnil = csnil(1).name(1)
    # ~ = '@';
    # if csnil & & ~isempty(
    #         inputname(3)), assignin('caller', 'FUntemp', td); evalin('caller', [inputname(3), '=FUntemp;']); end
    # if csnil & & ~isempty(
    #         inputname(2)), assignin('caller', 'FUntemp', t); evalin('caller', [inputname(2), '=FUntemp;']); end
    # if csnil & & ~isempty(
    #         inputname(4)), assignin('caller', 'FUntemp', q); evalin('caller', [inputname(4), '=FUntemp;']); end
    # if csnil & & ~isempty(
    #     inputname(1)), assignin('caller', 'FUntemp', p); evalin('caller', [inputname(1), '=FUntemp;']); end
    return gettheresult, p, t, td, q


def get_temp(press, data_press, data_temp):
    f = interp1d(data_press, data_temp)
    temp = f(press)
    return temp


def get_dewp(press, data_press, data_dewp):
    f = interp1d(data_press, data_dewp)
    dewp = f(press)
    return dewp


def get_speed(press, data_press, data_speed):
    f = interp1d(data_press, data_speed)
    speed = f(press)
    return speed


def get_direct(press, data_press, data_direct):
    f = interp1d(data_press, data_direct)
    direct = f(press)
    if direct > 360:
        return None
    return direct


def get_TT(data_press, data_temp, data_dewp):
    """
    计算全总指数(TT)
    :param data_press: 气压数据(列表)
    :param data_temp: 温度数据(列表)
    :param data_dewp: 露点数据(列表)
    """
    T500 = get_temp(500, data_press, data_temp)
    T850 = get_temp(850, data_press, data_temp)
    Td850 = get_dewp(850, data_press, data_dewp)
    TT = T850 + Td850 - 2 * T500
    return round(TT, 3)


def TT_Index(data_press, data_temp, data_dewp):
    """全总指数"""
    try:
        data = get_TT(data_press, data_temp, data_dewp)
        return data
    except Exception:
        return None


def get_r(e, p):
    """
    获得水汽混合比 （g/g）  约等于比湿q      （天气分析 P9）
    :param e:水汽压
    :param p:气压
    """
    return 0.622 * e / (p - e)


def get_rh(t, dewp, pres):
    """
    给定温度、露点、气压获得相对湿度
    :param t: 温度（℃）
    :param dewp: 露点温度（℃）
    :param pres: 气压（mb）
    """
    if t == 9999 or dewp == 9999 or pres == 9999:
        return 9999

    es = get_es(t)  # 根据温度(°C)得到饱和水汽压(mb) 如果用露点代替温度则可得到水汽压 （天气分析 P5）
    if t > 0:
        A = 2.35 * (10 ** 9)
        B = 5420
    else:
        A = 3.41 * (10 ** 10)
        B = 6130
    w = A * 0.622 * math.exp(-B / T(dewp)) / pres
    ws = get_r(es, pres)  # 获得水汽混合比 （g/g）
    return 100 * w / ws


def get_tw(t, Td, pres):
    """
    获得湿球温度（℃）
    :param t: 气温 ℃
    :param Td: 露点温度  ℃
    :param pres: 气压  hpa
    """
    # 王师兄的方法
    # E = get_es(t)  # 饱和水汽压
    #
    # q = get_r((E * get_rh(t, Td, pres)) / 100, pres)  # 比湿
    # h = 1.01 * t + (2500 + 1.84 * t) * q  # 干空气比晗
    #
    # Tw = -100
    # hw = -100
    # while abs(hw - h) > 0.1:
    #     Ew = get_es(Tw)
    #     qw = get_r(Ew, pres)
    #     hw = 1.01 * Tw + (2500 + 1.84 * Tw) * qw
    #     Tw = Tw + 0.01
    # return Tw

    # 直接公式计算法
    # rh = (math.e ** ((17.625 * Td) / (243.104 + Td))) / (math.e ** ((17.625 * t) / (243.104 + t))) * 100
    # Tw = t * math.atan(0.151977 * (rh + 8.313659) ** (1 / 2)) + math.atan(t + rh) - \
    #      math.atan(rh - 1.676331) + 0.00391838 * (rh ** (3 / 2)) * math.atan(0.023101 * rh) - 4.686035
    # return Tw

    # 焦老师提供的方法
    T = t + 273.15
    L = get_L(t)
    r = get_r(get_es(t), pres)
    Tw = 200
    while Tw < 373.15:
        es = get_es(Tw - 273.15)
        right = T - L * (0.622 * es / (pres - es) - r) / 1004.675
        if abs(Tw - right) < 0.3:
            return Tw - 273.15
        else:
            Tw += 0.1
    else:
        return None


def get_theta(t, P):
    """
    获得未饱和湿位温(°C)  位温定义为空气沿干绝热线过程变化到气压p=1000hPa时的温度（天气分析 P16）
    :param t: 气温(°C)
    :param P: 气压(mb)
    """
    theta = (t + 273.15) * ((1000 / P) ** 0.286) - 273.15  # (°C) C#的公式减273.15
    # theta = T(t) * ((1000 / P) ** 0.286)    # (K)书上的公式
    return theta


def JI_Index(data_press, data_temp, data_dewp):
    """
    计算杰弗森指数（JI）
    杰弗森指数, 计算公式为JI=1.6*WBPT850-T500-0.5*(T700-Td700)-8, WBPT为湿球潜温
    :param data_press:压强数据
    :param data_temp:温度数据
    :param data_dewp:露点数据
    """
    try:
        T500 = get_temp(500, data_press, data_temp)
        T700 = get_temp(700, data_press, data_temp)
        Td700 = get_dewp(700, data_press, data_dewp)
        # T850 = get_temp(850, data_press, data_temp)
        # Td850 = get_dewp(850, data_press, data_dewp)
        # thetaw850 = get_theta(get_tw(T850, Td850, 850), 850)

        T900 = get_temp(900, data_press, data_temp)
        Td900 = get_dewp(900, data_press, data_dewp)
        # thetaw900 = get_theta(get_tw(T900, Td900, 900), 900)
        Tw900 = get_tw(T900, Td900, 900)
        Thetaw900 = Tw900 * ((1000 / 900) ** 0.286)
        JI = 1.6 * Thetaw900 - T500 - 0.5 * (T700 - Td700) - 8
        # JI = 1.6 * Thetaw900 - T500 - 11
        return round(JI, 3)
    except Exception:
        return None


def A_index(data_press, data_temp, data_dewp):
    """
     给定高度数据集和温度数据集以及相对湿度数据集，计算A指数
    :param data_press: 气压数据集
    :param data_temp: 温度数据集
    :param data_dewp: 露点数据集
    """
    try:
        T500 = get_temp(500, data_press, data_temp)
        Td500 = get_dewp(500, data_press, data_dewp)
        T850 = get_temp(850, data_press, data_temp)
        Td850 = get_dewp(850, data_press, data_dewp)
        T700 = get_temp(700, data_press, data_temp)
        Td700 = get_dewp(700, data_press, data_dewp)
        A = (T850 - T500) - ((T850 - Td850) + (T700 - Td700) + (T500 - Td500))
        return round(A, 3)
    except Exception:
        return None


def K_index(data_press, data_temp, data_dewp):
    """
    给定高度数据集和温度数据集以及相对湿度数据集，计算K指数
    :param data_press: 气压数据集
    :param data_temp: 温度数据集
    :param data_dewp: 露点数据集
    :return:K指数
    """
    try:
        T500 = get_temp(500, data_press, data_temp)
        T700 = get_temp(700, data_press, data_temp)
        Td700 = get_dewp(700, data_press, data_dewp)
        T850 = get_temp(850, data_press, data_temp)
        Td850 = get_dewp(850, data_press, data_dewp)  # 850压强下的露点温度
        K = T850 - T500 + Td850 - T700 + Td700
        return round(K, 3)
    except Exception:
        return None


def SI_index(data_press, data_temp, data_dewp):
    """
    给定高度数据集和温度数据集以及相对湿度数据集，求出沙瓦特指数
    :param data_press: 气压数据集
    :param data_temp: 温度数据集
    :param data_dewp: 露点数据集
    :return: 沙瓦特指数
    """
    try:
        T850 = get_temp(850, data_press, data_temp)
        Td850 = get_dewp(850, data_press, data_dewp)
        tc = get_tc(T850, Td850)  # 抬升凝结温度
        pc = get_pc(T850, 850, tc)  # 抬升凝结压强
        T500 = get_temp(500, data_press, data_temp)
        T_se_start = T_se(pc, tc)  # 抬升凝结高度处对应的初始假相当位温
        while pc > 500:
            pc = pc - 1
            ret_T_1 = T_1(pc, tc, T_se_start)
            tc = ret_T_1
            if pc <= 500:
                return round(T500 - tc, 3)
        else:
            return None
    except Exception:
        return None


def get_press(data_hei):
    """
    给定高度求出压强
    :param data_hei: 高度(列表或者单个值)
    :return: 压强
    """
    if isinstance(data_hei, list):
        data_press = []
        for hei in data_hei:
            press = 1013.25 * (1 - hei / 44331) ** (10000 / 1903)
            data_press.append(round(press, 3))
        return data_press
    else:
        press = 1013.25 * (1 - data_hei / 44331) ** (10000 / 1903)
        return round(press, 3)


def LI_index(data_press, data_temp, data_dewp):
    """
    给定高度数据集和温度数据集以及相对湿度数据集，求出抬升指数 单位：℃
    :param data_press: 气压数据集
    :param data_temp: 温度数据集
    :param data_dewp: 露点数据集
    2019.7.6(此处的假相当位温仍就是未订正前的)
    """
    try:
        if data_press[0] > 1200:
            return None

        P900 = height2prs(900)  # 得到900m高度的气压
        T = get_temp(P900, data_press, data_temp)
        Td = get_dewp(P900, data_press, data_dewp)

        tc = get_tc(T, Td)  # 抬升凝结温度
        pc = get_pc(T, P900, tc)  # 抬升凝结压强
        T500 = get_temp(500, data_press, data_temp)
        T_se_start = T_se(pc, tc)
        while pc > 500:
            pc = pc - 1
            ret_T_1 = T_1(pc, tc, T_se_start)
            tc = ret_T_1
            if pc <= 500:
                return round(T500 - tc, 3)
        else:
            return None
    except Exception:
        return None


def S_index(data_press, data_temp, data_dewp):
    """
    计算S指数，用于表明雷暴发生潜力指数，常用于4月份到9月份雷暴发生可能性预测；
    :param data_press: 气压数据集
    :param data_temp: 温度数据集
    :param data_dewp: 露点数据集
    """
    T500 = get_temp(500, data_press, data_temp)
    T700 = get_temp(700, data_press, data_temp)
    T850 = get_temp(850, data_press, data_temp)
    Td700 = get_dewp(700, data_press, data_dewp)
    TT = get_TT(data_press, data_temp, data_dewp)
    VT = T850 - T500
    if VT >= 25:
        K = 0
    elif 22 < VT < 25:
        K = 2
    else:
        K = 6
    S = TT - (T700 - Td700) - K
    return round(S, 3)


def TQ_index(data_press, data_temp, data_dewp):
    """
    计算TQ指数，用于评估底层潜在对流情况，
    :param data_press: 气压数据集
    :param data_temp: 温度数据集
    :param data_dewp: 露点数据集
    """
    T700 = get_temp(700, data_press, data_temp)
    T850 = get_temp(850, data_press, data_temp)
    Td850 = get_dewp(850, data_press, data_dewp)
    TQ = (T850 + Td850) - 1.7 * T700
    return round(TQ, 3)


def CT_index(data_press, data_temp, data_dewp):
    """
    计算CT交叉总指数，用于重要天气过程指数
    :param data_press: 气压数据集
    :param data_temp: 温度数据集
    :param data_dewp: 露点数据集
    """
    T500 = get_temp(500, data_press, data_temp)
    Td850 = get_dewp(850, data_press, data_dewp)
    CT = Td850 - T500
    return round(CT, 3)


def VT_index(data_press, data_temp):
    """
    计算VT垂直总指数，通过温度和湿度来预测重要天气过程
    :param data_press: 气压数据集
    :param data_temp: 温度数据集
    """
    T500 = get_temp(500, data_press, data_temp)
    T850 = get_temp(850, data_press, data_temp)
    VT = T850 - T500
    return round(VT, 3)


def get_L(t):
    """
    给定温度得到潜热L（J/kg)
    :param tempc:温度（℃）
    """
    L = (2500 - 2.4 * t) * 1000
    return L


def ThetaSe(t, td, P):  # 2019.05.30 修改后的假相当位温计算函数
    """
    :return: 假相当位温
    """
    e = get_es(td)  # 水汽压
    r = get_r(e, P)  # 起点混合比
    theta = get_theta(t, P) + 273.15
    tc = get_tc(t, td) + 273.15

    # 计算抬升凝结混合比
    # pc = get_pc(P, t, tc)
    # tdc = get_dewp(pc, data_press, data_dewp)   # 抬升凝结露点
    # ec = get_es(tdc)
    # rc = get_r(ec, P)

    Lw = get_L(tc)
    # Lw = 2.501 * pow(10, -6)
    thetase = theta * exp((r * Lw) / (1004.675 * tc)) * (1 + 0.46 * r)
    return round(thetase, 3)


def KO_index(data_press, data_temp, data_dewp):
    """
    计算KO指数，用于评估雷暴发生指数       # 根据C#程序翻译的
    :param data_press: 压强列表
    :param data_temp: 温度列表
    :param data_dewp: 露点列表
    """
    T1000 = get_temp(1000, data_press, data_temp)
    Td1000 = get_dewp(1000, data_press, data_dewp)
    T850 = get_temp(850, data_press, data_temp)
    Td850 = get_dewp(850, data_press, data_dewp)
    T700 = get_temp(700, data_press, data_temp)
    Td700 = get_dewp(700, data_press, data_dewp)
    T500 = get_temp(500, data_press, data_temp)
    Td500 = get_dewp(500, data_press, data_dewp)
    thetae1000 = ThetaSe(T1000, Td1000, 1000)
    thetae850 = ThetaSe(T850, Td850, 850)
    thetae700 = ThetaSe(T700, Td700, 700)
    thetae500 = ThetaSe(T500, Td500, 500)
    KOI = 0.5 * (thetae500 + thetae700) - 0.5 * (thetae850 + thetae1000)
    return round(KOI, 3)


def prs2height(prs):
    return round(44331 * (1 - (prs / 1013.25) ** 0.1903), 3)

def bat_prs2height(pressures):
    """
    气压列表批量转化为高度列表
    @param pressures: 气压列表 单位：hPa
    @return: 高度列表 单位：m
    """
    heights = []
    for prs in pressures:
        heights.append(prs2height(prs))
    return heights

def height2prs(height):
    return round(pow(1 - height / 44331, 1.0 / 0.1903) * 1013.25, 3)

def bat_height2prs(heights):
    """
    高度列表批量转化为气压列表
    @param heights: 高度列表 单位：m
    @return: 气压列表 单位：hPa
    """
    pressures = []
    for height in heights:
        pressures.append(height2prs(height))
    return pressures

# def LCL_index(data_press, data_temp, data_dewp):
#     """
#     :param data_press: 压强列表
#     :param data_temp: 温度列表
#     :param data_dewp: 露点列表
#     """
#     # temp_0 = data_temp[0]
#     # dewp_0 = data_dewp[0]
#     # old = 123 * (temp_0 - dewp_0)
#     # print('old LCL', old)
#
#     T1000 = get_temp(1000, data_press, data_temp)
#     Td1000 = get_dewp(1000, data_press, data_dewp)
#     lcl = 123 * (T1000 - Td1000)
#     return height2prs(lcl)

def LCL_index(data_press, data_temp, data_dewp):
    tc = get_tc(data_temp[0], data_dewp[0])
    return get_pc(data_temp[0], data_press[0], tc)


def TI_index(data_press, data_temp, data_dewp):
    """
    给定高度数据集和温度数据集以及相对湿度数据集，计算TI指数(Thompson指数)
    :param data_press: 微波辐射计高度数据集
    :param data_temp: 微波辐射计温度数据集
    :param data_dewp: 微波辐射计相对湿度数据集
    """
    KI = K_index(data_press, data_temp, data_dewp)
    LI = LI_index(data_press, data_temp, data_dewp)
    TI = KI - LI
    return round(TI, 3)


def IC_index(data_press, data_temp, data_dewp):
    """
    计算条件性稳定度指数      计算公式中是假相当位温的标志，C#程序中计算的是相当位温，此处先用相当位温计算的

    2019.7.6（王磊） 使用订正后的假相当位温计算

    :param data_press: 微波辐射计高度数据集
    :param data_temp: 微波辐射计温度数据集
    :param data_dewp: 微波辐射计相对湿度数据集
    """
    try:
        T850 = get_temp(850, data_press, data_temp)
        Td850 = get_dewp(850, data_press, data_dewp)
        T500 = get_temp(500, data_press, data_temp)
        Td500 = get_dewp(500, data_press, data_dewp)
        Thetase850 = ThetaSe(T850, Td850, 850)
        Thetase500 = ThetaSe(T500, Td500, 500)
        Icondi = Thetase500 - Thetase850
        return round(Icondi, 3)
    except Exception:
        return None


def BIC_index(data_press, data_temp, data_dewp):
    """
    最大对流稳定度指数
    :param data_press: 微波辐射计高度数据集
    :param data_temp: 微波辐射计温度数据集
    :param data_dewp: 微波辐射计相对湿度数据集
    """
    Thetase_list = []
    prs = data_press[0]
    while prs >= 200:
        T_prs = get_temp(prs, data_press, data_temp)
        Td_prs = get_dewp(prs, data_press, data_dewp)
        Thetase = ThetaSe(T_prs, Td_prs, prs)
        prs -= 20
        Thetase_list.append(Thetase)

    min_Thetase = min(Thetase_list)
    max_Thetase = max(Thetase_list)
    Icondi = max_Thetase - min_Thetase
    return round(Icondi, 3)


def DCI_index(data_press, data_temp, data_dewp):
    """
    计算深对流指数
    :param data_press: 微波辐射计高度数据集
    :param data_temp: 微波辐射计温度数据集
    :param data_dewp: 微波辐射计相对湿度数据集
    """
    try:
        T850 = get_temp(850, data_press, data_temp)
        Td850 = get_dewp(850, data_press, data_dewp)
        LI = LI_index(data_press, data_temp, data_dewp)
        DCI = T850 + Td850 - LI
        return round(DCI, 3)
    except Exception:
        return None


def SWEAT_index(data_press, data_temp, data_dewp, data_speed, data_direct):
    """
    计算强天气威胁指数，乔辰炜实现
    :param data_press: 微波辐射计高度数据集
    :param data_temp: 微波辐射计温度数据集
    :param data_dewp: 微波辐射计相对湿度数据集
    :param data_speed: 微波辐射计风速数据集 单位：海里/小时
    :param data_direct: 微波辐射计风向数据集
    """
    try:
        Td850 = get_dewp(850, data_press, data_dewp)
        # 海里转换m/s，乘2
        V850 = get_speed(850, data_press, data_speed) * 2
        V500 = get_speed(500, data_press, data_speed) * 2
        if V500 > 9000 or V850 > 9000:
            return None
        dd850 = get_direct(850, data_press, data_direct)
        dd500 = get_direct(500, data_press, data_direct)
        if dd850 > 9000 or dd500 > 9000:
            return None

        TT = TT_Index(data_press, data_temp, data_dewp)
        # if TT < 49:
        #     SWEAT = 12 * Td850 + 2 * V850 + V500 + 125 * (math.sin(dd500 - dd850) + 0.2)
        # else:
        #     SWEAT = 12 * Td850 + 20 * (TT - 49) + 2 * V850 + V500 + 125 * (math.sin(dd500 - dd850) + 0.2)
        if Td850 < 0:
            Td850 = 0
        if TT < 49:
            TT = 49
        last = 125 * (math.sin(dd500 - dd850) + 0.2)
        if dd500 < 210 or dd500 > 310:
            last = 0
        if dd850 < 130 or dd850 > 250:
            last = 0
        if dd500 <= dd850 or V500 < 15 or V850 < 15:
            last = 0
        SWEAT = 12 * Td850 + 20 * (TT - 49) + 2 * V850 + V500 + last
        return SWEAT
    except Exception:
        return None


# def SSI_index(data_press, data_temp, data_dewp, data_height, data_wspeed):
#     """计算风暴强度指数"""
#     if len(data_height) != len(data_wspeed):
#         raise Exception("高度与风速的数组长度不一致")
#     CAPE = CAPE_index(data_press, data_temp, data_dewp)
#     f = interp1d(data_height, data_wspeed)
#     # TODO
#
#
def KYI_index(data_press, data_temp, data_dewp, data_wspeed, data_wdirect, lat):
    """
    计算山崎指数
    :param data_press: 微波辐射计高度数据集
    :param data_temp: 微波辐射计温度数据集
    :param data_dewp: 微波辐射计相对湿度数据集
    :param data_wspeed: 风速数据集
    :param data_wdirect: 风向数据集
    :param lat: unknown
    """
    try:
        SI = SI_index(data_press, data_temp, data_dewp)
        if SI is None:
            return None

        # 计算温度平流TA
        # ((0.601*power(10,-4))*((850_s*500_s) + power(((850_s-500_s)/3),2)) * (850_d-500_d) * sin(lat/180*ACOS(-1)))/3600 * 100000
        w_speed850 = get_speed(850, data_press, data_wspeed)
        w_speed500 = get_speed(500, data_press, data_wspeed)
        w_direct850 = get_direct(850, data_press, data_wdirect)
        w_direct500 = get_direct(500, data_press, data_wdirect)
        TA = (0.601 * pow(10, -4)) * ((w_speed850 * w_speed500) + pow((w_speed850 - w_speed500) / 3, 2) * (
                w_direct850 - w_direct500) * math.sin(lat/180 * math.acos(-1))) / 3600 * 100000

        T850 = get_temp(850, data_press, data_temp)
        Td850 = get_dewp(850, data_press, data_temp)
        if TA <= SI:
            KYI = 0
        else:
            KYI = (TA - SI) / (1 + (T850 - Td850))
        return round(KYI, 3)
    except Exception:
        return None


def inver_height(height, data_temp):
    """逆温层"""
    all_point = []  # 出现逆温层的索引列表
    for index, val in enumerate(data_temp):
        if index == len(data_temp) - 1:
            break
        if val <= data_temp[index + 1]:
            all_point.append(index)
            # if val==data[index+1]:
            #     all_point.append(index+1)

    con_digits = continuous_digits(all_point)
    delta_h, delta_t, h0 = 0, 0, None
    for i in con_digits:
        diff_val = round(data_temp[max(i)] - data_temp[min(i)], 2)
        if diff_val > 1 and height[min(i)] < 1000:
            delta_t += diff_val
            delta_h += (height[max(i)] - height[min(i)])
            if h0 is None:
                h0 = height[min(i)]

    return {
        "delta_h": delta_h,
        "delta_t": delta_t,
        "h0": h0
    }


def continuous_digits(lst):
    list_all = []
    fun = lambda x: x[1] - x[0]
    for k, g in groupby(enumerate(lst), fun):
        l1 = [j for i, j in g]  # 连续数字的列表
        if len(l1) > 1:
            list_all.append(l1)
    return list_all
