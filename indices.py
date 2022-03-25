import math

import numpy as np
from numpy import zeros
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
from sympy import *

microwave_height = (0, 10, 25, 50, 75, 100, 130, 160, 190, 220, 250, 280, 310, 340, 370, 400, 430, 460, 490, 520,
                    560, 600, 640, 680, 720, 760, 800, 840, 880, 920, 960, 1000, 1040, 1080, 1120, 1160, 1200, 1260,
                    1320,
                    1380, 1440, 1500, 1560, 1620, 1680, 1740, 1800, 1890, 1980, 2170, 2260, 2350, 2430, 2500, 2600,
                    2700,
                    2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3650, 3800, 3950, 4100, 4250, 4400, 4550, 4600,
                    4800,
                    5000, 5200, 5400, 5600, 5800, 6000, 6300, 6600, 6900, 7200, 7500, 7800, 8100, 8400, 8700, 9000,
                    9300,
                    9600, 9800, 10000)
K = 0.286  # python常量：约定俗成，不可更改，全部是大写字母
C_pd = 1004


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


def get_es(t):
    """
    求饱和水汽压，分冰面和水面,根据温度(°C)得到饱和水汽压(mb) 如果用露点代替温度则可得到水汽压 （天气分析 P5）
    """
    # if t >= 0:     # 水面饱和水汽压的计算公式
    #     return 6.11 * 10 ** ((7.5 * t) / (237.3 + t))
    # else:        # 冰面饱和水汽压的计算公式
    #     return 6.11 * 10 ** ((9.5 * t) / (265.5 + t))
    es = 6.112 * math.exp((17.67 * t) / (t + 243.5))  # 修正的公式
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


def moist_adiabat(pc, tc):
    """
    湿绝热线
    :param tc: 抬升凝结温度
    :param pc: 抬升凝结压强
    :return: 湿绝热线数据集
    """
    T_se_start = T_se(pc, tc)
    T_list = [tc, ]
    P_list = [pc, ]
    while pc > 210:
        t_i = tc
        p_i = pc - 10
        ret_T_1 = T_1(p_i, t_i, T_se_start)
        tc = ret_T_1
        pc = pc - 10
        P_list.append(round(pc, 3))
        T_list.append(round(tc, 3))
    return P_list, T_list


def inter_un(xi, yi, xnew):
    """
    插值，范围外插,(xi必须严格单调递增)
    """
    f = InterpolatedUnivariateSpline(xi, yi, k=1)
    return [round(i, 3) for i in list(f(xnew))]


def insert_moist_adiabat_data(moist_adiabat_P, moist_adiabat_T, press):
    """得到插值后的湿绝热线温度值，与压强一一对应"""
    xi = sorted(moist_adiabat_P)
    yi = sorted(moist_adiabat_T)
    moist_adiabat_data = inter_un(xi, yi, press)
    return moist_adiabat_data


def _rect_inter_inner(x1, x2):
    n1 = x1.shape[0] - 1
    n2 = x2.shape[0] - 1
    X1 = np.c_[x1[:-1], x1[1:]]
    X2 = np.c_[x2[:-1], x2[1:]]
    S1 = np.tile(X1.min(axis=1), (n2, 1)).T
    S2 = np.tile(X2.max(axis=1), (n1, 1))
    S3 = np.tile(X1.max(axis=1), (n2, 1)).T
    S4 = np.tile(X2.min(axis=1), (n1, 1))
    return S1, S2, S3, S4


def _rectangle_intersection_(x1, y1, x2, y2):
    S1, S2, S3, S4 = _rect_inter_inner(x1, x2)
    S5, S6, S7, S8 = _rect_inter_inner(y1, y2)
    C1 = np.less_equal(S1, S2)
    C2 = np.greater_equal(S3, S4)
    C3 = np.less_equal(S5, S6)
    C4 = np.greater_equal(S7, S8)
    ii, jj = np.nonzero(C1 & C2 & C3 & C4)
    return ii, jj


def intersection(x1, y1, x2, y2):
    ii, jj = _rectangle_intersection_(x1, y1, x2, y2)
    n = len(ii)

    dxy1 = np.diff(np.c_[x1, y1], axis=0)
    dxy2 = np.diff(np.c_[x2, y2], axis=0)

    T = np.zeros((4, n))
    AA = np.zeros((4, 4, n))
    AA[0:2, 2, :] = -1
    AA[2:4, 3, :] = -1
    AA[0::2, 0, :] = dxy1[ii, :].T
    AA[1::2, 1, :] = dxy2[jj, :].T

    BB = np.zeros((4, n))
    BB[0, :] = -x1[ii].ravel()
    BB[1, :] = -x2[jj].ravel()
    BB[2, :] = -y1[ii].ravel()
    BB[3, :] = -y2[jj].ravel()

    for i in range(n):
        try:
            T[:, i] = np.linalg.solve(AA[:, :, i], BB[:, i])
        except:
            T[:, i] = np.NaN

    in_range = (T[0, :] >= 0) & (T[1, :] >= 0) & (T[0, :] <= 1) & (T[1, :] <= 1)

    xy0 = T[2:, in_range]
    xy0 = xy0.T
    return xy0[:, 0], xy0[:, 1]


def CAPE(press, temp, mois, LFC_P=None, LFC_T=None, El_P=None, EL_T=None):
    '''得到交点之间的数据集<=LFC_P and >= El_p'''
    p = press[:]  # 进行列表复制，否则使用del会导致原来变量被删除
    t = temp[:]
    m = mois[:]
    if LFC_P == None and LFC_T == None:
        return 0
    elif El_P == None and EL_T == None:
        for i in range(len(p)):
            if p[i] <= LFC_P:
                p.insert(i, LFC_P)
                t.insert(i, LFC_T)
                m.insert(i, LFC_T)
                del p[:i]
                del t[:i]
                del m[:i]
                break
    else:
        for i in range(len(p)):
            if p[i] <= LFC_P:
                p.insert(i, LFC_P)
                t.insert(i, LFC_T)
                m.insert(i, LFC_T)
                del p[:i]
                del t[:i]
                del m[:i]
                break
        for i in range(1, len(p)):
            if p[-i] >= El_P:
                p.insert(-i + 1, El_P)
                t.insert(-i + 1, EL_T)
                m.insert(-i + 1, EL_T)
                del p[-i + 1:]
                del t[-i + 1:]
                del m[-i + 1:]
                break

    '''使用复合梯形规则沿给定轴积分'''
    area1 = np.trapz(m, p)
    # area1 = np.trapz(press, mois)
    area2 = np.trapz(t, p)
    # area2 = np.trapz(press, temp)
    area = abs(area1 - area2)
    return round(area, 3)


def get_LfcAndElAndCape(press, temp, mois):
    """
    计算温度和湿绝热线交点的对应高度LFC和EL的压强
    """
    try:
        index = temp.index(None) + 1
    except:
        index = 0
    f1 = interp1d(press[index:], temp[index:], )
    f2 = interp1d(press[index:], mois[index:], )
    x1 = np.array(temp[index:])
    y1 = np.array(press[index:])
    x2 = np.array(mois[index:])
    intersection_T, intersection_P = intersection(x1, y1, x2, y1)
    data_P = []
    data_T = []
    if len(intersection_T) != 0:  # 判断交点距离，若距离过小，默认交点为列表最后一个交点
        x = max(intersection_T) - min(intersection_T)
        if x <= 1:
            data_P.append(intersection_P[-1])
            data_T.append(intersection_T[-1])
        else:
            data_P = intersection_P
            data_T = intersection_T
    # LFC_T, EL_T = data_T
    # LFC_P, EL_P = data_P
    """
    1、交点数大于等于2
        a、取出交点列表中的最后2个交点（压强值）
        b、判断两交点之间所有的点对应的f2-f1的值是否大于0
        c、大于0返回True,小于0返回False
        d、若f2-f1都大于0，则两交点分别是El和LFC
        e、若f2-f1不是都大于0，则最后一个交点为LFC
    2、交点数等于1
        a、判断该交点到小于其200压强下点之间所有点的对应f2-f1的值是否大于0
        b、大于0返回True,小于0返回False
        c、若都大于0，则点为LFC
        d、若都小于0，则点为EL
    """
    if len(data_P) >= 2:
        first = data_P[-1]
        second = data_P[-2]
        test_data = []
        while first < second - 2:
            first += 2
            if f2(first) - f1(first) >= 0:
                test_data.append(True)
            else:
                test_data.append(False)
        if False in test_data:
            cape = CAPE(press, temp, mois, data_P[-1], data_T[-1])
            height = {
                'LFC': data_P[-1],
                'cape': cape
            }
        else:
            cape = CAPE(press, temp, mois, data_P[-2], data_T[-2], data_P[-1], data_T[-1])
            height = {
                'LFC': data_P[-2],
                "EL": data_P[-1],
                'cape': cape
            }
        return height
    elif len(data_P) == 1:
        first = data_P[0]
        test_data = []
        end = first - 100
        while first >= end:
            first -= 2
            if f2(first) - f1(first) >= 0:
                test_data.append(True)
            else:
                test_data.append(False)
        if False not in test_data:
            cape = CAPE(press, temp, mois, data_P[0], data_T[0])
            height = {
                'LFC': data_P[0],
                'cape': cape
            }
            return height
        if True not in test_data:
            height = {
                'EL': data_P[0],
                'cape': 0
            }
            return height
        else:
            height = {
                'LFC': 9999,
                "EL": 9999,
                'cape': 0
            }
            return height
    else:
        height = {
            'LFC': 9999,
            "EL": 9999,
            'cape': 0
        }
        return height


def cape(press, dewp, t):
    """对流有效位能"""
    result = {
        "unit": 'J/kg',
        "label": "对流有效位能,根据温度列表和相对湿度列表数据，首先计算出湿绝热线，再根据其找出自由对流高度和平衡高度，取两个高度之间的积分和即cape值"
    }
    data = []
    tc = get_tc(t[0], dewp)
    pc = get_pc(t[0], press[0], tc)
    moist_adiabat_P, moist_adiabat_T = moist_adiabat(pc, tc)
    moist_adiabat_data = insert_moist_adiabat_data(moist_adiabat_P, moist_adiabat_T, press)  # 插值后的湿绝热线温度值
    hei_important_data = get_LfcAndElAndCape(press, t, moist_adiabat_data)
    meta = {
        # "x": time,
        "y": hei_important_data.get('cape', 0)
    }
    data.append(meta)
    result["data"] = data
    return data[0]['y']


def getqvs(p, t):
    getqvsresult = []
    es = 0
    eps = 287.04 / 461.5
    es = 611.2 * exp(17.67 * (t - 273.15) / (t - 29.65))
    getqvsresult = eps * es / (p - es)
    # csnil = dbstack(1);
    # csnil = csnil(1).name(1)
    # ~ = '@';
    # if csnil & & ~isempty(
    #         inputname(2)), assignin('caller', 'FUntemp', t); evalin('caller', [inputname(2), '=FUntemp;']); end
    # if csnil & & ~isempty(
    #         inputname(1)), assignin('caller', 'FUntemp', p); evalin('caller', [inputname(1), '=FUntemp;']); end
    return getqvsresult


def getqvi(p, t):
    getqviresult = []
    es = 0
    eps = 287.04 / 461.5
    es = 611.2 * math.exp(21.8745584 * (t - 273.15) / (t - 7.66))
    getqviresult = eps * es / (p - es)
    # csnil=dbstack(1);
    # csnil=csnil(1).name(1)~='@';
    # if csnil & & ~isempty(
    #         inputname(2)), assignin('caller', 'FUntemp', t); evalin('caller', [inputname(2), '=FUntemp;']); end
    # if csnil & & ~isempty(
    #         inputname(1)), assignin('caller', 'FUntemp', p); evalin('caller', [inputname(1), '=FUntemp;']); end
    return getqviresult


def getthe(p, t, td, q):
    tlcl = 0
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


def new_CAPE(p_in, t_in, td_in):
    """
    计算CAPE与CIN（对流抑制）
        p_in - 一维压力数组（mb）（实数）
        t_in - 温度 (C) 的一维数组（实数）
        td_in - 露点温度 (C) 的一维数组（实数）
    """
    source = 1
    nk = len(p_in)
    pinc = 100.0
    # source = 2
    ml_depth = 200.0
    adiabat = 1
    # doit = False
    # ice = False
    # cloud = False
    # not_converged = False
    # k = 0
    kmax = 0
    # n = 0
    # nloop = 0
    # i = 0
    # orec = 0
    p = zeros(1, nk)
    t = zeros(1, nk)
    td = zeros(1, nk)
    pi = zeros(1, nk)
    q = zeros(1, nk)
    th = zeros(1, nk)
    thv = zeros(1, nk)
    z = zeros(1, nk)
    # pt = zeros(1, nk)
    # pb = zeros(1, nk)
    # pc = zeros(1, nk)
    # pn = zeros(1, nk)
    # ptv = zeros(1, nk)
    # the = 0
    # maxthe = 0
    # parea = 0
    # narea = 0
    # lfc = 0
    # th1 = 0
    # p1 = 0
    # t1 = 0
    # qv1 = 0
    # ql1 = 0
    # qi1 = 0
    # b1 = 0
    # pi1 = 0
    # thv1 = 0
    # qt = 0
    # dp = 0
    # dz = 0
    # ps = 0
    # frac = 0
    th2 = 0
    p2 = 0
    t2 = 0
    qv2 = 0
    # ql2 = 0
    # qi2 = 0
    b2 = 0
    pi2 = 0
    thv2 = 0
    # thlast = 0
    # fliq = 0
    # fice = 0
    # tbar = 0
    # qvbar = 0
    # qlbar = 0
    # qibar = 0
    # lhv = 0
    # lhs = 0
    # lhf = 0
    # rm = 0
    # cpm = 0
    avgth = 0
    avgqv = 0
    g = 9.81
    p00 = 100000.0
    cp = 1005.7
    rd = 287.04
    rv = 461.5
    xlv = 2501000.0
    xls = 2836017.0
    t0 = 273.15
    cpv = 1875.0
    cpl = 4190.0
    cpi = 2118.636
    lv1 = xlv + (cpl - cpv) * t0
    lv2 = cpl - cpv
    ls1 = xls + (cpi - cpv) * t0
    ls2 = cpi - cpv
    rp00 = 1.0 / p00
    # eps = rd / rv
    reps = rv / rd
    rddcp = rd / cp
    cpdrd = cp / rd
    cpdg = cp / g
    converge = 0.0002
    debug_level = 0

    # 将 p,t,td 转换为 mks 单位； 得到 pi,q,th,thv
    for k in range(nk):
        p[k] = 100.0 * p_in[k]
        t[k] = 273.15 + t_in[k]
        td[k] = 273.15 + td_in[k]
        pi[k] = (p[k] * rp00) ** rddcp
        q[k] = getqvs(p[k], td[k])
        th[k] = t[k] / pi[k]
        thv[k] = th[k] * (1.0 + reps * q[k]) / (1.0 + q[k])
    k = np.fix(nk + 1)

    # 使用流体静力方程获得高度
    z[0] = 0.0
    for k in range(1, nk):
        dz = -cpdg * 0.5 * (thv[k] + thv[k - 1]) * (pi[k] - pi[k - 1])
        z[k] = z[k - 1] + dz
    k = np.fix(nk + 1)

    if source == 1:  # use surface parcel
        kmax = 1
    elif source == 2:
        if p[0] < 50000.0:
            kmax = 1
            maxthe, p[0], t[0], td[0], q[0] = getthe(p[0], t[0], td[0], q[0])
        else:
            maxthe = 0.0
            for k in range(nk):
                if p[k] >= 50000.0:
                    the, p[k], t[k], td[k], q[k] = getthe(p[k], t[k], td[k], q[k])
                    if the > maxthe:
                        maxthe = the
                        kmax = np.fix(k)
            k = np.fix(nk + 1)

        if debug_level >= 100:
            print('  kmax,maxthe = ', kmax, maxthe)
    elif source == 3:
        if (z[1] - z[0]) > ml_depth:
            avgth = th[1]
            avgqv = q[1]
            kmax = 1
        elif z[nk] < ml_depth:
            avgth = th(nk)
            avgqv = q(nk)
            kmax = np.fix(nk)
        else:
            avgth = 0.0
            avgqv = 0.0
            k = 2
            if debug_level >= 100:
                print('  ml_depth = ', ml_depth)
            if debug_level >= 100:
                print('  k,z,th,q:')
            if debug_level >= 100:
                print(1, z[1], th[1], q[1])
            while (z[k] <= ml_depth) & (k <= nk):
                if debug_level >= 100:
                    print(k, z[k], th[k], q[k])
                avgth = avgth + 0.5 * (z[k] - z[k - 1]) * (th[k] + th[k - 1])
                avgqv = avgqv + 0.5 * (z[k] - z[k - 1]) * (q[k] + q[k - 1])
                k = np.fix(k + 1)

            th2 = th[k - 1] + (th[k] - th[k - 1]) * (ml_depth - z[k - 1]) / (z[k] - z[k - 1])
            qv2 = q[k - 1] + (q[k] - q[k - 1]) * (ml_depth - z[k - 1]) / (z[k] - z[k - 1])
            if debug_level >= 100:
                print(999, ml_depth, th2, qv2)
            avgth = avgth + 0.5 * (ml_depth - z[k - 1]) * (th2 + th[k - 1])
            avgqv = avgqv + 0.5 * (ml_depth - z[k - 1]) * (qv2 + q[k - 1])
            if debug_level >= 100:
                print(k, z(k), th(k), q(k))
            avgth = avgth / ml_depth
            avgqv = avgqv / ml_depth
            kmax = 1
        if debug_level >= 100:
            print(avgth, avgqv)
    else:
        # writef(1, ['%0.15g \n']);
        # writef(1, ['%s \n'], '  Unknown value for source');
        # writef(1, ['%0.15g \n']);
        # writef(1, ['%s %0.15g \n'], '  source = ', source);
        # writef(1, ['%0.15g \n']);
        # error(['stop encountered in original fortran code  ', char(10), ';']);
        print("error")

    narea = 0.0
    if (source == 1) or (source == 2):
        kmax = int(kmax) - 1    # matlab索引从1开始 
        k = np.fix(kmax)
        th2 = th[kmax]
        pi2 = pi[kmax]
        p2 = p[kmax]
        t2 = t[kmax]
        thv2 = thv[kmax]
        qv2 = q[kmax]
        b2 = 0.0
    elif source == 3:
        k = np.fix(kmax)
        th2 = avgth
        qv2 = avgqv
        thv2 = th2 * (1.0 + reps * qv2) / (1.0 + qv2)
        pi2 = pi[kmax]
        p2 = p[kmax]
        t2 = th2 * pi2
        b2 = g * (thv2 - thv[kmax]) / thv[kmax]

    ql2 = 0.0
    qi2 = 0.0
    qt = qv2
    cape = 0.0
    cin = 0.0
    lfc = 0.0
    doit = True
    cloud = False
    if adiabat == 1 or adiabat == 2:
        ice = False
    else:
        ice = True
    t2_orig = t2
    the, p2, t2, dumvar4, qv2 = getthe(p2, t2, t2, qv2)
    # t2[dumvar4 != t2_orig] = dumvar4[dumvar4 != t2_orig]
    if debug_level >= 100:
        print('  the = ', the)
    # if(debug_level>=100)
    # writef(1,['%s \n'], '  Start loop:');
    # writef(1,['%s %0.15g %0.15g %0.15g \n'], '  p2,th2,qv2 = ',p2,th2,qv2);
    # end;

    while doit & (k < nk-1):
        k = int(np.fix(k + 1))
        b1 = b2
        dp = p[k - 1] - p[k]
        if dp < pinc:
            nloop = 1
        else:
            nloop = np.fix(1 + np.fix(dp / pinc))
            dp = dp / nloop
        for n in range(nloop):
            p1 = p2
            t1 = t2
            pi1 = pi2
            th1 = th2
            qv1 = qv2
            ql1 = ql2
            qi1 = qi2
            thv1 = thv2
            p2 = p2 - dp
            pi2 = (p2 * rp00) ** rddcp
            thlast = th1
            i = 0
            not_converged = True
            while not_converged:
                i = np.fix(i + 1)
                t2 = thlast * pi2
                if ice:
                    fliq = max(min((t2 - 233.15) / (273.15 - 233.15), 1.0), 0.0)
                    fice = 1.0 - fliq
                else:
                    fliq = 1.0
                    fice = 0.0
                qv2 = min(qt, fliq * getqvs(p2, t2) + fice * getqvi(p2, t2))
                qi2 = max(fice * (qt - qv2), 0.0)
                ql2 = max(qt - qv2 - qi2, 0.0)
                tbar = 0.5 * (t1 + t2)
                qvbar = 0.5 * (qv1 + qv2)
                qlbar = 0.5 * (ql1 + ql2)
                qibar = 0.5 * (qi1 + qi2)
                lhv = lv1 - lv2 * tbar
                lhs = ls1 - ls2 * tbar
                lhf = lhs - lhv
                rm = rd + rv * qvbar
                cpm = cp + cpv * qvbar + cpl * qlbar + cpi * qibar
                th2 = th1 * math.exp(
                    lhv * (ql2 - ql1) / (cpm * tbar) + lhs * (qi2 - qi1) / (cpm * tbar) + (rm / cpm - rd / cp) * log(
                        p2 / p1))
                # if i > 90:
                #     print(i, th2, thlast, th2 - thlast)
                # if i > 100:
                #     writef(1, ['%0.15g \n'])
                #     writef(1, ['%s \n'], '  Error:  lack of convergence')
                #     writef(1, ['%0.15g \n'])
                #     writef(1, ['%s \n'], '  ... stopping iteration ')
                #     writef(1, ['%0.15g \n'])
                #     error(['stop encountered in original fortran code  ', char(10), ' 1001;'])
                if abs(th2 - thlast) > converge:
                    thlast = thlast + 0.3 * (th2 - thlast)
                else:
                    not_converged = False

            if ql2 >= 1.0e-10:
                cloud = True
            if adiabat == 1 or adiabat == 3:
                qt = qv2
                ql2 = 0.0
                qi2 = 0.0
            # elif adiabat <= 0 or adiabat >= 5:
            # writef(1, ['%0.15g \n']);
            # writef(1, ['%s \n'], '  Undefined adiabat');
            # writef(1, ['%0.15g \n']);
            # error(['stop encountered in original fortran code  ', char(10), ' 10000;']);
        n = np.fix(nloop + 1)
        thv2 = th2 * (1.0 + reps * qv2) / (1.0 + qv2 + ql2 + qi2)
        b2 = g * (thv2 - thv[k]) / thv[k]
        dz = -cpdg * 0.5 * (thv[k] + thv[k - 1]) * (pi[k] - pi[k - 1])
        t2_orig = t2
        the, p2, t2, dumvar4, qv2 = getthe(p2, t2, t2, qv2)
        # t2[dumvar4 != t2_orig] = dumvar4[dumvar4 != t2_orig]
        if (b2 >= 0.0) and (b1 < 0.0):
            ps = p[k - 1] + (p[k] - p[k - 1]) * (0.0 - b1) / (b2 - b1)
            frac = b2 / (b2 - b1)
            parea = 0.5 * b2 * dz * frac
            narea = narea - 0.5 * b1 * dz * (1.0 - frac)
            # if (debug_level >= 200)
            #     writef(1, ['%s %0.15g %0.15g \n'], '      b1,b2 = ', b1, b2);
            # writef(1, ['%s %0.15g %0.15g %0.15g \n'], '      p1,ps,p2 = ', p(k - 1), ps, p(k));
            # writef(1, ['%s %0.15g \n'], '      frac = ', frac);
            # writef(1, ['%s %0.15g \n'], '      parea = ', parea);
            # writef(1, ['%s %0.15g \n'], '      narea = ', narea);
            cin = cin + narea
            narea = 0.0
        elif (b2 < 0.0) and (b1 > 0.0):
            ps = p[k - 1] + (p[k] - p[k - 1]) * (0.0 - b1) / (b2 - b1)
            frac = b1 / (b1 - b2)
            parea = 0.5 * b1 * dz * frac
            narea = -0.5 * b2 * dz * (1.0 - frac)
            # if debug_level >= 200:
            #     writef(1, ['%s %0.15g %0.15g \n'], '      b1,b2 = ', b1, b2)
            #     writef(1, ['%s %0.15g %0.15g %0.15g \n'], '      p1,ps,p2 = ', p(k - 1), ps, p(k))
            #     writef(1, ['%s %0.15g \n'], '      frac = ', frac)
            #     writef(1, ['%s %0.15g \n'], '      parea = ', parea)
            #     writef(1, ['%s %0.15g \n'], '      narea = ', narea)
        elif b2 < 0.0:
            parea = 0.0
            narea = narea - 0.5 * dz * (b1 + b2)
        else:
            parea = 0.5 * dz * (b1 + b2)
            narea = 0.0
        cape = cape + max(0.0, parea)
        # if (debug_level >= 200)
        #     writef(1, [repmat(['%13.4f'], 1, 5), repmat(' ', 1, 2), '%1f' ' \n'], p2, b1, b2, cape, cin, cloud);
        # % format(5(f13
        # .4), 2
        # x, l1);
        # end;
        if (p[k] <= 10000.0) and (b2 < 0.0):
            doit = False
    return round(cape, 3)


def get_temp(press, data_press, data_temp):
    f = interp1d(data_press, data_temp, bounds_error=False, fill_value='extrapolate')
    temp = f(press)
    return temp


def get_dewp(press, data_press, data_dewp):
    f = interp1d(data_press, data_dewp, bounds_error=False, fill_value='extrapolate')
    dewp = f(press)
    return dewp


def get_TT(data_press, data_temp, data_dewp):
    """
    计算全总指数(TT)
    :param data_press: 气压数据(列表)
    :param data_temp: 温度数据(列表)
    :param data_dewp: 露点数据(列表)
    """
    try:
        T500 = get_temp(500, data_press, data_temp)
        T850 = get_temp(850, data_press, data_temp)
        Td850 = get_dewp(850, data_press, data_dewp)
        TT = T850 + Td850 - 2 * T500
        return round(TT, 3)
    except ValueError:
        return 'error'


def TT_Index(data_press, data_temp, data_dewp):
    """全总指数"""
    # height = list(microwave_height)
    # result = {
    #     "unit": "℃",
    #     "label": "全总指数，根据温度列表和相对湿度列表数据计算T850(850压强下温度),Td850（850压强下露点温度），T500（500压强下温度），"
    #              "计算公式为TT=T850+Td850-2*T500，常用于重要天气过程的指示，TT越大越容易发生强对流天气。"
    # }
    # list_data = []
    # for i in range(len(lv2_11)):
    data = get_TT(data_press, data_temp, data_dewp)
    # meta = {
    #     "x": time,
    #     "y": data
    # }
    # list_data.append(meta)
    # result["data"] = list_data
    return data


def get_r(e, p):
    """
    获得水汽混合比 （g/g）        （天气分析 P9）
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
    :param t: 气温
    :param Td: 露点温度
    :param pres: 气压
    """
    E = get_es(t)  # 饱和水汽压

    q = get_r((E * get_rh(t, Td, pres)) / 100, pres)  # 比湿
    h = 1.01 * t + (2500 + 1.84 * t) * q  # 干空气比晗

    Tw = -100
    hw = -100
    while abs(hw - h) > 0.1:
        Ew = get_es(Tw)
        qw = get_r(Ew, pres)
        hw = 1.01 * Tw + (2500 + 1.84 * Tw) * qw
        Tw = Tw + 0.01

    return Tw


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
        T850 = get_temp(850, data_press, data_temp)
        Td850 = get_dewp(850, data_press, data_dewp)
    except ValueError:
        return 'error'
    thetaw850 = get_theta(get_tw(T850, Td850, 850), 850)
    JI = 1.6 * thetaw850 - T500 - 0.5 * (T700 - Td700) - 8
    # result = {
    #     "label": "杰弗森指数, 计算公式为JI=1.6*WBPT850-T500-0.5*(T700-Td700)-8, WBPT为湿球潜温",
    #     'data': [{'x': time, 'y': round(JI, 3)}]
    # }
    return round(JI, 3)


def A_index(data_press, data_temp, data_dewp):
    """
     给定高度数据集和温度数据集以及相对湿度数据集，计算A指数
    :param data_press: 气压数据集
    :param data_temp: 温度数据集
    :param data_dewp: 露点数据集
    """
    T500 = get_temp(500, data_press, data_temp)
    Td500 = get_dewp(500, data_press, data_dewp)
    T850 = get_temp(850, data_press, data_temp)
    Td850 = get_dewp(850, data_press, data_dewp)
    T700 = get_temp(700, data_press, data_temp)
    Td700 = get_dewp(700, data_press, data_dewp)
    A = (T850 - T500) - ((T850 - Td850) + (T700 - Td700) + (T500 - Td500))
    return round(A, 3)


def K_index(data_press, data_temp, data_dewp):
    """
    给定高度数据集和温度数据集以及相对湿度数据集，计算K指数
    :param data_press: 气压数据集
    :param data_temp: 温度数据集
    :param data_dewp: 露点数据集
    :return:K指数
    """
    T500 = get_temp(500, data_press, data_temp)
    T700 = get_temp(700, data_press, data_temp)
    Td700 = get_dewp(700, data_press, data_dewp)
    T850 = get_temp(850, data_press, data_temp)
    Td850 = get_dewp(850, data_press, data_dewp)  # 850压强下的露点温度
    K = T850 - T500 + Td850 - T700 + Td700
    return round(K, 3)


def SI_index(data_press, data_temp, data_dewp):
    """
    给定高度数据集和温度数据集以及相对湿度数据集，求出沙瓦特指数
    :param data_press: 气压数据集
    :param data_temp: 温度数据集
    :param data_dewp: 露点数据集
    :return: 沙瓦特指数
    """
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
    给定高度数据集和温度数据集以及相对湿度数据集，求出抬升指数
    :param data_press: 气压数据集
    :param data_temp: 温度数据集
    :param data_dewp: 露点数据集
    2019.7.6(此处的假相当位温仍就是未订正前的)
    """
    P900 = get_press(900)  # 得到900m高度的气压
    T = get_temp(P900, data_press, data_temp)
    Td = get_dewp(P900, data_press, data_dewp)
    # except ValueError:
    #     return 'error'
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
    :param t: 温度
    :param T: 绝对温度
    :param P_d: 系统干空气分压
    :param r: 饱和混合比
    :param L_v: 汽化潜热
    :return: 假相当位温
    """
    es = get_es(t)  # 饱和水汽压
    r = get_r(es, P)  # 将饱和水汽压代入得到饱和混合比
    theta = get_theta(t, P) + 273.15
    tc = get_tc(t, td) + 273.15
    Lw = get_L(tc)
    thetase = theta * exp((r * Lw) / (1004 * tc))
    return round(thetase, 3)


def KO_index(data_press, data_temp, data_dewp):
    """
    计算KO指数，用于评估雷暴发生指数       # 根据C#程序翻译的
    :param data_press: 微波辐射计高度数据集
    :param data_temp: 微波辐射计温度数据集
    :param data_dewp: 微波辐射计相对湿度数据集
    """
    try:
        T1000 = get_temp(1000, data_press, data_temp)
        Td1000 = get_dewp(1000, data_press, data_dewp)
        T850 = get_temp(850, data_press, data_temp)
        Td850 = get_dewp(850, data_press, data_dewp)
        T700 = get_temp(700, data_press, data_temp)
        Td700 = get_dewp(700, data_press, data_dewp)
        T500 = get_temp(500, data_press, data_temp)
        Td500 = get_dewp(500, data_press, data_dewp)
    except ValueError:
        return 'error'
    thetae1000 = ThetaSe(T1000, Td1000, 1000)
    thetae850 = ThetaSe(T850, Td850, 850)
    thetae700 = ThetaSe(T700, Td700, 700)
    thetae500 = ThetaSe(T500, Td500, 500)
    KOI = 0.5 * (thetae500 + thetae700) - 0.5 * (thetae850 + thetae1000)
    return round(KOI, 3)


def LCL_index(data_temp, data_dewp):
    temp_0 = data_temp[0]
    dewp_0 = data_dewp[0]
    lcl = 123 * (temp_0 - dewp_0)
    return round(lcl, 3)


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
    T850 = get_temp(850, data_press, data_temp)
    Td850 = get_dewp(850, data_press, data_dewp)
    T500 = get_temp(500, data_press, data_temp)
    Td500 = get_dewp(500, data_press, data_dewp)
    Thetase850 = ThetaSe(T850, Td850, 850)
    Thetase500 = ThetaSe(T500, Td500, 500)
    Icondi = Thetase500 - Thetase850
    return round(Icondi, 3)


def DCI_index(data_press, data_temp, data_dewp):
    """
    计算深对流指数
    :param data_press: 微波辐射计高度数据集
    :param data_temp: 微波辐射计温度数据集
    :param data_dewp: 微波辐射计相对湿度数据集
    """
    T850 = get_temp(850, data_press, data_temp)
    Td850 = get_dewp(850, data_press, data_dewp)
    LI = LI_index(data_press, data_temp, data_dewp)
    DCI = T850 + Td850 - LI
    return round(DCI, 3)
