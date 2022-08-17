# -*- coding:utf-8 -*-
__author__ = 'wanglei'

import json
import math
from itertools import groupby
import os
import numpy as np
import pandas as pd
from sympy import *
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline

microwave_height = (0, 10, 25, 50, 75, 100, 130, 160, 190, 220, 250, 280, 310, 340, 370, 400, 430, 460, 490, 520,
                    560, 600, 640, 680, 720, 760, 800, 840, 880, 920, 960, 1000, 1040, 1080, 1120, 1160, 1200, 1260,
                    1320, 1380, 1440, 1500, 1560, 1620, 1680, 1740, 1800, 1890, 1980, 2170, 2260, 2350, 2430, 2500, 2600,
                    2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3650, 3800, 3950, 4100, 4250, 4400, 4550, 4600,
                    4800, 5000, 5200, 5400, 5600, 5800, 6000, 6300, 6600, 6900, 7200, 7500, 7800, 8100, 8400, 8700, 9000,
                    9300, 9600, 9800, 10000)


def inter_un(xi, yi, xnew):
    """
    插值，范围外插,(xi必须严格单调递增)
    """
    f = InterpolatedUnivariateSpline(xi, yi, k=1)
    return [round(i, 3) for i in list(f(xnew))]


def get_y(data_hei, data, x):
    """
    给定高度，根据插值函数f求对应高度的数据
    :param data_hei: 高度数据（递增）
    :param x: 已知高度值
    :return: 已知高度值对应的数据（温度，湿度等）
    """
    f = interp1d(data_hei, data, kind='cubic')
    y = f(x)
    return round(float(y), 3)


def get_height(press):
    """
    给定压强计算高度
    :param press: 压强
    :return: 高度
    """
    height = 44331 * (1 - (press / 1013.25) ** 0.1903)
    return round(height, 3)


def get_hei(press):
    """
    给定压强计算高度
    :param press: 压强
    :return: 高度
    """
    return round(44331 * (1 - (press / 1013.25) ** 0.1903), 3)


def get_temp(press, data_hei, data_temp):
    """
    给定压强，反算出对应温度(℃）
    :param press: 压强
    :param data_hei: 高度数据集
    :param data_temp: 温度数据集
    :return: 温度
    """
    x = get_height(press)
    temp = get_y(data_hei, data_temp, x)
    return round(temp, 3)


def get_rh(press, data_hei, data_rh):
    """
    给定压强，反算出对应相对湿度(%)
    :param press: 压强
    :param data_hei: 高度数据集
    :param data_rh: 相对湿度数据集
    :return: 相对湿度
    """
    x = get_height(press)
    rh = get_y(data_hei, data_rh, x)
    return round(rh, 3)


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


def get_dewp(data_temp, data_rh):
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


def T(t):
    """
    将摄氏温度转换为开氏温度
    """
    return t + 273.15


def r_s(P, t):
    """
    获得饱和混合比(g/g)
    :param P:气压
    :param t:气温(℃)
    :return:饱和混合比
    """
    # 622 * e_s(t) / P    单位为(g/kg)
    return 0.622 * get_es(t) / P


def get_r(e, p):
    """
    获得水汽混合比 （g/g）        （天气分析 P9）
    :param e:水汽压
    :param p:气压
    """
    return 0.622 * e / (p - e)


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


def get_Lv(t):
    """
    汽化潜热
    :param t: 温度
    :return: 汽化潜热
    """
    return 2.5 * 10 ** 6 - 2323 * t  # 大气物理学P129-130


def get_L(t):
    """
    给定温度得到潜热L（J/kg)
    :param tempc:温度（℃）
    """
    L = (2500 - 2.4 * t) * 1000
    return L


def P_d(P, t):
    """
    获得系统干空气分压
    :param P: 气压
    :param t: 气温
    """
    return P - get_es(t)


def getrh(t, dewp, pres):
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

    q = get_r((E * getrh(t, Td, pres)) / 100, pres)  # 比湿
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


K = 0.286  # python常量：约定俗成，不可更改，全部是大写字母
C_pd = 1004


def e_s(t):
    if t > 0:
        return 6.11 * 10 ** ((7.5 * t) / (237.3 + t))
    else:
        return 6.11 * 10 ** ((9.5 * t) / (265.5 + t))


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


def lcl(t, dewp):
    """
    获得抬升凝结高度
    :param tewp:地面温度
    :param dewp:地面露点
    """
    z = 123 * (t - dewp)
    return z


def thi(temp, dewp):
    """
    给定温度和露点获得人体舒适度指数
    :param temp: 温度（℃）
    :param dewp: 露点（℃）
    """
    THI = temp - 1.0799 * exp(0.03755 * temp) * (1 - exp(0.0801 * (dewp - 14)))
    return THI


def get_S(data_height, data_temp, data_rh):
    """
    计算S指数，用于表明雷暴发生潜力指数，常用于4月份到9月份雷暴发生可能性预测；
    :param data_height: 高度数据(列表)
    :param data_temp: 温度数据(列表)
    :param data_rh: 湿度数据(列表)
    """
    T500 = get_temp(500, data_height, data_temp)
    T700 = get_temp(700, data_height, data_temp)
    T850 = get_temp(850, data_height, data_temp)
    f700 = get_rh(850, data_height, data_rh)
    Td700 = get_dewp(T700, f700)
    TT = get_TT(data_height, data_temp, data_rh)
    VT = T850 - T500
    if VT >= 25:
        K = 0
    elif VT > 22 and VT < 25:
        K = 2
    elif VT <= 22:
        K = 6
    S = TT - (T700 - Td700) - K
    return round(S, 3)


def get_TT(data_height, data_temp, data_rh):
    """
    计算全总指数(TT)
    :param data_height: 高度数据(列表)
    :param data_temp: 温度数据(列表)
    :param data_rh: 湿度数据(列表)
    """
    T500 = get_temp(500, data_height, data_temp)
    T850 = get_temp(850, data_height, data_temp)
    f850 = get_rh(850, data_height, data_rh)
    Td850 = get_dewp(T850, f850)
    TT = T850 + Td850 - 2 * T500
    return round(TT, 3)


def get_TQ(data_height, data_temp, data_rh):
    """
    计算TQ指数，用于评估底层潜在对流情况，
    :param data_height:高度数据(列表)
    :param data_temp:温度数据(列表)
    :param data_rh:湿度数据(列表)
    """
    T700 = get_temp(700, data_height, data_temp)
    T850 = get_temp(850, data_height, data_temp)
    f850 = get_rh(850, data_height, data_rh)
    Td850 = get_dewp(T850, f850)
    TQ = (T850 + Td850) - 1.7 * T700
    return round(TQ, 3)


def get_CT(data_height, data_temp, data_rh):
    """
    计算CT交叉总指数，用于重要天气过程指数
    :param data_height: 高度数据集
    :param data_temp: 温度数据集
    :param data_rh: 湿度数据集
    """
    T500 = get_temp(500, data_height, data_temp)
    T850 = get_temp(850, data_height, data_temp)
    f850 = get_rh(850, data_height, data_rh)
    Td850 = get_dewp(T850, f850)
    CT = Td850 - T500
    return round(CT, 3)


def get_JI(data_height, data_temp, data_rh):
    """
    计算杰弗森指数（JI）
    :param data_height:高度数据
    :param data_temp:温度数据
    :param data_rh:湿度数据
    """
    T500 = get_temp(500, data_height, data_temp)
    T700 = get_temp(700, data_height, data_temp)
    f700 = get_rh(700, data_height, data_rh)
    Td700 = get_dewp(T700, f700)
    T850 = get_temp(850, data_height, data_temp)
    f850 = get_rh(850, data_height, data_rh)
    Td850 = get_dewp(T850, f850)
    thetaw850 = get_theta(get_tw(T850, Td850, 850), 850)
    JI = 1.6 * thetaw850 - T500 - 0.5 * (T700 - Td700) - 8
    return round(JI, 3)


def get_VT(data_height, data_temp):
    """
    计算VT垂直总指数，通过温度和湿度来预测重要天气过程
    :param data_height: 高度数据集
    :param data_temp: 温度数据集
    """
    T500 = get_temp(500, data_height, data_temp)
    T850 = get_temp(850, data_height, data_temp)
    VT = T850 - T500
    return round(VT, 3)


def get_k(data_hei, data_temp, data_rh):
    """
    给定高度数据集和温度数据集以及相对湿度数据集，计算K指数
    :param T850: 800hPa压强下温度        mb = mbar毫巴 = hPa百帕
    :param Td850:800hPa压强下露点温度
    :param T700:700hPa压强下温度
    :param Td700:700hPa压强下露点温度
    :param T500:500hPa压强下温度
    :return:K指数
    """
    T500 = get_temp(500, data_hei, data_temp)
    T700 = get_temp(700, data_hei, data_temp)
    f700 = get_rh(700, data_hei, data_rh)
    Td700 = get_dewp(T700, f700)
    T850 = get_temp(850, data_hei, data_temp)
    f850 = get_rh(850, data_hei, data_rh)  # 850压强下的相对湿度
    Td850 = get_dewp(T850, f850)  # 850压强下的露点温度
    K = T850 - T500 + Td850 - T700 + Td700
    return round(K, 3)


def get_si(data_hei, data_temp, data_rh):
    """
    给定高度数据集和温度数据集以及相对湿度数据集，求出沙瓦特指数
    :param data_hei: 高度数据集
    :param data_temp: 温度数据集
     :param data_rh: 相对湿度数据集
    :return: 沙瓦特指数
    """
    T850 = get_temp(850, data_hei, data_temp)
    f850 = get_rh(850, data_hei, data_rh)
    Td850 = get_dewp(T850, f850)
    tc = get_tc(T850, Td850)  # 抬升凝结温度
    pc = get_pc(T850, 850, tc)  # 抬升凝结压强
    T500 = get_temp(500, data_hei, data_temp)
    T_se_start = T_se(pc, tc)  # 抬升凝结高度处对应的初始假相当位温
    while pc > 500:
        pc = pc - 1
        ret_T_1 = T_1(pc, tc, T_se_start)
        tc = ret_T_1
        if pc <= 500:
            return round(T500 - tc, 3)


def get_A(data_hei, data_temp, data_rh):
    """
     给定高度数据集和温度数据集以及相对湿度数据集，计算A指数
    :param data_hei: 高度数据集
    :param data_temp: 温度数据集
    :param data_rh: 相对湿度数据集
    """
    T500 = get_temp(500, data_hei, data_temp)
    f500 = get_rh(500, data_hei, data_rh)
    Td500 = get_dewp(T500, f500)
    T850 = get_temp(850, data_hei, data_temp)
    f850 = get_rh(850, data_hei, data_rh)
    Td850 = get_dewp(T850, f850)
    T700 = get_temp(700, data_hei, data_temp)
    f700 = get_rh(700, data_hei, data_rh)
    Td700 = get_dewp(T700, f700)
    A = (T850 - T500) - ((T850 - Td850) + (T700 - Td700) + (T500 - Td500))
    return round(A, 3)


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


def insert_moist_start_data(press, moist_adiabat_P, height, t):
    """将湿绝热线的起始压强和温度添加到压强列表和温度列表中以保证数据都是相同压强下的数据"""
    if press[0] > moist_adiabat_P[0]:  # 将湿绝热线的起始压强和温度添加到压强列表中
        for i in range(len(height)):
            if press[i] < moist_adiabat_P[0]:
                press.insert(i, moist_adiabat_P[0])
                t.insert(i, None)
                break
            else:
                pass
    elif press[0] < moist_adiabat_P[0]:
        press.insert(0, moist_adiabat_P[0])
        t.insert(0, None)
    else:
        pass


def insert_moist_adiabat_data(moist_adiabat_P, moist_adiabat_T, press):
    """得到插值后的湿绝热线温度值，与压强一一对应"""
    xi = sorted(moist_adiabat_P)
    yi = sorted(moist_adiabat_T)
    moist_adiabat_data = inter_un(xi, yi, press)
    return moist_adiabat_data


def get_LfcAndElAndCape(press, temp, mois):
    """
    计算温度和湿绝热线交点的对应高度LFC和EL的压强
    """
    try:
        index = temp.index(None) + 1
    except:
        index = 0
    f1 = interp1d(press[index:], temp[index:], "cubic")
    f2 = interp1d(press[index:], mois[index:], "cubic")
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


def get_LI(data_hei, data_temp, data_rh):
    """
    给定高度数据集和温度数据集以及相对湿度数据集，求出抬升指数
    :param data_hei: 高度数据集
    :param data_temp: 温度数据集
    :param data_rh: 相对湿度数据集
    2019.7.6(此处的假相当位温仍就是未订正前的)
    """
    P900 = get_press(900)  # 得到900m高度的气压
    T = get_temp(P900, data_hei, data_temp)
    f = get_rh(P900, data_hei, data_rh)
    Td = get_dewp(T, f)
    tc = get_tc(T, Td)  # 抬升凝结温度
    pc = get_pc(T, P900, tc)  # 抬升凝结压强
    T500 = get_temp(500, data_hei, data_temp)
    T_se_start = T_se(pc, tc)
    while pc > 500:
        pc = pc - 1
        ret_T_1 = T_1(pc, tc, T_se_start)
        tc = ret_T_1
        if pc <= 500:
            return round(T500 - tc, 3)


def get_TI(data_hei, data_temp, data_rh):
    """
    给定高度数据集和温度数据集以及相对湿度数据集，计算TI指数(Thompson指数)
    :param data_hei:高度数据集
    :param data_temp:温度数据集
    :param data_rh:相对湿度数据集
    """
    KI = get_k(data_hei, data_temp, data_rh)
    LI = get_LI(data_hei, data_temp, data_rh)
    TI = KI - LI
    return round(TI, 3)


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


def DelThetaE(data_hei, data_temp, data_rh):
    # def DelThetaE(t1, dewp1, p1, data_hei, data_temp, data_rh):
    """
    根据温度、露点、气压，计算相当位温差, 地面相当位温和中层相当位温的差值, 中层为300hPa高度层
    :param t1: 地面温度(°C)
    :param dewp1: 地面露点(°C)
    :param p1: 地面气压(mb)   1013.25
    :param t2: 300hPa高度层的温度(°C)
    :param dewp2: 300hPa高度层的露点(°C)
    :param p2: 300hPa

    2019.7.12 将相当位温修改为假相当位温，数据貌似合理。

    """
    t0 = data_temp[0]
    f0 = data_rh[0]
    td0 = get_dewp(t0, f0)
    # thetae0 = Thetae(t0, td0, 1013.25)
    thetae0 = ThetaSe(t0, td0, 1013.25)
    t300 = get_temp(300, data_hei, data_temp)
    f300 = get_rh(300, data_hei, data_rh)
    td300 = get_dewp(t300, f300)
    # thetae1 = Thetae(t300, td300, 300)
    thetae1 = ThetaSe(t300, td300, 300)
    delthetae = abs(thetae0 - thetae1)
    return round(delthetae, 3)


def get_GOES(data_hei, data_temp, data_rh):
    """
    给定高度数据集和温度数据集以及相对湿度数据集，计算GOES混合微下击暴流指数
    :param data_hei: 高度数据集
    :param data_temp: 温度数据集
    :param data_rh: 相对湿度数据集
    """
    T670 = get_temp(670, data_hei, data_temp)
    f670 = get_rh(670, data_hei, data_rh)
    Td670 = get_dewp(T670, f670)
    T850 = get_temp(850, data_hei, data_temp)
    f850 = get_rh(850, data_hei, data_rh)
    Td850 = get_dewp(T850, f850)
    h850 = get_hei(850)  # 1457.585
    h670 = get_hei(670)  # 3355.742      h670- h850=1898.157
    h = (h670 - h850) / 1000  # 高度差，将单位m换算成km
    G = (T850 - T670) / h
    HMI = G + (T850 - Td850) - (T670 - Td670)
    return round(HMI, 3)


def get_DCI(data_hei, data_temp, data_rh):
    """
    计算深对流指数
    :param data_hei: 高度数据集
    :param data_temp: 温度数据集
    :param data_rh: 相对湿度数据集
    """
    T850 = get_temp(850, data_hei, data_temp)
    f850 = get_rh(850, data_hei, data_rh)
    Td850 = get_dewp(T850, f850)
    LI = get_LI(data_hei, data_temp, data_rh)
    DCI = T850 + Td850 - LI
    return round(DCI, 3)


def get_MVV(cape):
    """
    计算最大垂直速度
    :param cape:cape正值和
    """
    MVV = sqrt(2 * abs(cape))
    return round(MVV, 3)


def get_Icondi(data_height, data_temp, data_rh):
    """
    计算条件性稳定度指数      计算公式中是假相当位温的标志，C#程序中计算的是相当位温，此处先用相当位温计算的

    2019.7.6（王磊） 使用订正后的假相当位温计算

    :param data_height:微波辐射计高度数据集
    :param data_temp: 微波辐射计温度数据集
    :param data_rh: 微波辐射计相对湿度数据集
    """
    T850 = get_temp(850, data_height, data_temp)
    f850 = get_rh(850, data_height, data_rh)
    Td850 = get_dewp(T850, f850)
    T500 = get_temp(500, data_height, data_temp)
    f500 = get_rh(500, data_height, data_rh)
    Td500 = get_dewp(T500, f500)
    Thetase850 = ThetaSe(T850, Td850, 850)
    Thetase500 = ThetaSe(T500, Td500, 500)
    Icondi = Thetase500 - Thetase850
    return round(Icondi, 3)


def get_KOI(data_height, data_temp, data_rh):
    """
    计算KO指数，用于评估雷暴发生指数       # 根据C#程序翻译的
    :param data_height: 微波辐射计高度数据集
    :param data_temp: 微波辐射计温度数据集
    :param data_rh: 微波辐射计相对湿度数据集
    """
    T1000 = get_temp(1000, data_height, data_temp)
    f1000 = get_rh(1000, data_height, data_rh)
    Td1000 = get_dewp(T1000, f1000)
    T850 = get_temp(850, data_height, data_temp)
    f850 = get_rh(850, data_height, data_rh)
    Td850 = get_dewp(T850, f850)
    T700 = get_temp(700, data_height, data_temp)
    f700 = get_rh(700, data_height, data_rh)
    Td700 = get_dewp(T700, f700)
    T500 = get_temp(500, data_height, data_temp)
    f500 = get_rh(500, data_height, data_rh)
    Td500 = get_dewp(T500, f500)
    thetae1000 = ThetaSe(T1000, Td1000, 1000)
    thetae850 = ThetaSe(T850, Td850, 850)
    thetae700 = ThetaSe(T700, Td700, 700)
    thetae500 = ThetaSe(T500, Td500, 500)
    KOI = 0.5 * (thetae500 + thetae700) - 0.5 * (thetae850 + thetae1000)
    return round(KOI, 3)


def get_MDPI(data_height, data_temp, data_rh):
    """
    计算微下击暴流潜势日指数
    :param data_height: 微波辐射计高度数据集
    :param data_temp: 微波辐射计温度数据集
    :param data_rh: 微波辐射计相对湿度数据集
    """
    data_pres = get_press(data_height)  # 根据微波辐射计的高度列表，求出压强列表
    thetae1 = []  # 存放近地层150mb内的相当位温列表
    thetae2 = []  # 存放650至500mb中的相当位温列表
    for each in data_pres:
        if each >= 850:  # 压强大于850，属于近地层150mb内，计算最大相当位温
            T = get_temp(each, data_height, data_temp)
            rh = get_rh(each, data_height, data_rh)
            Td = get_dewp(T, rh)
            thetae = ThetaSe(T, Td, each)
            thetae1.append(thetae)
        elif each >= 500 and each <= 650:
            T = get_temp(each, data_height, data_temp)
            rh = get_rh(each, data_height, data_rh)
            Td = get_dewp(T, rh)
            thetae = ThetaSe(T, Td, each)
            thetae2.append(thetae)
    maxthetae = max(thetae1)
    minthetae = min(thetae2)
    MDPI = (maxthetae - minthetae) / 20
    return round(MDPI, 3)


def get_WMSI(data_height, data_temp, data_rh, cape):
    """
    计算湿下击暴流严重性指数
    :param data_height: 微波辐射计高度数据集
    :param data_temp: 微波辐射计温度数据集
    :param data_rh: 微波辐射计相对湿度数据集
    :param cape: 对流有效位能
    """
    delthetae = DelThetaE(data_height, data_temp, data_rh)
    WMSI = (abs(cape) * delthetae) / 1000
    return round(WMSI, 3)


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


# =============================================================================================
# =============================================================================================
# =============================================================================================


def get_df(df_lv2, qc):
    if qc:
        global microwave_height
        microwave_height = (
            0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250,
            3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000, 6250, 6500, 6750, 7000, 7250, 7500, 7750,
            8000, 8250, 8500, 8750, 9000, 9250, 9500, 9750, 10000)
        col = list(df_lv2.columns)
        col = col[:11]
        new_col = col + [str(i) for i in microwave_height]
        df_lv2.columns = new_col
        # 0.000km不能写死，防止文件表头变了，可以重新设置表头,选取0km开始的列
        df_lv2.loc[df_lv2["10"] == 11, "0":] = df_lv2.loc[df_lv2["10"] == 11, "0":] - 273.15
        df_lv2 = df_lv2.round(2)


    lv2_11 = df_lv2.loc[df_lv2["10"] == 11].values.tolist()
    lv2_12 = df_lv2.loc[df_lv2["10"] == 12].values.tolist()
    lv2_13 = df_lv2.loc[df_lv2["10"] == 13].values.tolist()
    return lv2_11, lv2_12, lv2_13


def handle_df(df_11, df_12, df_13, col):
    df_11.drop(columns=col, inplace=True)
    df_12.drop(columns=col, inplace=True)
    df_13.drop(columns=col, inplace=True)


def create_json_file(result, filepath):
    with open(filepath, 'w') as f:
        json.dump(result, f, ensure_ascii=False)


def S_Index(lv2_11, lv2_13):
    """S指数"""
    height = list(microwave_height)
    result = {
        "unit": None,
        "label": "S指数，表明雷暴发生潜力指数,根据温度列表和相对湿度列表数据计算TT、T850（850压强下温度）、T500（500压强下温度）、T700（700压强下温度）、"
                 "TD700（700压强下露点温度）等，计算公式为S=TT-(T700-TD700)-K，其中，当VT≥25时K=0，当22<VT<25,时K=2，当VT≤22,K=6, VT= T_850-T_500。"
    }
    list_data = []
    for i in range(len(lv2_11)):
        data = get_S(height, lv2_11[i][11:], lv2_13[i][11:])
        meta = {
            "x": lv2_11[i][1],
            "y": data
        }
        list_data.append(meta)
    result["data"] = list_data
    return result


def TT_Index(lv2_11, lv2_13):
    """全总指数"""
    height = list(microwave_height)
    result = {
        "unit": "℃",
        "label": "全总指数，根据温度列表和相对湿度列表数据计算T850(850压强下温度),Td850（850压强下露点温度），T500（500压强下温度），"
                 "计算公式为TT=T850+Td850-2*T500，常用于重要天气过程的指示，TT越大越容易发生强对流天气。"
    }
    list_data = []
    for i in range(len(lv2_11)):
        data = get_TT(height, lv2_11[i][11:], lv2_13[i][11:])
        meta = {
            "x": lv2_11[i][1],
            "y": data
        }
        list_data.append(meta)
    result["data"] = list_data
    return result


def TQ_Index(lv2_11, lv2_13):
    """TQ指数"""
    height = list(microwave_height)
    result = {
        "unit": "℃",
        "label": "TQ指数，用于评估底层潜在对流情况，计算公式为TQ=(T850+Td850)-1.7*T700。T850（850压强下温度），Td850（850压强下露点），T700（700压强下温度）"
    }
    list_data = []
    for i in range(len(lv2_11)):
        data = get_TQ(height, lv2_11[i][11:], lv2_13[i][11:])
        meta = {
            "x": lv2_11[i][1],
            "y": data
        }
        list_data.append(meta)
    result["data"] = list_data
    return result


def CT_index(lv2_11, lv2_13):
    """交叉总指数"""
    height = list(microwave_height)
    result = {
        "unit": "℃",
        "label": "交叉总指数Cross Totals(CT)Index，常用于重要天气过程指数，计算公式为CT=Td850-T500, Td850(850压强下露点温度)，T500（500压强下温度）"
    }
    list_data = []
    for i in range(len(lv2_11)):
        data = get_CT(height, lv2_11[i][11:], lv2_13[i][11:])
        meta = {
            "x": lv2_11[i][1],
            "y": data
        }
        list_data.append(meta)
    result["data"] = list_data
    return result


def get_index(index_function, lv2_11=None, lv2_12=None, lv2_13=None, unit=None, label=None):
    height = list(microwave_height)
    result = {
        "unit": unit,
        "label": label
    }
    list_data = []
    if lv2_11 is not None and lv2_13 is not None and lv2_12 is None:
        for i in range(len(lv2_11)):
            data = index_function(height, lv2_11[i][11:], lv2_13[i][11:])
            meta = {
                "x": lv2_11[i][1],
                "y": data
            }
            list_data.append(meta)
        result["data"] = list_data
        return result
    elif lv2_11 is not None and lv2_12 is None and lv2_13 is None:
        for i in range(len(lv2_11)):
            data_vt = get_VT(height, lv2_11[i][11:])
            meta = {
                "x": lv2_11[i][1],
                "y": data_vt
            }
            list_data.append(meta)
        result["data"] = list_data
        return result


def JI_index(lv2_11, lv2_13):
    """杰弗森指数"""
    """无量纲"""
    label = "杰弗森指数, 计算公式为JI=1.6*WBPT850-T500-0.5*(T700-Td700)-8, WBPT为湿球潜温"
    JI = get_index(get_JI, lv2_11=lv2_11, lv2_13=lv2_13, label=label)
    return JI


def VT_index(lv2_11):
    """垂直总指数"""
    unit = "℃"
    label = "垂直总指数Vertical Totals(VT) Index，通过温度和湿度来预测重要天气过程，计算公式为VT=T850-T500。"
    VT = get_index(get_VT, lv2_11=lv2_11, label=label, unit=unit)
    return VT


def K_index(lv2_11, lv2_13):
    """K指数"""
    unit = "℃"
    label = "K Index（K指数），用于雷暴潜在发生指数，反映大气的层结稳定度情况，K指数越大，层结越不稳定，K=(T850-T500)+Td850-(T700-Td700)。"
    K_index = get_index(get_k, lv2_11=lv2_11, lv2_13=lv2_13, label=label, unit=unit)
    return K_index


def SI_index(lv2_11, lv2_13):
    """沙瓦特指数"""
    unit = "℃"
    label = "肖瓦特指数, 定义为850hPa等压面上的湿空气块沿干绝热线抬升，到达抬升凝结高度后再沿湿绝热线上升至500bPa时具有的气块温度与500bPa等压面上的环境温度的差值"
    SI = get_index(get_si, lv2_11=lv2_11, lv2_13=lv2_13, label=label, unit=unit)
    return SI


def A_index(lv2_11, lv2_13):
    """A指数"""
    unit = "℃"
    label = "A指数，计算公式为A=(T850-T500 )-[(T850-Td850 )+(T700-Td700 )+(T500-Td500)]，A 值越大, 表明大气越不稳定或对流层中下层饱和程度越高, 越有利于产生降水。"
    A = get_index(get_A, lv2_11=lv2_11, lv2_13=lv2_13, label=label, unit=unit)
    return A


def cape(lv2_11, lv2_13):
    """对流有效位能"""
    height = list(microwave_height)
    result = {
        "unit": 'J/kg',
        "label": "对流有效位能,根据温度列表和相对湿度列表数据，首先计算出湿绝热线，再根据其找出自由对流高度和平衡高度，取两个高度之间的积分和即cape值"
    }
    data = []
    for i in range(len(lv2_11)):
        press = get_press(height)
        t = lv2_11[i][11:]
        r = lv2_13[i][11:]
        time = lv2_11[i][1]
        dewp = get_dewp(t[0], r[0])
        tc = get_tc(t[0], dewp)
        pc = get_pc(t[0], press[0], tc)
        moist_adiabat_P, moist_adiabat_T = moist_adiabat(pc, tc)
        moist_adiabat_data = insert_moist_adiabat_data(moist_adiabat_P, moist_adiabat_T, press)  # 插值后的湿绝热线温度值
        hei_important_data = get_LfcAndElAndCape(press, t, moist_adiabat_data)
        meta = {
            "x": time,
            "y": hei_important_data.get('cape', 0)
        }
        data.append(meta)
    result["data"] = data
    return result


def LI_index(lv2_11, lv2_13):
    """抬升指数"""
    unit = "℃"
    label = "抬升指数，表明条件性稳定度指数，定义为平均气块从修正的低层900m高度沿干绝热线上升，到达凝结高度后再沿湿绝热线上升至500hPa时所具有的温度与500hPa等压面上的环境温度的差值"
    LI = get_index(get_LI, lv2_11=lv2_11, lv2_13=lv2_13, label=label, unit=unit)
    return LI


def TI_index(lv2_11, lv2_13):
    """Thompson指数"""
    label = "TI指数，常用于指示雷暴发生的可能性,TI=K-LI，K（K指数），LI（抬升指数）"
    TI = get_index(get_TI, lv2_11=lv2_11, lv2_13=lv2_13, label=label, unit="℃")
    return TI


def delThetae(lv2_11, lv2_13):
    """相当位温差"""
    unit = "K"
    label = "相当位温差，用于评估“湿下击暴流”的潜在发生指数，是地面相当位温和中层相当位温的差值，该方法定义，中层为300hPa高度层。"
    delThetaE = get_index(DelThetaE, lv2_11=lv2_11, lv2_13=lv2_13, unit=unit, label=label)
    return delThetaE


def GOES_index(lv2_11, lv2_13):
    """GOES混合微下击暴流指数"""
    label = "GOES混合微下击暴流指数, HMI=G+(T850-Td850)-(T670-Td670), 其中G表示850mb到670mb每千米温度递减率"
    GOES = get_index(get_GOES, lv2_11=lv2_11, lv2_13=lv2_13, label=label)
    return GOES


def DCI_index(lv2_11, lv2_13):
    """深对流指数"""
    label = "深对流指数（DCI）,DCI=T850+Td850-LI, T850：850压强下的温度；Td850：850压强下的露点温度；LI：抬升指数"
    DCI = get_index(get_DCI, lv2_11=lv2_11, lv2_13=lv2_13, label=label, unit="℃")
    return DCI


def MVV_index(lv2_11, lv2_13):
    """最大垂直速度"""
    result = {
        "unit": "m/s",
        "label": "最大垂直速度(MVV), 潜在对流上升气流的最大垂直速度, 计算公式math.sqrt(2*cape)"
    }
    height = list(microwave_height)
    data = []
    for i in range(len(lv2_11)):
        press = get_press(height)
        t = lv2_11[i][11:]
        r = lv2_13[i][11:]
        time = lv2_11[i][1]
        dewp = get_dewp(t[0], r[0])
        tc = get_tc(t[0], dewp)
        pc = get_pc(t[0], press[0], tc)
        moist_adiabat_P, moist_adiabat_T = moist_adiabat(pc, tc)
        moist_adiabat_data = insert_moist_adiabat_data(moist_adiabat_P, moist_adiabat_T, press)  # 插值后的湿绝热线温度值
        hei_important_data = get_LfcAndElAndCape(press, t, moist_adiabat_data)
        cape = hei_important_data.get('cape', 0)
        MVV = get_MVV(cape)
        meta = {
            "x": time,
            "y": MVV
        }
        data.append(meta)
    result['data'] = data
    return result


def Icondi_index(lv2_11, lv2_13):
    """条件性稳定度指数"""
    label = "条件性稳定度指数，I_Condi=θ_se500-θ_se850， 850hPa代表起始高度，500hPa代表上层"
    Icondi = get_index(get_Icondi, lv2_11=lv2_11, lv2_13=lv2_13, label=label, unit="℃")
    return Icondi


def KOI_index(lv2_11, lv2_13):
    """KO指数"""
    label = "KO指数, KO=(θe_500+θe_700)/2-(θe_850+θe_1000)/2, θe为等潜温"
    KOI = get_index(get_KOI, lv2_11=lv2_11, lv2_13=lv2_13, label=label)
    return KOI


def MDPI_index(lv2_11, lv2_13):
    """微下击暴流潜势日指数"""
    label = "微下击暴流潜势日指数,是近地层150mb内（850压强内）的最大相当位温和650至500mb中最小的相当位温的函数，若MDPI≥1则认为有微下击暴流。"
    MDPI = get_index(get_MDPI, lv2_11=lv2_11, lv2_13=lv2_13, label=label)
    return MDPI


def WMSI_index(lv2_11, lv2_13):
    """湿下击暴流严重性指数"""
    height = list(microwave_height)
    result = {
        "unit": None,
        "label": "湿下击暴流严重性指数,WMSI=(CAPE*DeltaePT)/1000"
    }
    data = []
    for i in range(len(lv2_11)):
        press = get_press(height)
        t = lv2_11[i][11:]
        r = lv2_13[i][11:]
        time = lv2_11[i][1]
        dewp = get_dewp(t[0], r[0])
        tc = get_tc(t[0], dewp)
        pc = get_pc(t[0], press[0], tc)
        moist_adiabat_P, moist_adiabat_T = moist_adiabat(pc, tc)
        moist_adiabat_data = insert_moist_adiabat_data(moist_adiabat_P, moist_adiabat_T, press)  # 插值后的湿绝热线温度值
        hei_important_data = get_LfcAndElAndCape(press, t, moist_adiabat_data)
        cape = hei_important_data.get("cape", 0)
        WMSI = get_WMSI(height, t, r, cape)
        meta = {
            'x': time,
            "y": WMSI
        }
        data.append(meta)
    result['data'] = data
    return result


def lcl_index(lv2_11, lv2_13):
    result = {"unit": "m",
              "label": "抬升凝结高度,根据LV2数据微文件中的温度和相对湿度数据计算出地面的露点和温度，再求出lcl"}
    data = []
    for i in range(len(lv2_11)):
        dewp = get_dewp(lv2_11[i][11:], lv2_13[i][11:])
        temp_0 = lv2_11[i][11]
        dewp_0 = dewp[0]
        lcl = 123 * (temp_0 - dewp_0)
        meta = {
            "x": lv2_11[i][1],
            "y": round(lcl, 3)
        }
        data.append(meta)
    result["data"] = data
    return result


def thi_index(lv2_11, lv2_13):
    result = {"unit": None,
              "label": "人体舒适度指数,分别取地面的温度和露点温度计算"}
    list_data = []
    for i in range(len(lv2_11)):
        dewp = get_dewp(lv2_11[i][11:], lv2_13[i][11:])
        temp_0 = lv2_11[i][11]
        dewp_0 = dewp[0]
        data = thi(temp_0, dewp_0)
        meta = {
            "x": lv2_11[i][1],
            "y": round(data, 3)
        }
        list_data.append(meta)
    result["data"] = list_data
    return result


def pblh(LV1_path):
    """边界层高度"""
    if LV1_path is None:
        ex = Exception('lack of LV1 file path')
        raise ex

    col = [1, 8, 11, 17, 19, 22, 27, 28, 30, 33, 35, 38, 39, 41]
    coe = (-24.69, 49.80, -123.93, 40.14, 22.45, 95.61, -45.3, -29.95, 96.63, -252.47, -41.84, 1.54, 215.59)
    b = 5497.46
    try:
        df_lv1 = pd.read_csv(LV1_path, encoding="GBK", usecols=col)
    except:
        df_lv1 = pd.read_csv(LV1_path, encoding="utf-8", usecols=col)
    df_lv1["pblh"] = np.sum(df_lv1.iloc[:, 1:] * np.array(coe), axis=1) + b
    result = {"unit": "m",
              "label": "边界层高度，利用LV1数据文件14个通道亮温值计算"}
    data = []
    for index in df_lv1.index:
        time = df_lv1.loc[index][0]
        pblh_val = df_lv1.loc[index][-1]
        meta = {
            'x': time,
            'y': round(pblh_val, 3)
        }
        data.append(meta)
    result["data"] = data
    return result


def inver_height(LV2_11):
    """逆温层高度"""
    result = []
    height = list(microwave_height)
    for i in LV2_11:
        time = i[1]
        data = i[11:]
        meta = {
            'x': time,
            'y': None
        }
        for index, val in enumerate(data):
            if index > 0 and index < (len(LV2_11) - 1):
                if val > data[index - 1] and val >= data[index + 1]:
                    meta['y'] = height[index]
                    break
        result.append(meta)
    return result


def inver_height_new(LV2_11):
    """逆温层(新)"""
    result = {"unit": "m",
              "label": "逆温层，利用微波辐射计LV2温度数据进行计算，分别取温度随高度增加而递增的高度起始点和终止点"}
    data = []
    height = list(microwave_height)
    for i in LV2_11:
        time = i[1]
        data_temp = i[11:]
        print(data_temp)
        meta = {
            'x': time,
            'y': []
        }
        # for index, val in enumerate(data):
        #     location = 0
        #     if index > 0 and index < (len(data) - 1):
        # while location==0:
        #     if val > data[index - 1] and val >= data[index + 1]:
        #         # location = index
        #         for index, val in enumerate(data[location:]):
        #             if val < data[index - 1] and val <= data[index + 1]:
        #                 if abs(height[location]-val)>1:
        #                     hei_list = [height[location], height[index]]
        #                     hei_all.append(hei_list)
        #                     location = index

        # hei_list = []
        # count = 0
        # max_sign=0
        # min_sign=0
        # max_hei = []
        # min_hei = []
        # for index, val in enumerate(data):
        #     if index > 0 and index < (len(data) - 1):
        #         if val > data[index - 1] and val >= data[index + 1]:
        #             max_hei.append(index)
        #             count+=1
        #             max_sign=count
        #             print("max_sin: ", max_sign)
        #             print("da", index)
        #         if val < data[index - 1] and val <= data[index + 1]:
        #             min_hei.append(index)
        #             count += 1
        #             min_sign=count
        #             print("min_sin: ", min_sign)
        #             print("xiao", index)
        # # print(hei_list)
        # print(max_hei)
        # print(min_hei)
        # return hei_list
        all_point = []  # 出现逆温层的索引列表
        for index, val in enumerate(data_temp):
            if index == len(data_temp) - 1:
                break
            if val <= data_temp[index + 1]:
                all_point.append(index)
                # if val==data[index+1]:
                #     all_point.append(index+1)

        con_digits = continuous_digits(all_point)
        y_list = []
        delta_h, delta_t, h0 = 0, 0, None
        for i in con_digits:
            diff_val = round(data_temp[max(i)] - data_temp[min(i)], 2)
            if height[min(i)] < 1000:
                delta_t += diff_val
                delta_h += (height[max(i)] - height[min(i)])
                if h0 is None:
                    h0 = height[min(i)]

                y_list_meta = [height[min(i)], height[max(i)]]
                y_list.append(y_list_meta)

        # meta['y'] = y_list
        meta['y'] = [delta_h, delta_t, h0]
        data.append(meta)
    result['data'] = data
    return result


def continuous_digits(lst):
    list_all = []
    fun = lambda x: x[1] - x[0]
    for k, g in groupby(enumerate(lst), fun):
        l1 = [j for i, j in g]  # 连续数字的列表
        if len(l1) > 1:
            list_all.append(l1)
    return list_all


def product(product_name=None, df_11=None, df_12=None, df_13=None):
    product_all = ("S", "TT", "TQ", "CT", "JI", "VT", "K", "SI", "A", "CAPE", "LI", "TI", "GOES",
                   "DelThetae", "DCI", "MVV", "ICONDI", "KOI", "MDPI", "WMSI", "LCL", "THI")

    function = (S_Index, TT_Index, TQ_Index, CT_index, JI_index, VT_index, K_index, SI_index, A_index, cape,
                LI_index, TI_index, GOES_index, delThetae, DCI_index, MVV_index, Icondi_index, KOI_index,
                MDPI_index, WMSI_index, lcl_index, thi_index)

    index_function = product_all.index(product_name)
    if all([df_11 is not None, df_13 is not None, df_12 is None]):
        result = function[index_function](df_11, df_13)
        return result
    if all([df_11 is not None, df_13 is None, df_12 is None]):
        result = function[index_function](df_11)
        return result


products = ("S", "TT", "TQ", "CT", "JI", "VT", "K", "SI", "A", "CAPE", "LI", "TI", "GOES",
            "DelThetae", "DCI", "MVV", "ICONDI", "KOI", "MDPI", "WMSI", "LCL", "THI", 'PBLH', 'Inv_Hei')

funcs = (S_Index, TT_Index, TQ_Index, CT_index, JI_index, VT_index, K_index, SI_index, A_index, cape,
         LI_index, TI_index, GOES_index, delThetae, DCI_index, MVV_index, Icondi_index, KOI_index,
         MDPI_index, WMSI_index, lcl_index, thi_index, pblh, inver_height_new)


def calc_index(lv2_path, product_name_list, save_path, lv1_path=None, qc=False):
    """
    S:S指数               TT：全总指数              TQ：TQ指数
    CT：交叉总指数         JI：杰弗森指数             VT：垂直总指数
    K：K指数              SI：沙瓦特指数             A：A指数
    CAPE：对流有效位能      LI：抬升指数              TI：Thompson指数
    GOES：混合微下击暴流指数 DelThetae：相当位温差      DCI：深对流指数
    MVV：最大垂直速度       ICONDI:条件性稳定度指数     KOI:KO指数
    MDPI:微下击暴流潜势日指数                         LCL：抬升凝结高度
    WMSI：湿下击暴流严重性指数                        THI：人体舒适度指数
    PBLH:边界层高度        Inv_Hei:逆温层高度


    传入LV2文件路径，和需要生成的指数列表，如["TT", "TQ", "CT"],生成json文件
    """
    try:
        df_lv2 = pd.read_csv(lv2_path, encoding="GBK")
    except:
        df_lv2 = pd.read_csv(lv2_path, encoding="utf-8")
    lv2_11, lv2_12, lv2_13 = get_df(df_lv2, qc)
    lv2_name = os.path.splitext(os.path.basename(lv2_path))[0]

    for i in product_name_list:
        try:
            if i == "VT":
                result = VT_index(lv2_11)

            elif i == 'PBLH':
                result = pblh(lv1_path)

            elif i == 'Inv_Hei':
                result = inver_height_new(lv2_11)

            else:
                index_func = products.index(i)
                result = funcs[index_func](lv2_11, lv2_13)
            # create_json_file(result, os.path.join(save_path, f"{lv2_name}-{i}.json"))
            print(result)
        except Exception as e:
            print(f"{i}指数计算出错：{e}")


if __name__ == '__main__':
    # product_all = ("S", "TT", "TQ", "CT", "JI", "VT", "K", "SI", "A", "CAPE", "LI", "TI", "GOES",
    #                "DelThetae", "DCI", "MVV", "ICONDI", "KOI", "MDPI", "WMSI", "LCL", "THI", "PBLH", "Inv_Hei")
    product_all = ("Inv_Hei",)
    # ret = main(r"54399-AD-2019-05-23LV2.csv", product_all)
    calc_index("54424-AD-2019-07-02LV2.csv", product_all, "./json_file")

