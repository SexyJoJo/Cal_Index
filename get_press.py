import math


def get_press(z2, t2, z1, p1, t1):
    """
    给定高度求出压强
    :param z2: 高度
    :param t2: z2对应温度
    :param z1: 站台海拔
    :param p1: 站台近地面气压
    :param t1: 站台近地面温度
    :return: z2对应气压
    """
    a = 1 / 273
    t = (t1 + t2) / 2
    p2 = p1 / (math.e ** ((z2 - z1) / (8000 * 1 + a * t)))
    return round(p2, 3)