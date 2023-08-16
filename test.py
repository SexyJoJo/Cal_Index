import matplotlib.pyplot as plt
import pandas as pd
from indices import *
import indices
from get_press import get_press
import numpy as np

# df = pd.read_csv(
#     r"C:\Users\JOJO\Documents\WeChat Files\wxid_o318fo1t9tk622\FileStorage\File\2023-07\54511_2023-04-28.csv")
# results = []
# for index, row in df.iterrows():
#     height = np.array(eval(row['height'])) * 1000
#     tem = eval(row['hTamb'])
#     rhu = eval(row['hRh'])
#     dewp = []
#     prs = []
#     for i in range(len(tem)):
#         dewp.append(cal_dewp(tem[i], rhu[i]))
#         prs.append(get_press(height[i], tem[i], 0, 1013.25, tem[0]))
#
#     parcel_prof = cal_parcel_prof(prs, tem, dewp)
#     results.append(round(LI_index(prs, tem, parcel_prof), 2))
# print(results)
# plt.plot(results)
# plt.show()

def get_press_by_height(height):
    return 1013.25 * math.pow((1 - height / 44331.0), (1.0 / 0.1903))


# temperature = [31.9, 32.01, 32.02, 31.95729035074671, 31.68, 31.52, 31.43, 31.23, 31.1, 31.11, 30.97, 30.75, 30.59, 30.49, 30.11, 29.82, 29.64, 29.33, 29.06, 28.87, 28.76, 28.33, 28.02, 27.48, 26.93, 26.42, 26.09, 25.57, 25.24, 25.02, 24.96664621790144, 24.66, 24.29, 24.09, 23.8, 23.21, 22.8, 22.38, 21.9, 21.37, 21.201578353166234, 21.05, 20.53, 20.16, 19.74, 19.36, 19.12, 18.87, 18.46, 18.16, 17.95, 17.42, 16.49, 15.39, 13.93, 12.2, 10.541172705368655, 10.41, 8.76, 7.07, 5.46, 3.91, 2.9139239472601726, 2.5, 1.05, 0.0, -0.39, -1.88, -3.43, -4.0, -5.15, -5.459144835572779, -6.88, -8.66, -9.0, -10.0, -10.33, -12.08, -13.8, -15.56, -16.653985184054175, -17.31, -19.14, -20.0, -21.03, -22.91, -24.71, -26.6, -28.54, -30.4, -31.43306320633803, -32.31, -34.21, -35.69, -36.79]

temperature = [20.5, 20.84, 20.74, 20.53, 20.42, 20.29, 20.16, 19.96, 19.84, 19.68, 19.48, 19.33, 19.18, 19.01, 18.81, 18.65, 18.47, 18.28, 18.07, 17.86, 17.59, 17.3, 16.99, 16.63, 16.28, 15.9, 15.49, 15.06, 14.68, 14.29, 13.95, 13.63, 13.3, 12.96, 12.61, 12.24, 11.9, 11.57, 11.22, 10.92, 10.56, 10.18, 9.79, 9.41, 8.99, 8.63, 8.28, 7.95, 7.65, 7.1, 6.32, 5.36, 4.25, 2.86, 1.52, 0.21, -1.12, -2.49, -3.84, -5.24, -6.67, -8.11, -9.6, -11.07, -12.6, -14.13, -15.72, -17.32, -19.01, -20.68, -22.43, -24.18, -25.99, -27.79, -29.61, -31.49, -33.43, -35.35, -37.35, -39.38, -41.37, -42.95, -44.18]
rhu = [80.3, 85.55, 87.02, 88.6, 88.8, 89.08, 89.44, 89.9, 90.31, 90.21, 90.0, 89.84, 89.51, 88.83, 88.52, 88.3, 87.89, 87.63, 88.1, 88.11, 88.18, 88.43, 88.43, 87.93, 87.78, 87.96, 87.87, 88.2, 88.53, 88.61, 88.1, 87.36, 86.49, 85.46, 84.06, 83.28, 82.38, 80.9, 79.69, 78.62, 76.92, 75.13, 74.01, 72.5, 71.3, 70.14, 69.68, 68.96, 68.45, 68.58, 70.93, 75.31, 81.78, 89.43, 95.0, 95.0, 95.0, 95.0, 88.76, 81.41, 74.19, 66.41, 60.21, 54.43, 49.08, 44.54, 41.06, 37.28, 34.34, 32.32, 29.94, 27.85, 26.46, 24.33, 22.26, 20.11, 18.12, 15.59, 13.56, 11.4, 9.48, 7.74, 6.65]
height = np.array([0.0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0]) * 1000 + 2380.9
dew_point = []
pressure = []
for i in range(len(temperature)):
    dew_point.append(cal_dewp(temperature[i], rhu[i]))
    pressure.append(get_press(height[i], temperature[i], 0, 1013.25, temperature[0]))
    # pressure.append(get_press_by_height(height[i]))

# 杰弗森指数
JI = indices.JI_Index(pressure, temperature, dew_point)
print('JI', JI)

# A指数
A = indices.A_index(pressure, temperature, dew_point)
print('A', A)

# S指数
S = indices.S_index(pressure, temperature, dew_point)
print('S', S)

# 条件性稳定度指数
IC = indices.IC_index(pressure, temperature, dew_point)
print('IC', IC)

# 深对流指数
DCI = indices.DCI_index(pressure, temperature, dew_point)
print('DCI', DCI)

# 自由对流高度
LFC = indices.LFC_index(pressure, temperature, dew_point)
print('LFC', LFC)

# # 强天气威胁指数
# SWEAT = indices.SWEAT_index(pressure, temperature, dew_point, w_speed, w_direct)
# print('SWEAT', SWEAT)
#
# # 山崎指数
# KYI = indices.KYI_index(pressure, temperature, dew_point, w_speed, w_direct, lat=100)   # lat为站点经度
# print("KYI", KYI)

# 逆温层
inv_hei = indices.inver_height(height, temperature)
print("inv_hei", inv_hei)

# 露点温度
dewp = indices.cal_dewp(15, 50)
print("dewp", dewp)

# 高度转气压
# 参数:(待求气压对应高度， 待求气压对应温度， 站点海拔高度， 站点近地面气压， 站点近地面温度)
p = get_press(100, 20, 0, 1013.25, 22)
print("prs", p)

# 抬升凝结气压，抬升凝结温度
lcl_pressure, lcl_temperature = indices.LCL_index(pressure, temperature, dew_point)
print('lcl_pressure', lcl_pressure)
print('lcl_temperature', lcl_temperature)

# 状态曲线
parcel_prof = indices.cal_parcel_prof(pressure, temperature, dew_point)
print("parcel_prof", parcel_prof)

# 全总指数
TT = indices.TT_Index(pressure, temperature, dew_point)
print('TT', TT)

# 沙氏指数
SI = indices.SI_index(pressure, temperature, dew_point)
print('SI', SI)

# K指数
K = indices.K_index(pressure, temperature, dew_point)
print('K', K)

# 抬升指数
LI = indices.LI_index(pressure, temperature, parcel_prof)
print('LI', LI)

# 对流有效位能
CAPE, CIN = indices.CAPE_index(pressure, temperature, dew_point)
print('CAPE', CAPE)
print('CIN', CIN)