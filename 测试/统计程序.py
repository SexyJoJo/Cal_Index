import math
import os
import sympy
from matplotlib import pyplot as plt
import numpy as np
import parse_000
import indices
import pandas as pd
PATH = r"D:\文档\微波辐射计\北京指数\3月7-14日TlogP"


def get_stations(time):
    stations = []
    filename = str(time) + '.000'
    with open(os.path.join(PATH, filename), 'r', encoding='gbk') as f:
        for line in f:
            if line.startswith('5'):
                stations.append(int(line[:5]))
    stations.sort()
    result = [str(station) for station in stations]
    return result

def get_micaps_indices(column):
    # df = pd.read_table(r"data\Micaps_out\Tph_21070220.dat", sep=r'\s+', skiprows=1)
    df = pd.read_table(r"../data/Micaps_out/Tph_22031420.dat", sep=r'\s+', skiprows=1)
    return df[column].tolist()

def get_MSE(A, B):
    if len(A) == len(B):
        return sum([(x - y) ** 2 for x, y in zip(A, B)]) / len(A)

def get_RMSE(A, B):
    mse = get_MSE(A, B)
    rmse = math.sqrt(mse)
    return round(rmse, 2)

def R(A, B):
    return 1 - get_MSE(A, B) / np.var(np.array(A))


if __name__ == '__main__':
    time = 20210702200000
    # time = 22031420
    # time = 22030708
    stations = get_stations(time)

    result = []
    for station in stations:
        parse_000.get_data(station, time)
        df = parse_000.to_df()
        with open("temp.txt", "r") as f:
            meta = f.readline().split()    # 读取表头

        # CAPE_value = indices.CAPE(
        #     df['pressure'].tolist(),
        #     df['temperature'].tolist(),
        #     df['dew_point'].tolist(),
        # )
        # result.append(CAPE_value)

        # TT = indices.TT_Index(
        #     df['pressure'].tolist(),
        #     df['temperature'].tolist(),
        #     df['dew_point'].tolist(),
        # )
        # result.append(TT)

        # JI = indices.JI_Index(
        #     df['pressure'].tolist(),
        #     df['temperature'].tolist(),
        #     df['dew_point'].tolist(),
        # )
        # result.append(JI)

        # A = indices.A_index(
        #     df['pressure'].tolist(),
        #     df['temperature'].tolist(),
        #     df['dew_point'].tolist(),
        # )
        # result.append(A)

        # K = indices.K_index(
        #     df['pressure'].tolist(),
        #     df['temperature'].tolist(),
        #     df['dew_point'].tolist(),
        # )
        # result.append(K)

        # SI = indices.SI_index(
        #     df['pressure'].tolist(),
        #     df['temperature'].tolist(),
        #     df['dew_point'].tolist(),
        # )
        # result.append(SI)

        # LCL = indices.LCL_index(
        #     df['pressure'].tolist(),
        #     df['temperature'].tolist(),
        #     df['dew_point'].tolist(),
        # )
        # result.append(LCL)

        # IC = indices.IC_index(
        #     df['pressure'].tolist(),
        #     df['temperature'].tolist(),
        #     df['dew_point'].tolist(),
        # )
        # result.append(IC)
        #
        # BIC = indices.BIC_index(
        #     df['pressure'].tolist(),
        #     df['temperature'].tolist(),
        #     df['dew_point'].tolist(),
        # )
        # result.append(BIC)
        #
        # LI = indices.LI_index(
        #     df['pressure'].tolist(),
        #     df['temperature'].tolist(),
        #     df['dew_point'].tolist(),
        # )
        # result.append(LI)
        #
        # DCI = indices.DCI_index(
        #     df['pressure'].tolist(),
        #     df['temperature'].tolist(),
        #     df['dew_point'].tolist(),
        # )
        # result.append(DCI)

        SWEAT = indices.SWEAT_index(
            df['pressure'].tolist(),
            df['temperature'].tolist(),
            df['dew_point'].tolist(),
            df['w_speed'].tolist(),
            df['w_direct'].tolist(),
        )
        result.append(SWEAT)

    #     KYI = indices.KYI_index(
    #         df['pressure'].tolist(),
    #         df['temperature'].tolist(),
    #         df['dew_point'].tolist(),
    #         df['w_speed'].tolist(),
    #         df['w_direct'].tolist(),
    #         lat=float(meta[2])
    #     )
    #     result.append(KYI)
    # print("KYI", result)

    # micaps_result = get_micaps_indices("湿对流有效位能")
    # micaps_result = get_micaps_indices("总指数")
    # micaps_result = get_micaps_indices("修正杰弗逊指数")
    # micaps_result = get_micaps_indices("A指数")
    # micaps_result = get_micaps_indices("K指数")
    # micaps_result = get_micaps_indices("沙氏指数")
    # micaps_result = get_micaps_indices("抬升凝结高度")
    # micaps_result = get_micaps_indices("对流稳定度指数")
    # micaps_result = get_micaps_indices("最大对流稳定度指数")
    # micaps_result = get_micaps_indices("抬升指数")
    # micaps_result = get_micaps_indices("修正对流指数")
    micaps_result = get_micaps_indices("强天气威胁指数")

    # 去除异常值
    drop_index = []
    if len(result) == len(micaps_result):
        for i in range(len(micaps_result)):
            if micaps_result[i] > 9000:
                drop_index.append(i)
                continue
            try:
                if result[i] is None:
                    drop_index.append(i)
            except TypeError:
                if isinstance(result[i], sympy.core.numbers.NaN):
                    drop_index.append(i)

        for i in drop_index[::-1]:
            del result[i]
            del micaps_result[i]
            del stations[i]

        drop_index = []
        for i in range(len(result)):
            if result[i] > 9000 or result[i] < -9000:
                drop_index.append(i)
        for i in drop_index[::-1]:
            del result[i]
            del micaps_result[i]
            del stations[i]

        # 绘图
        RMSE = get_RMSE(result, micaps_result)  # 均方根误差
        COR = R(result, micaps_result)  # 相关系数

        print(micaps_result)
        print(result)

        # for i in range(len(micaps_result)):
        #     if abs(micaps_result[i] - result[i]) > 20:
        #         print(stations[i])

        plt.rcParams['font.sans-serif'] = ['SimHei']
        plt.rcParams['axes.unicode_minus'] = False
        plt.scatter(micaps_result, result)
        plt.xlabel("Micaps")
        plt.ylabel("指数程序")

        # plt.title("CAPE")
        # plt.xlim(0, 2000)
        # plt.ylim(0, 2000)
        # plt.plot([0, 2000], [0, 2000], color='black')
        # plt.text(60, 1850, f'COR={COR}')
        # plt.text(60, 1750, f'RMSE={RMSE}')
        # plt.text(60, 1650, f'站点数={len(result)}')

        # plt.title("总指数TT")
        # plt.xlim(0, 100)
        # plt.ylim(0, 100)
        # plt.plot([0, 100], [0, 100], color='black')
        # plt.text(5, 70, f'站点数={len(result)}')
        # plt.text(5, 80, f'RMSE={RMSE}')
        # plt.text(5, 90, f'COR={COR}')

        # plt.title("修正杰弗逊指数JI")
        # plt.xlim(0, 200)
        # plt.ylim(0, 200)
        # plt.plot([0, 200], [0, 200], color='black')
        # plt.text(5, 125, f'站点数={len(result)}')
        # plt.text(5, 150, f'RMSE={RMSE}')
        # plt.text(5, 175, f'COR={COR}')

        # plt.title("A指数")
        # plt.xlim(0, 50)
        # plt.ylim(0, 50)
        # plt.plot([0, 50], [0, 50], color='black')
        # plt.text(1, 35, f'站点数={len(result)}')
        # plt.text(1, 40, f'RMSE={RMSE}')
        # plt.text(1, 45, f'COR={COR}')

        # plt.title("K指数")
        # plt.xlim(0, 50)
        # plt.ylim(0, 50)
        # plt.plot([0, 50], [0, 50], color='black')
        # plt.text(1, 35, f'站点数={len(result)}')
        # plt.text(1, 40, f'RMSE={RMSE}')
        # plt.text(1, 45, f'COR={COR}')

        # plt.title("沙氏指数SI")
        # plt.xlim(-10, 50)
        # plt.ylim(-10, 50)
        # plt.plot([-10, 50], [-10, 50], color='black')
        # plt.text(-8, 35, f'站点数={len(result)}')
        # plt.text(-8, 40, f'RMSE={RMSE}')
        # plt.text(-8, 45, f'COR={COR}')

        # plt.title("抬升凝结高度LCL")
        # plt.xlim(0, 1200)
        # plt.ylim(0, 1200)
        # plt.plot([0, 1200], [0, 1200], color='black')
        # plt.text(10, 900, f'站点数={len(result)}')
        # plt.text(10, 1000, f'RMSE={RMSE}')
        # plt.text(10, 1100, f'COR={COR}')

        # plt.title("对流稳定度指数IC")
        # plt.xlim(-50, 50)
        # plt.ylim(-50, 50)
        # plt.plot([-50, 50], [-50, 50], color='black')
        # plt.text(-48, 35, f'站点数={len(result)}')
        # plt.text(-48, 40, f'RMSE={RMSE}')
        # plt.text(-48, 45, f'COR={COR}')

        # plt.title("最大对流稳定度指数BIC")
        # plt.xlim(-50, 100)
        # plt.ylim(-50, 100)
        # plt.plot([-50, 100], [-50, 100], color='black')
        # plt.text(-48, 35, f'站点数={len(result)}')
        # plt.text(-48, 40, f'RMSE={RMSE}')
        # plt.text(-48, 45, f'COR={COR}')

        # plt.title("抬升指数LI")
        # plt.xlim(-50, 100)
        # plt.ylim(-50, 100)
        # plt.plot([-50, 100], [-50, 100], color='black')
        # plt.text(-47, 80, f'站点数={len(result)}')
        # plt.text(-47, 87, f'RMSE={RMSE}')
        # plt.text(-47, 94, f'COR={COR}')

        # plt.title("深对流指数DCI")
        # plt.xlim(-50, 100)
        # plt.ylim(-50, 100)
        # plt.plot([-50, 100], [-50, 100], color='black')
        # plt.text(-47, 80, f'站点数={len(result)}')
        # plt.text(-47, 87, f'RMSE={RMSE}')
        # plt.text(-47, 94, f'COR={COR}')

        plt.title("SWEAT")
        plt.xlim(0, 500)
        plt.ylim(0, 500)
        plt.plot([0, 500], [0, 500], color='black')
        plt.text(60, 480, f'COR={COR}')
        plt.text(60, 460, f'RMSE={RMSE}')
        plt.text(60, 440, f'站点数={len(result)}')
        plt.show()


