from itertools import groupby

microwave_height = (0, 10, 25, 50, 75, 100, 130, 160, 190, 220, 250, 280, 310, 340, 370, 400, 430, 460, 490, 520,
                    560, 600, 640, 680, 720, 760, 800, 840, 880, 920, 960, 1000, 1040, 1080, 1120, 1160, 1200, 1260,
                    1320, 1380, 1440, 1500, 1560, 1620, 1680, 1740, 1800, 1890, 1980, 2170, 2260, 2350, 2430, 2500, 2600,
                    2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3650, 3800, 3950, 4100, 4250, 4400, 4550, 4600,
                    4800, 5000, 5200, 5400, 5600, 5800, 6000, 6300, 6600, 6900, 7200, 7500, 7800, 8100, 8400, 8700, 9000,
                    9300, 9600, 9800, 10000)


def inver_height_new(LV2_11):
    """逆温层(新)"""
    result = {"unit": "m",
              "label": "逆温层，利用微波辐射计LV2温度数据进行计算，分别取温度随高度增加而递增的高度起始点和终止点"}
    data = []
    height = list(microwave_height)
    for i in LV2_11:
        time = i[1]
        data_temp = i[11:]
        meta = {
            'x': time,
            'y': []
        }

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
            if diff_val > 1 and height[min(i)] < 1000:
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