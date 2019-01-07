import numpy as np
import csv
np.random.seed(0)


def linear_data_generation(path, size, length):
    features = np.random.rand(length, size)
    w = np.random.rand(size)
    b = np.random.rand(1)
    noise = 0.1 * np.random.rand(length)
    labels = np.add(np.add(np.dot(features, np.transpose(w)), b), np.transpose(noise))
    with open(path, "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(np.hstack((features, labels.reshape(length, 1))))
    return features, labels


def linear_iid(data_len, num_users):
    num_items = int(data_len/num_users)
    dict_users, all_idxs = {}, [i for i in range(data_len)]
    for i in range(num_users):
        dict_users[i] = set(np.random.choice(all_idxs, num_items, replace=False))
        all_idxs = list(set(all_idxs) - dict_users[i])
    return dict_users


def collect_data(data, idx):
    user_data = []
    for i in idx:
        user_data.append(data[i])
    return user_data