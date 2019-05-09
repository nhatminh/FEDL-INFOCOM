#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python version: 3.6

import matplotlib
import sys
import pylab
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import copy
import numpy as np
from torchvision import datasets, transforms
from tqdm import tqdm
import torch
import torch.nn.functional as F
from torch import autograd
from tensorboardX import SummaryWriter
import h5py

from sampling import mnist_iid, mnist_noniid, cifar_iid, cifar_noniid
from options import args_parser
from Update import LocalFSVGRUpdate, LocalUpdate, LocalFedProxUpdate
from FedNets import MLP1, CNNMnist, CNNCifar
from averaging import *

def test(net_g, data_loader, args):
    # testing
    test_loss = 0
    correct = 0
    l = len(data_loader)
    for idx, (data, target) in enumerate(data_loader):
        if args.gpu != -1:
            data, target = data.cuda(), target.cuda()
        data, target = autograd.Variable(data), autograd.Variable(target)
        log_probs = net_g(data)
        test_loss += F.nll_loss(log_probs, target, size_average=False).data[0] # sum up batch loss
        y_pred = log_probs.data.max(1, keepdim=True)[1] # get the index of the max log-probability
        correct += y_pred.eq(target.data.view_as(y_pred)).long().cpu().sum()

    test_loss /= len(data_loader.dataset)
    print('\nTest set: Average loss: {:.4f} \nAccuracy: {}/{} ({:.2f}%)\n'.format(
        test_loss, correct, len(data_loader.dataset),
        100. * correct / len(data_loader.dataset)))
    return correct, test_loss


def calculate_avg_grad(users_g):
    avg_grad = np.zeros(users_g[0][1].shape)
    total_size = np.sum([u[0] for u in users_g]) #Total number of samples of all users
    for i in range(len(users_g)):
        avg_grad = np.add(avg_grad, users_g[i][0] * users_g[i][1])#+= users_g[i][0] * users_g[i][1]   <=> n_k * grad_k
    avg_grad = np.divide(avg_grad, total_size)#/= total_size
    # print("avg_grad:", avg_grad)
    return avg_grad

def main(loc_ep, weighted, alg):
    # parse args
    args = args_parser()
    algo_list = ['fedavg', 'fedprox', 'fsvgr']
    # define paths
    path_project = os.path.abspath('..')

    summary = SummaryWriter('local')
    args.gpu = 0  # -1 (CPU only) or GPU = 0
    args.lr = 0.002  # 0.001 for cifar dataset
    args.model = 'mlp'  # 'mlp' or 'cnn'
    args.dataset = 'mnist'  # 'cifar' or 'mnist'
    args.num_users = 5
    args.epochs = 30  # numb of global iters
    args.local_ep = loc_ep  # numb of local iters
    args.local_bs = 1201  # Local Batch size (>=1200 = full dataset size of a user for mnist, 2000 for cifar)
    args.algorithm = alg  # 'fedavg', 'fedprox', 'fsvgr'
    args.iid = False
    args.verbose = False
    print("dataset:", args.dataset, " num_users:", args.num_users, " epochs:", args.epochs, "local_ep:", args.local_ep)

    # load dataset and split users
    dict_users = {}
    dataset_train = []
    if args.dataset == 'mnist':
        dataset_train = datasets.MNIST('../data/mnist/', train=True, download=True,
                                       transform=transforms.Compose([
                                           transforms.ToTensor(),
                                           transforms.Normalize((0.1307,), (0.3081,))
                                       ]))
        # sample users
        if args.iid:
            dict_users = mnist_iid(dataset_train, args.num_users)
        else:
            dict_users = mnist_noniid(dataset_train, args.num_users)
    elif args.dataset == 'cifar':
        transform = transforms.Compose(
            [transforms.ToTensor(),
             transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))])
        dataset_train = datasets.CIFAR10('../data/cifar', train=True, transform=transform, target_transform=None,
                                         download=True)
        if args.iid:
            dict_users = cifar_iid(dataset_train, args.num_users)
        else:
            dict_users = cifar_noniid(dataset_train, args.num_users)
            # exit('Error: only consider IID setting in CIFAR10')
    else:
        exit('Error: unrecognized dataset')
    img_size = dataset_train[0][0].shape

    # build model
    net_glob = None
    if args.model == 'cnn' and args.dataset == 'cifar':
        if args.gpu != -1:
            torch.cuda.set_device(args.gpu)
            net_glob = CNNCifar(args=args).cuda()
        else:
            net_glob = CNNCifar(args=args)
    elif args.model == 'cnn' and args.dataset == 'mnist':
        if args.gpu != -1:
            torch.cuda.set_device(args.gpu)
            net_glob = CNNMnist(args=args).cuda()
        else:
            net_glob = CNNMnist(args=args)
    elif args.model == 'mlp':
        len_in = 1
        for x in img_size:
            len_in *= x
        if args.gpu != -1:
            torch.cuda.set_device(args.gpu)
            # net_glob = MLP1(dim_in=len_in, dim_hidden=128, dim_out=args.num_classes).cuda()
            net_glob = MLP1(dim_in=len_in, dim_hidden=256, dim_out=args.num_classes).cuda()
        else:
            # net_glob = MLP1(dim_in=len_in, dim_hidden=128, dim_out=args.num_classes)
            net_glob = MLP1(dim_in=len_in, dim_hidden=256, dim_out=args.num_classes)
    else:
        exit('Error: unrecognized model')
    print("Nerual Net:", net_glob)

    net_glob.train()  # Train() does not change the weight values
    # copy weights
    w_glob = net_glob.state_dict()

    # w_size = 0
    # for k in w_glob.keys():
    #     size = w_glob[k].size()
    #     if (len(size) == 1):
    #         nelements = size[0]
    #     else:
    #         nelements = size[0] * size[1]
    #     w_size += nelements * 4
    #     # print("Size ", k, ": ",nelements*4)
    # print("Weight Size:", w_size, " bytes")
    # print("Weight & Grad Size:", w_size * 2, " bytes")
    # print("Each user Training size:", 784 * 8 / 8 * args.local_bs, " bytes")
    # print("Total Training size:", 784 * 8 / 8 * 60000, " bytes")
    # # training
    global_grad = []
    user_grads = []
    loss_test = []
    acc_test = []
    cv_loss, cv_acc = [], []
    val_loss_pre, counter = 0, 0
    net_best = None
    val_acc_list, net_list = [], []
    # print(dict_users.keys())

    rs_avg_acc, rs_avg_loss, rs_glob_acc, rs_glob_loss= [], [], [], []

    ###  FedAvg Aglorithm  ###
    if args.algorithm == 'fedavg':
        # for iter in tqdm(range(args.epochs)):
        for iter in range(args.epochs):
            w_locals, loss_locals, acc_locals, num_samples_list = [], [], [], []
            for idx in range(args.num_users):
                if(args.local_bs>1200):
                    local = LocalUpdate(args=args, dataset=dataset_train, idxs=dict_users[idx], tb=summary, bs=200*(4+idx)) #Batch_size bs = full data
                else:
                    local = LocalUpdate(args=args, dataset=dataset_train, idxs=dict_users[idx], tb=summary, bs=args.local_bs)
                num_samples, w, loss, acc = local.update_weights(net=copy.deepcopy(net_glob))
                num_samples_list.append(num_samples)
                w_locals.append(copy.deepcopy(w))
                loss_locals.append(copy.deepcopy(loss))
                # print("User ", idx, " Acc:", acc, " Loss:", loss)
                acc_locals.append(copy.deepcopy(acc))
            # update global weights
            if(weighted):
                w_glob = weighted_average_weights(w_locals, num_samples_list)
            else:
                w_glob = average_weights(w_locals)


            # copy weight to net_glob
            net_glob.load_state_dict(w_glob)
            # global test
            list_acc, list_loss = [], []
            net_glob.eval()
            # for c in tqdm(range(args.num_users)):
            for c in range(args.num_users):
                if (args.local_bs > 1200):
                    net_local = LocalUpdate(args=args, dataset=dataset_train, idxs=dict_users[c], tb=summary, bs=200*(4+c)) #Batch_size bs = full data
                else:
                    net_local = LocalUpdate(args=args, dataset=dataset_train, idxs=dict_users[c], tb=summary, bs=args.local_bs)
                acc, loss = net_local.test(net=net_glob)
                list_acc.append(acc)
                list_loss.append(loss)
            print("\nEpoch: {}, Global test loss {}, Global test acc: {:.2f}%".format(iter,
                                                                                      sum(list_loss) / len(list_loss),
                                                                                      100. * sum(list_acc) / len(
                                                                                          list_acc)))

            # print loss
            loss_avg = sum(loss_locals) / len(loss_locals)
            acc_avg = sum(acc_locals) / len(acc_locals)
            if args.epochs % 1 == 0:
                print('\nUsers Train Average loss:', loss_avg)
                print('\nTrain Train Average accuracy', acc_avg)
            # loss_test.append(sum(list_loss) / len(list_loss))
            # acc_test.append(sum(list_acc) / len(list_acc))

            rs_avg_acc.append(acc_avg)
            rs_avg_loss.append(loss_avg)
            rs_glob_acc.append(sum(list_acc) / len(list_acc))
            rs_glob_loss.append(sum(list_loss) / len(list_loss))
            # if (acc_avg >= 0.89):
            #     return iter+1

    ###  FedProx Aglorithm  ###
    elif args.algorithm == 'fedprox':
        args.mu = 0.005  ### change mu 0.001
        args.limit = 0.3
        # for iter in tqdm(range(args.epochs)):
        for iter in range(args.epochs):
            w_locals, loss_locals, acc_locals, num_samples_list = [], [], [], []
            # m = max(int(args.frac * args.num_users), 1)
            # idxs_users = np.random.choice(range(args.num_users), m, replace=False)
            for idx in range(args.num_users):
                if(args.local_bs>1200):
                    local = LocalFedProxUpdate(args=args, dataset=dataset_train, idxs=dict_users[idx], tb=summary, bs=200*(4+idx)) #Batch_size bs = full data
                else:
                    local = LocalFedProxUpdate(args=args, dataset=dataset_train, idxs=dict_users[idx], tb=summary, bs=args.local_bs)

                num_samples, w, loss, acc = local.update_FedProx_weights(net=copy.deepcopy(net_glob))
                num_samples_list.append(num_samples)
                w_locals.append(copy.deepcopy(w))
                loss_locals.append(copy.deepcopy(loss))
                # print("User ", idx, " Acc:", acc, " Loss:", loss)
                acc_locals.append(copy.deepcopy(acc))
            # update global weights
            if(weighted):
                w_glob = weighted_average_weights(w_locals, num_samples_list)
            else:
                w_glob = average_weights(w_locals)

            # copy weight to net_glob
            net_glob.load_state_dict(w_glob)
            # global test
            list_acc, list_loss = [], []
            net_glob.eval()
            # for c in tqdm(range(args.num_users)):
            for c in range(args.num_users):
                if (args.local_bs > 1200):
                    net_local = LocalFedProxUpdate(args=args, dataset=dataset_train, idxs=dict_users[c], tb=summary,
                                            bs=200 * (4 + c))  # Batch_size bs = full data
                else:
                    net_local = LocalFedProxUpdate(args=args, dataset=dataset_train, idxs=dict_users[c], tb=summary,
                                            bs=args.local_bs)

                acc, loss = net_local.test(net=net_glob)
                list_acc.append(acc)
                list_loss.append(loss)
            print("\nEpoch: {}, Global test loss {}, Global test acc: {:.2f}%".format(iter,
                                                                                      sum(list_loss) / len(list_loss),
                                                                                      100. * sum(list_acc) / len(
                                                                                          list_acc)))

            # print loss
            loss_avg = sum(loss_locals) / len(loss_locals)
            acc_avg = sum(acc_locals) / len(acc_locals)
            if args.epochs % 1 == 0:
                print('\nUsers train average loss:', loss_avg)
                print('\nUsers train average accuracy', acc_avg)
            # loss_test.append(sum(list_loss) / len(list_loss))
            # acc_test.append(sum(list_acc) / len(list_acc))

            rs_avg_acc.append(acc_avg)
            rs_avg_loss.append(loss_avg)
            rs_glob_acc.append(sum(list_acc) / len(list_acc))
            rs_glob_loss.append(sum(list_loss) / len(list_loss))
            # if (acc_avg >= 0.89):
            #     return iter+1

    ###  FSVGR Aglorithm  ###
    elif args.algorithm == 'fsvgr':
        args.ag_scalar = 1.  # 0.001 or 0.1
        args.lg_scalar = 1
        args.threshold = 0.001
        # for iter in tqdm(range(args.epochs)):
        for iter in range(args.epochs):

            print("=========Global epoch {}=========".format(iter))
            w_locals, loss_locals, acc_locals = [], [], []
            """
            First communication round: server send w_t to client --> client calculate gradient and send
            to sever --> server calculate average global gradient and send to client
            """
            for idx in range(args.num_users):
                local = LocalFSVGRUpdate(args=args, dataset=dataset_train, idxs=dict_users[idx], tb=summary)
                num_sample, grad_k = local.calculate_global_grad(net=copy.deepcopy(net_glob))
                user_grads.append([num_sample, grad_k])
            global_grad = calculate_avg_grad(user_grads)

            """
            Second communication round: client update w_k_t+1 and send to server --> server update global w_t+1
            """
            for idx in range(args.num_users):
                print("Training user {}".format(idx))
                local = LocalFSVGRUpdate(args=args, dataset=dataset_train, idxs=dict_users[idx], tb=summary)
                num_samples, w_k, loss, acc = local.update_FSVGR_weights(global_grad, idx, copy.deepcopy(net_glob),
                                                                         iter)
                w_locals.append(copy.deepcopy([num_samples, w_k]))
                print("Global_Epoch ", iter, "User ", idx, " Acc:", acc, " Loss:", loss)
                loss_locals.append(copy.deepcopy(loss))
                acc_locals.append(copy.deepcopy(acc))

            # w_t = net_glob.state_dict()
            w_glob = average_FSVRG_weights(w_locals, args.ag_scalar, copy.deepcopy(net_glob), args.gpu)

            # copy weight to net_glob
            net_glob.load_state_dict(w_glob)

            # global test
            list_acc, list_loss = [], []
            net_glob.eval()
            # for c in tqdm(range(args.num_users)):
            for c in range(args.num_users):
                net_local = LocalFSVGRUpdate(args=args, dataset=dataset_train, idxs=dict_users[c], tb=summary)
                acc, loss = net_local.test(net=net_glob)
                list_acc.append(acc)
                list_loss.append(loss)

            print("\nTest Global Weight:", list_acc)
            print("\nEpoch: {}, Global test loss {}, Global test acc: {:.2f}%".format(iter,
                                                                                      sum(list_loss) / len(list_loss),
                                                                                      100. * sum(list_acc) / len(
                                                                                          list_acc)))

            # print loss
            loss_avg = sum(loss_locals) / len(loss_locals)
            acc_avg = sum(acc_locals) / len(acc_locals)
            if iter % 1 == 0:
                print('\nEpoch: {}, Users train average loss: {}'.format(iter, loss_avg))
                print('\nEpoch: {}, Users train average accuracy: {}'.format(iter, acc_avg))
            loss_test.append(sum(list_loss) / len(list_loss))
            acc_test.append(sum(list_acc) / len(list_acc))

            if (acc_avg >= 0.89):
                return
    if(weighted):
        alg=alg+'1'
    simple_save_data(loc_ep, alg, rs_avg_acc, rs_avg_loss, rs_glob_acc, rs_glob_loss)
    plot_rs(loc_ep, alg)

    # # plot loss curve
    # plt.figure(1)
    # plt.subplot(121)
    # plt.plot(range(len(loss_test)), loss_test)
    # plt.ylabel('train_loss')
    # plt.xlabel('num_epoches')
    # plt.subplot(122)
    # plt.plot(range(len(acc_test)), acc_test)
    # plt.ylabel('train_accuracy')
    # plt.xlabel('num_epoches')
    # plt.savefig(
    #     '../save/new_fed_{}_{}_{}_{}_C{}_iid{}.png'.format(args.algorithm, args.dataset, args.model, args.epochs,
    #                                                        args.frac, args.iid))


def simple_save_data(loc_ep,alg,rs_avg_acc, rs_avg_loss, rs_glob_acc, rs_glob_loss):
    with h5py.File('data_{}_{}.h5'.format(alg,loc_ep), 'w') as hf:
        hf.create_dataset('rs_avg_acc', data=rs_avg_acc)
        hf.create_dataset('rs_avg_loss', data=rs_avg_loss)
        hf.create_dataset('rs_glob_acc', data=rs_glob_acc)
        hf.create_dataset('rs_glob_loss', data=rs_glob_loss)

        hf.close()

def simple_read_data(loc_ep,alg):
    hf = h5py.File('data_{}_{}.h5'.format(alg,loc_ep), 'r')
    rs_avg_acc = np.array(hf.get('rs_avg_acc')[:])
    rs_avg_loss = np.array(hf.get('rs_avg_loss')[:])
    rs_glob_acc = np.array(hf.get('rs_glob_acc')[:])
    rs_glob_loss = np.array(hf.get('rs_glob_loss')[:])
    return rs_avg_acc, rs_avg_loss, rs_glob_acc, rs_glob_loss

def plot_rs(loc_ep,alg):
    rs_avg_acc, rs_avg_loss, rs_glob_acc, rs_glob_loss = simple_read_data(loc_ep,alg)

    plt.figure(1)
    plt.plot(rs_avg_acc)
    plt.ylabel('avg_acc')
    plt.xlabel('num_epochs')
    plt.savefig('avg_acc_{}_{}.png'.format(alg,loc_ep))

    plt.figure(2)
    plt.plot(rs_avg_loss)
    plt.ylabel('avg_loss')
    plt.xlabel('num_epochs')
    plt.savefig('avg_loss_{}_{}.png'.format(alg,loc_ep))

    plt.figure(3)
    plt.plot(rs_glob_acc)
    plt.ylabel('glob_acc')
    plt.xlabel('num_epochs')
    plt.savefig('glob_acc_{}_{}.png'.format(alg, loc_ep))

    plt.figure(4)
    plt.plot(rs_glob_loss)
    plt.ylabel('glob_loss')
    plt.xlabel('num_epochs')
    plt.savefig('glob_loss_{}_{}.png'.format(alg, loc_ep))

def plot_summary():

    Numb_Glob_Iters=20
    Numb_Algs = 2

    avg_acc, avg_acc1     = np.zeros((3,30)), np.zeros((3,30))
    avg_loss, avg_loss1   = np.zeros((3, 30)), np.zeros((3, 30))
    glob_acc, glob_acc1   = np.zeros((3, 30)), np.zeros((3, 30))
    glob_loss, glob_loss1 = np.zeros((3, 30)), np.zeros((3, 30))


    avg_acc[0, :], avg_loss[0, :], glob_acc[0, :], glob_loss[0, :] = simple_read_data(20, 'fedavg')
    avg_acc[1, :], avg_loss[1, :], glob_acc[1, :], glob_loss[1, :] = simple_read_data(20, 'fedprox')
    avg_acc[2, :], avg_loss[2, :], glob_acc[2, :], glob_loss[2, :] = simple_read_data(20, 'fedprox1')

    avg_acc1[0, :], avg_loss1[0, :], glob_acc1[0, :], glob_loss1[0, :] = simple_read_data(30, 'fedavg')
    avg_acc1[1, :], avg_loss1[1, :], glob_acc1[1, :], glob_loss1[1, :] = simple_read_data(30, 'fedprox')
    avg_acc1[2, :], avg_loss1[2, :], glob_acc1[2, :], glob_loss1[2, :] = simple_read_data(30, 'fedprox1')
    algs_lbl = ["FedAvg - 20", "FedProx - 20","FedProx1 - 20"]
    algs_lbl1 = ["FedAvg - 30", "FedProx - 30", "FedProx1 - 30"]

    plt.figure(1)
    for i in range(Numb_Algs):
        plt.plot(avg_acc[i, :Numb_Glob_Iters], linestyle=":", label=algs_lbl[i])
        plt.plot(avg_acc1[i, :Numb_Glob_Iters], label=algs_lbl1[i])
    plt.legend(loc='best')
    plt.ylabel('Average Accuracy')
    plt.xlabel('Number of Global Iterations')
    plt.savefig('avg_acc.png')

    plt.figure(2)
    for i in range(Numb_Algs):
        plt.plot(avg_loss[i, :Numb_Glob_Iters], linestyle=":", label=algs_lbl[i])
        plt.plot(avg_loss1[i, :Numb_Glob_Iters], label=algs_lbl1[i])
    plt.legend(loc='best')
    plt.ylabel('Average Loss')
    plt.xlabel('Number of Global Iterations')
    plt.savefig('avg_loss.png')

    plt.figure(3)
    for i in range(Numb_Algs):
        plt.plot(glob_acc[i, :Numb_Glob_Iters], linestyle=":", label=algs_lbl[i])
        plt.plot(glob_acc1[i, :Numb_Glob_Iters], label=algs_lbl1[i])
    plt.legend(loc='best')
    plt.ylabel('Global Accuracy')
    plt.xlabel('num_epochs')
    plt.savefig('glob_acc.png')

    plt.figure(4)
    for i in range(Numb_Algs):
        plt.plot(glob_loss[i, :Numb_Glob_Iters], linestyle=":", label=algs_lbl[i])
        plt.plot(glob_loss1[i, :Numb_Glob_Iters], label=algs_lbl1[i])
    plt.legend(loc='best')
    plt.ylabel('Global Loss')
    plt.xlabel('Number of Global Iterations')
    plt.savefig('glob_loss.png')



if __name__ == '__main__':
   # local_eps = range(5, 45, 5)
   # numb_of_glob_iters = []
   # for k in local_eps:
   #     print("\n== TESTING with local_eps=", k)
   #     numb_of_glob_iters.append(main(k))
   # print("-- FINISH --")
   # print(numb_of_glob_iters)

    SUMARRY = True
    if(SUMARRY):
        plot_summary()
    else:
       main(loc_ep=10, weighted=False, alg='fedprox') #'fedavg', 'fedprox', 'fsvgr'
       print("-- FINISH -- :",)


# # testing
# list_acc, list_loss = [], []
# net_glob.eval()
# for c in tqdm(range(args.num_users)):
#     net_local = LocalFSVGRUpdate(args=args, dataset=dataset_train, idxs=dict_users[c], tb=summary)
#     acc, loss = net_local.test(net=net_glob)
#     list_acc.append(acc)
#     list_loss.append(loss)
# print("average acc: {:.2f}%".format(100.*sum(list_acc)/len(list_acc)))
#
