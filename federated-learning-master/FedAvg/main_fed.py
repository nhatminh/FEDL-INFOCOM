#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python version: 3.6

import matplotlib
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

from sampling import mnist_iid, mnist_noniid, cifar_iid
from options import args_parser
from Update import LocalFSVGRUpdate
from FedNets import MLP, CNNMnist, CNNCifar
from averaging import average_FSVRG_weights


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
    total_size = np.sum([u[0] for u in users_g])
    for i in range(len(users_g)):
        avg_grad = np.add(avg_grad, users_g[i][0] * users_g[i][1])#+= users_g[i][0] * users_g[i][1]
    avg_grad = np.divide(avg_grad, total_size)#/= total_size
    return avg_grad


if __name__ == '__main__':
    # parse args
    args = args_parser()

    # define paths
    path_project = os.path.abspath('..')

    summary = SummaryWriter('local')
    #Defaults: 100, 10, 10
    args.ag_scalar = 1
    args.model = 'mlp'
    args.dataset = 'mnist'  #  'cifar' or 'mnist'
    args.num_users = 5
    args.frac = 1.          # fraction number of users will be selected to update
    args.epochs = 10        # numb of global iters
    args.local_ep = 100       # numb of local iters
    args.local_bs = 10      # Local Batch size
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
        dataset_train = datasets.CIFAR10('../data/cifar', train=True, transform=transform, target_transform=None, download=True)
        if args.iid:
            dict_users = cifar_iid(dataset_train, args.num_users)
        else:
            exit('Error: only consider IID setting in CIFAR10')
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
            net_glob = MLP(dim_in=len_in, dim_hidden=64, dim_out=args.num_classes).cuda()
        else:
            net_glob = MLP(dim_in=len_in, dim_hidden=64, dim_out=args.num_classes)
    else:
        exit('Error: unrecognized model')
    print(net_glob)
    net_glob.train()

    # copy weights
    w_glob = net_glob.state_dict()

    # training
    global_grad = []
    user_grads = []
    loss_train = []
    cv_loss, cv_acc = [], []
    val_loss_pre, counter = 0, 0
    net_best = None
    val_acc_list, net_list = [], []
    #print(dict_users.keys())
    for iter in tqdm(range(args.epochs)):
        print("=========Global epoch {}=========".format(iter))
        w_locals, loss_locals, acc_locals = [], [], []
        m = max(int(args.frac * args.num_users), 1)
        idxs_users = np.random.choice(range(args.num_users), m, replace=False)
        """
        First communication round: server send w_t to client --> client calculate gradient and send
        to sever --> server calculate average global gradient and send to client
        """
        for idx in idxs_users:
            #print(idx)
            local = LocalFSVGRUpdate(args=args, dataset=dataset_train, idxs=dict_users[idx], tb=summary)
            num_sample, grad_k = local.calculate_global_grad(net=copy.deepcopy(net_glob))
            user_grads.append((num_sample, grad_k))
        global_grad = calculate_avg_grad(user_grads)

        """
        Second communication round: client update w_k_t+1 and send to server --> server update global w_t+1
        """
        for idx in idxs_users:
            print("Training user {}".format(idx))
            local = LocalFSVGRUpdate(args=args, dataset=dataset_train, idxs=dict_users[idx], tb=summary)
            num_samples, w_k, loss, acc = local.update_FSVGR_weights(global_grad, net=copy.deepcopy(net_glob))
            w_locals.append(copy.deepcopy((num_samples, w_k)))
            loss_locals.append(copy.deepcopy(loss))
            acc_locals.append(copy.deepcopy(acc))
        w_glob = average_FSVRG_weights(w_locals, args.ag_scalar, copy.deepcopy(net_glob))

        # copy weight to net_glob
        net_glob.load_state_dict(w_glob)

        # print loss
        loss_avg = sum(loss_locals) / len(loss_locals)
        acc_avg = sum(acc_locals) / len(acc_locals)
        if iter % 1 == 0:
            print('\nTrain loss:', loss_avg)
            print('\nTrain accuracy', acc_avg)
        loss_train.append(loss_avg)
        #acc, _ = test(net_glob)

    # plot loss curve
    plt.figure()
    plt.plot(range(len(loss_train)), loss_train)
    plt.ylabel('train_loss')
    plt.savefig('../save/fed_{}_{}_{}_C{}_iid{}.png'.format(args.dataset, args.model, args.epochs, args.frac, args.iid))

    # testing
    list_acc, list_loss = [], []
    net_glob.eval()
    for c in tqdm(range(args.num_users)):
        net_local = LocalFSVGRUpdate(args=args, dataset=dataset_train, idxs=dict_users[c], tb=summary)
        acc, loss = net_local.test(net=net_glob)
        list_acc.append(acc)
        list_loss.append(loss)
    print("average acc: {:.2f}%".format(100.*sum(list_acc)/len(list_acc)))

