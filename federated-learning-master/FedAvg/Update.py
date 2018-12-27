#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python version: 3.6

import torch
from torch import nn, autograd
from torch.utils.data import DataLoader, Dataset
import numpy as np
from sklearn import metrics
import copy
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt


class DatasetSplit(Dataset):
    def __init__(self, dataset, idxs):
        self.dataset = dataset
        self.idxs = list(idxs)

    def __len__(self):
        return len(self.idxs)

    def __getitem__(self, item):

        image, label = self.dataset[int(self.idxs[item])]
        return image, label


class LocalUpdate(object):
    def __init__(self, args, dataset, idxs, tb):
        self.args = args
        self.loss_func = nn.CrossEntropyLoss() #nn.NLLLoss()
        self.ldr_train, self.ldr_val, self.ldr_test = self.train_val_test(dataset, list(idxs))
        self.tb = tb

    def train_val_test(self, dataset, idxs):
        # split train, validation, and test
        idxs_train = idxs#[:420]
        idxs_val = idxs[420:480]
        idxs_test = idxs[480:]

        #If batch_size = 420 -> random permutation of all items
        train = DataLoader(DatasetSplit(dataset, idxs_train), batch_size=self.args.local_bs, shuffle=True)
        val = DataLoader(DatasetSplit(dataset, idxs_val), batch_size=int(len(idxs_val)/10), shuffle=True)
        test = DataLoader(DatasetSplit(dataset, idxs_train), batch_size=int(len(idxs_train)), shuffle=False)
        return train, val, test

    def update_weights(self, net):
        net.train()
        # train and update
        optimizer = torch.optim.SGD(net.parameters(), lr=self.args.lr, momentum=0.5)

        epoch_loss = []
        epoch_acc = []
        for iter in range(self.args.local_ep):
            batch_loss = []
            for batch_idx, (images, labels) in enumerate(self.ldr_train):
                if self.args.gpu != -1:
                    images, labels = images.cuda(), labels.cuda()
                images, labels = autograd.Variable(images), autograd.Variable(labels)
                net.zero_grad()
                log_probs = net(images)
                loss = self.loss_func(log_probs, labels)
                loss.backward()
                optimizer.step()
                if self.args.gpu != -1:
                    loss = loss.cpu()
                if self.args.verbose and batch_idx % 10 == 0:
                    print('Update Epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(
                        iter, batch_idx * len(images), len(self.ldr_train.dataset),
                               100. * batch_idx / len(self.ldr_train), loss.data.item()))
                self.tb.add_scalar('loss', loss.data.item())
                batch_loss.append(loss.data.item())
            epoch_loss.append(sum(batch_loss)/len(batch_loss))
            acc, _ = self.test(net)
            epoch_acc.append(acc)
        return net.state_dict(), sum(epoch_loss) / len(epoch_loss), sum(epoch_acc) / len(epoch_acc)

    def test(self, net):
        loss = 0
        log_probs = []
        labels = []
        for batch_idx, (images, labels) in enumerate(self.ldr_test):
            if self.args.gpu != -1:
                images, labels = images.cuda(), labels.cuda()
            images, labels = autograd.Variable(images), autograd.Variable(labels)
            net = net.float()
            log_probs = net(images)
            loss = self.loss_func(log_probs, labels)
        if self.args.gpu != -1:
            loss = loss.cpu()
            log_probs = log_probs.cpu()
            labels = labels.cpu()
        y_pred = np.argmax(log_probs.data, axis=1)
        acc = metrics.accuracy_score(y_true=labels.data, y_pred=y_pred)
        return acc, loss.data.item()


class LocalFSVGRUpdate(LocalUpdate):
    def __init__(self, args, dataset, idxs, tb):
        super(LocalFSVGRUpdate, self).__init__(args, dataset, idxs, tb)
        self.lg_scalar = args.lg_scalar
        self.lr = args.lr

    def calculate_global_grad(self, net):
        global_grad = [np.zeros(v.shape) for v in net.parameters()]
        total_size = 0
        for batch_idx, (images, labels) in enumerate(self.ldr_train):
            if self.args.gpu != -1:
                images, labels = images.cuda(), labels.cuda()
            total_size += len(images)
            images, labels = autograd.Variable(images), autograd.Variable(labels)
            net.zero_grad()
            log_probs = net(images)
            loss = self.loss_func(log_probs, labels)
            loss.backward()

            for i, param in enumerate(net.parameters()):
                # print("Parameter:",i, " size: ", len(param.grad.data.numpy()))
                grad_data = param.grad.data
                if self.args.gpu != -1:
                    grad_data = grad_data.cpu()
                global_grad[i] += grad_data.numpy()[0]
            if self.args.gpu != -1:
                loss = loss.cpu()
        #global_grad = np.divide(global_grad, total_size)#global_grad /= total_size => Check ?? (we don't need to divide size)
        global_grad = np.divide(global_grad, 1)
        return total_size, global_grad

    def fetch_grad(self, net, images, labels):
        grad = [np.zeros(v.shape) for v in net.parameters()]
        if self.args.gpu != -1:
            images, labels = images.cuda(), labels.cuda()
        images, labels = autograd.Variable(images), autograd.Variable(labels)
        net.zero_grad()
        log_probs = net(images)
        loss = self.loss_func(log_probs, labels)
        loss.backward()
        for i, param in enumerate(net.parameters()):
            if self.args.gpu != -1:
                grad[i] = param.grad.data.cpu()
            else: grad[i] = param.grad.data
        return grad

    def update_FSVGR_weights(self, server_avg_grad, uid, net):
        net.train()
        # train and update
        epoch_loss = []
        epoch_acc = []
        server_net = copy.deepcopy(net)
        total_size = 0
        for iter in range(self.args.local_ep):
            batch_loss = []
            for batch_idx, (images, labels) in enumerate(self.ldr_train):
                if self.args.gpu != -1:
                    images, labels = images.cuda(), labels.cuda()
                if iter == 0:
                    total_size += len(images)
                images, labels = autograd.Variable(images), autograd.Variable(labels)
                net = net.float()
                net.zero_grad()
                log_probs = net(images)
                loss = self.loss_func(log_probs, labels)
                loss.backward()
                client_w_grad = self.fetch_grad(net, images, labels)
                server_w_grad = self.fetch_grad(server_net, images, labels)
                # if iter == 50:
                #     print(client_w_grad)
                for i, param in enumerate(net.parameters()):
                    weights = param.data
                    if self.args.gpu != -1:
                        weights = weights.cpu()
                        param.data = np.subtract(weights, self.lr * (
                            np.add(self.lg_scalar * (np.subtract(client_w_grad[i], server_w_grad[i])),
                                   server_avg_grad[i]))).cuda()
                    else:
                        param.data = np.subtract(weights, self.lr * (
                            np.add(self.lg_scalar * (np.subtract(client_w_grad[i], server_w_grad[i])),
                                   server_avg_grad[i])))

                    #param.data -= self.lr * (self.lg_scalar * (np.subtract(client_w_grad[i], server_w_grad[i])) + server_avg_grad[i])
                if self.args.gpu != -1:
                    loss = loss.cpu()
                if self.args.verbose and batch_idx % 10 == 0:
                    print('Update Epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(
                        iter, batch_idx * len(images), len(self.ldr_train.dataset),
                              100. * batch_idx / len(self.ldr_train), loss.data.item()))
                self.tb.add_scalar('loss', loss.data.item())
                batch_loss.append(loss.data.item())
            epoch_loss.append(sum(batch_loss) / len(batch_loss))
            acc, _ = self.test(net)
            epoch_acc.append(acc)
        plt.figure()
        plt.subplot(211)
        plt.plot(range(len(epoch_loss)), epoch_loss)
        plt.ylabel('train_loss')
        plt.xlabel('num_local_epoches')
        plt.subplot(212)
        plt.plot(range(len(epoch_acc)), epoch_acc)
        plt.ylabel('train_accuracy')
        plt.xlabel('num_local_epoches')
        plt.savefig('../save/{}.png'.format(uid))

            #print('Local Epoch: {}, accuracy: {:.6f}'.format(iter, acc))
        return total_size, net.state_dict(), sum(epoch_loss) / len(epoch_loss), sum(epoch_acc) / len(epoch_acc)

