import matplotlib.pyplot as plt
import os
import copy
import numpy as np
import tensorflow as tf
from tqdm import tqdm
import tensorboard
from options import args_parser
from sample_data import linear_data_generation, linear_iid, collect_data
from model import ServerModel, ClientLocalModel, FederatedClientsModel
import time
if __name__ == '__main__':
    
    args = args_parser()
    args.gpu = -1  # -1 (CPU only) or GPU = 0
    args.lr = 0.0000035
    args.ag_scalar = 1.0
    args.model = 'linear'
    args.dataset = 'ld'
    args.num_users = 3
    args.epochs = 20  # numb of global iters
    args.local_ep = 20  # numb of local iters
    args.local_bs = 10  # Local Batch size (420 = full dataset size of a user)
    args.algorithm = 'fsvgr'
    args.lg_scalar = 1.0
    args.iid = True
    args.feature_dims = 10
    print("dataset:", args.dataset, " num_users:", args.num_users, " epochs:", args.epochs, "local_ep:", args.local_ep)
    
    train_data, train_label = [], []
    dict_users = {}
    if args.dataset == 'ld' and args.model == 'linear':
        train_data, train_label = linear_data_generation('./data/ld_train.csv', args.feature_dims, 100)
        train_label = np.expand_dims(train_label, axis=-1)
    if args.iid == True and args.model == 'linear':
        dict_users = linear_iid(len(train_data), args.num_users)
        
    if args.algorithm == 'fsvgr':
        tf.reset_default_graph()
        server = ServerModel(args.feature_dims, 1, args.ag_scalar)
        clients = FederatedClientsModel(dict_users, args.feature_dims, 1, train_data, train_label,
                                        args.local_ep, args.lr, args.lg_scalar, args.local_bs)

        with tf.name_scope("global_training"):
            with tf.Session() as sess:
                log_dir = './log'
                writer = tf.summary.FileWriter(log_dir+'/global_training', sess.graph)
                sess.run(tf.global_variables_initializer())
                global_loss = []
                avg_loss = []
                for g in range(args.epochs):
                    print("================Global Epoch {}================".format(g))
                    sw, sb = sess.run(server.server_params)
                    agg_clients_w, agg_clients_b, avg_users_loss = clients.users_locally_training(sw, sb, g)
                    sess.run(server.update, feed_dict={server.clients_agg_weights: agg_clients_w,
                                                       server.clients_agg_bias: agg_clients_b})
                    avg_loss.append(avg_users_loss)
                    #Evaluate server
                    eval_loss = sess.run(server.loss, feed_dict={server.features: train_data, server.labels: train_label})
                    global_loss.append(eval_loss)
        plt.figure()
        #print(global_loss)
        plt.subplot(121)
        plt.plot(range(len(global_loss)), global_loss)
        plt.ylabel('global_test_loss')
        plt.xlabel('global_epoch_num')
        plt.subplot(122)
        plt.plot(range(len(avg_loss)), avg_loss)
        plt.savefig('./save/global_and_users_avg_loss_{}_global_epochs_{}_local_epochs.png'.format(args.epochs, args.local_ep))



                
                
    
        
        
    
    
    
        