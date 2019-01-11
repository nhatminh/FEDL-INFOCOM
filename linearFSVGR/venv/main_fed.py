import matplotlib.pyplot as plt
import os
import copy
import numpy as np
import tensorflow as tf
from tqdm import tqdm
import tensorboard
from options import args_parser
from sample_data import linear_iid_data_generation, linear_non_iid_data_generation, \
    get_iris_classification_non_iid_data, linear_iid, collect_data
from model import ServerModel, ClientLocalModel, FederatedClientsModel
import time
if __name__ == '__main__':
    
    args = args_parser()
    args.lr = 0.000000001
    args.ag_scalar = 1.0
    args.task = 'classification' #or regression
    args.dataset = 'ld'
    args.num_users = 3
    args.epochs = 20  # numb of global iters
    args.local_ep = 20  # numb of local iters
    args.local_bs = 10  # Local Batch size (420 = full dataset size of a user)
    args.algorithm = 'fsvgr'
    args.lg_scalar = 1.0
    args.iid = False
    args.feature_dims = 4
    args.hidden_dims = 10
    args.user_classes = 2
    args.output_dims = 3
    args.is_reg = True
    print("dataset:", args.dataset, " num_users:", args.num_users, " epochs:", args.epochs, "local_ep:", args.local_ep)
    
    train_data, train_label = [], []
    dict_users = {}
    if args.iid == True and args.task == 'regression':
        train_data, train_label = linear_iid_data_generation('./data/ld_train_iid.csv', args.feature_dims, 100)
        train_label = np.expand_dims(train_label, axis=-1)
        dict_users = linear_iid(len(train_data), args.num_users)
    if not args.iid and args.task == 'regression':
        train_data, train_label, dict_users = linear_non_iid_data_generation('./data/ld_train_non_iid.csv', args.feature_dims,
                                                                             300, args.num_users)
        train_label = np.expand_dims(train_label, axis=-1)
    
    if not args.iid and args.task == 'classification':
        train_data, train_label, dict_users = get_iris_classification_non_iid_data('./data/iris_training.csv', 
                                                                                   args.num_users, args.user_classes)
        args.is_reg = False
        
    if args.algorithm == 'fsvgr':
        tf.reset_default_graph()
        server = ServerModel(args.feature_dims, args.hidden_dims, args.output_dims, args.ag_scalar, args.is_reg)
        clients = FederatedClientsModel(dict_users, args.feature_dims, args.hidden_dims, args.output_dims, train_data, train_label,
                                        args.local_ep, args.lr, args.lg_scalar, args.local_bs, args.is_reg)

        with tf.name_scope("global_training"):
            with tf.Session() as sess:
                log_dir = './log'
                writer = tf.summary.FileWriter(log_dir+'/global_training', sess.graph)
                sess.run(tf.global_variables_initializer())
                global_loss = []
                avg_loss = []
                for g in range(args.epochs):
                    print("================Global Epoch {}================".format(g))
                    sw_1, sb_1, sw_2, sb_2 = sess.run(server.server_params)

                    agg_clients_w1, agg_clients_b1, agg_clients_w2, agg_clients_b2, avg_users_loss = \
                        clients.users_locally_training(sw_1, sb_1, sw_2, sb_2, g)
                    sess.run(server.update, feed_dict={server.clients_agg_weights_1: agg_clients_w1,
                                                       server.clients_agg_bias_1: agg_clients_b1,
                                                       server.clients_agg_weights_2: agg_clients_w2,
                                                       server.clients_agg_bias_2: agg_clients_b2})
                    avg_loss.append(avg_users_loss)
                    #Evaluate server
                    eval_loss = np.mean(sess.run(server.loss, feed_dict={server.features: train_data, server.labels: train_label}))
                    print("global_test_loss: ", eval_loss)
                    global_loss.append(eval_loss)
        plt.figure()
        #print(global_loss)
        plt.subplot(121)
        plt.plot(range(len(global_loss)), global_loss)
        plt.ylabel('global_test_loss')
        plt.xlabel('global_epoch_num')
        plt.subplot(122)
        plt.plot(range(len(avg_loss)), avg_loss)
        plt.ylabel('users_avg_loss')
        plt.xlabel('global_epoch_num')
        plt.savefig('./save/{}_{}_global_and_users_avg_loss_{}_global_epochs_{}_local_epochs.png'.format(args.task, args.iid, args.epochs, args.local_ep))



                
                
    
        
        
    
    
    
        