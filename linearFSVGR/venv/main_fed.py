import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import os
import copy
import numpy as np
import tensorflow as tf
from tqdm import tqdm
import tensorboard
from sklearn import metrics
from options import args_parser
from sample_data import linear_iid_data_generation, linear_non_iid_data_generation, \
    get_iris_classification_non_iid_data, linear_iid, collect_data
from model import ServerModel, ClientLocalModel, FederatedClientsModel
import time
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
if __name__ == '__main__':
    
    args = args_parser()
    args.lr = 0.0000001 #0.0001 for regression task. but this local learning rate must be very small to converge for MLP model
    args.ag_scalar = 1 #global step size
    args.task = 'classification' # classification or regression
    args.dataset = 'iris' #iris or ld
    args.num_users = 7
    args.epochs = 10  # numb of global iters
    args.local_ep = 20 # numb of local iters
    args.local_bs = 120  # Local Batch size (120 = full dataset size of a user)
    args.algorithm = 'fsvgr'
    args.lg_scalar = 1.0
    args.iid = False
    args.feature_dims = 4 #change to 4 if implement classification task
    args.hidden_dims = 10 #for classification task only
    args.user_classes = 2 #for classfification task only
    args.output_dims = 3 #change to 3 for classfication task, 1 for regression task
    args.is_reg = True
    args.eval = 'local_iters' #local__iters or users_num

    users_num_list = [3, 8, 12, 20]
    local_iter_list = [5, 10, 15, 20]
    color_list = ['green', 'red', 'blue', 'yellow']
    iter_label_list = ['5 local iters', '10 local iters', '15 local iters', '20 local iters']
    user_label_list = ['5 users', '8 users', '12 users', '20 users']
    plt.figure()
    #plt.title('1-hidden-layer MLP for different num of users')
    for i in range(len(local_iter_list)):
        for k in range(len(users_num_list)):
            print("dataset:", args.dataset, " num_users:", users_num_list[k], " epochs:", args.epochs, "local_ep:", local_iter_list[i])

            train_data, train_label = [], []
            dict_users = {}
            if args.iid and args.task == 'regression':
                train_data, train_label = linear_iid_data_generation('./data/ld_train_iid.csv', args.feature_dims, 100)
                train_label = np.expand_dims(train_label, axis=-1)
                dict_users = linear_iid(len(train_data), users_num_list[k])
            if not args.iid and args.task == 'regression':
                train_data, train_label, dict_users = linear_non_iid_data_generation('./data/ld_train_non_iid.csv', args.feature_dims,
                                                                                     300, users_num_list[k])
                train_label = np.expand_dims(train_label, axis=-1)

            if not args.iid and args.task == 'classification':
                train_data, train_label, dict_users = get_iris_classification_non_iid_data('./data/iris_training.csv',
                                                                                           users_num_list[k], args.user_classes)
                args.is_reg = False

            if args.algorithm == 'fsvgr':
                tf.reset_default_graph()
                server = ServerModel(args.feature_dims, args.hidden_dims, args.output_dims, args.ag_scalar, args.is_reg)
                clients = FederatedClientsModel(dict_users, args.feature_dims, args.hidden_dims, args.output_dims, train_data, train_label,
                                                local_iter_list[i], args.lr, args.lg_scalar, args.local_bs, args.is_reg)

                with tf.name_scope("global_training"):
                    with tf.Session() as sess:
                        log_dir = './log'
                        writer = tf.summary.FileWriter(log_dir+'/global_training', sess.graph)
                        sess.run(tf.global_variables_initializer())
                        global_loss = []
                        global_acc = []
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
                            logits = sess.run(server.forward, feed_dict={server.features: train_data})
                            pred = np.argmax(logits, axis=1)
                            truth = np.argmax(train_label, axis=1)
                            # print(logits)
                            # print(pred)
                            # print(train_label)
                            print("global_test_loss: ", eval_loss)
                            if not args.is_reg:
                                eval_acc = metrics.accuracy_score(y_true=truth, y_pred=pred)
                                print("global_test_acc: ", eval_acc)
                                global_acc.append(eval_acc)
                            global_loss.append(eval_loss)

                        #print(global_loss)
                        #plt.subplot(121)
                        if args.eval == 'local_iters':
                            plt.plot(range(len(global_loss)), global_loss, color=color_list[i], label=iter_label_list[i])
                        if args.eval == 'users_num':
                            plt.plot(range(len(global_loss)), global_loss, color=color_list[k], label=user_label_list[k])

                        # if not args.is_reg:
                        #     plt.subplot(122)
                        #     plt.plot(range(len(global_acc)), global_acc)
                        #     plt.ylabel('global_test_acc')
                        #     plt.xlabel('global_epoch_num')
                        #plt.subplot(121)
                        # plt.plot(range(len(avg_loss)), avg_loss)
                        # plt.ylabel('users_avg_loss')
                        # plt.xlabel('global_epoch_num')
    plt.legend()
    plt.ylabel('global_test_loss')
    plt.xlabel('global_epoch_num')
    plt.savefig('./save/{}_{}_global_and_users_avg_loss_{}_global_epochs_{}_local_epochs.png'.format(args.task, args.iid, args.epochs, args.local_ep))



                
                
    
        
        
    
    
    
        
