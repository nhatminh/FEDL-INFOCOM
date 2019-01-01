import matplotlib.pyplot as plt
import os
import copy
import numpy as np
import tensorflow as tf
from tqdm import tqdm
import tensorboard
from options import args_parser
from sample_data import linear_data_generation, linear_iid, collect_data
from model import ServerLinearModel, ClientLinearModel
import time
if __name__ == '__main__':
    
    args = args_parser()
    args.gpu = -1  # -1 (CPU only) or GPU = 0
    args.lr = 0.0000035
    args.ag_scalar = 1.0
    args.model = 'linear'  # 'mlp' or 'cnn' or 'linear'
    args.dataset = 'ld'  # 'cifar' or 'mnist' or 'ld'
    args.num_users = 3
    args.frac = 1.  # fraction number of users will be selected to update
    args.epochs = 2  # numb of global iters
    args.local_ep = 20  # numb of local iters
    args.local_bs = 10  # Local Batch size (420 = full dataset size of a user)
    args.algorithm = 'fsvgr'
    args.lg_scalar = 1.0
    args.iid = True
    print("dataset:", args.dataset, " num_users:", args.num_users, " epochs:", args.epochs, "local_ep:", args.local_ep)
    
    train_data, train_label = [], []
    dict_users = {}
    if args.dataset == 'ld' and args.model == 'linear':
        train_data, train_label = linear_data_generation('./data/ld_train.csv', 1, 100)

    if args.iid == True and args.model == 'linear':
        dict_users = linear_iid(len(train_data), args.num_users)
        
    if args.algorithm == 'fsvgr':
        server = ServerLinearModel(1, 1, args.ag_scalar)
        global_loss = []
        for g in range(args.epochs):
            print("================Global Epoch {}================".format(g))
            clients_grads = []
            clients_states = []
            """
            First communication round: server send w_t to client --> client calculate gradient and send
            to sever --> server calculate average global gradient and send to client
            """
            server_state = server.get_state_dict()
            #t0 = time.time()
            for u in range(args.num_users):
                #t1 = time.time()
                user_data = collect_data(train_data, dict_users[u])
                user_label = collect_data(train_label, dict_users[u])
                #t2 = time.time()
                #print("user data collection time:", str(t2-t1))
                client_u = ClientLinearModel(user_data, args.local_bs, user_label, 1, 1, u, g)
                #t3 = time.time()
                #print("user model initialization time:", str(t3-t2))
                grad_w, grad_b = client_u.calculate_gradient_send_server(server_state[0], server_state[1])
                #t4 = time.time()
                #print("calculate gradient time:", str(t4-t3))
                client_u_grad = [grad_w, grad_b]
                clients_grads.append(client_u_grad)
            #t5 = time.time()
            #print("first communication time:", str(t5-t0))
            global_w_grad, global_b_grad = server.calculate_average_gradient(clients_grads)
            #t6 = time.time()
            #print("global average gradient calculation time:", str(t6-t5))
            """
            Second communication round: client update w_k_t+1 and send to server --> server update global w_t+1
            """
            for u in range(args.num_users):
                user_data = collect_data(train_data, dict_users[u])
                user_label = collect_data(train_label, dict_users[u])
                client_u = ClientLinearModel(user_data, args.local_bs, user_label, 1, 1, u, g)
                t7 = time.time()
                client_w, client_b = client_u.fsvgr_fit(args.local_ep, args.lr, args.lg_scalar, [global_w_grad, global_b_grad], server_state[0], server_state[1])
                t8 = time.time()
                print("user fit model time:", str(t8-t7))
                client_state = [client_w, client_b]
                clients_states.append(client_state)
            server.update(clients_states)

            #Evaluate server model
            test_data, test_label = linear_data_generation('./data/ld_test.csv', 1, 100)
            test_data = np.array(test_data, dtype=np.float32)
            predictions = server.forward(test_data)
            loss = server.evaluate(predictions, test_label)
            print("Global epoch: {}, server loss: {}".format(g, loss))
            global_loss.append(loss)
            
        plt.figure()
        print(global_loss)
        plt.plot(range(len(global_loss)), global_loss)
        plt.ylabel('global_test_loss')
        plt.xlabel('global_epoch_num')
        plt.savefig('./save/global_test_loss_{}_global_epochs_{}_local_epochs.png'.format(args.epochs, args.local_ep))



                
                
    
        
        
    
    
    
        