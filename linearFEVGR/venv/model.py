import tensorflow as tf
import numpy as np
import copy
import time
import matplotlib.pyplot as plt

class ServerLinearModel(object):
    def __init__(self, input_dim, output_dim, ag_scalar):
        self.sess = tf.Session()
        with tf.variable_scope("server_param", reuse=tf.AUTO_REUSE):
            self.w = tf.get_variable("local_weights", shape=[input_dim, output_dim],
                                dtype=tf.float32, trainable=True, initializer=tf.truncated_normal_initializer)
            self.b = tf.get_variable("local_bias", shape=[output_dim], trainable=True, dtype=tf.float32)
        self.ag_scalar = ag_scalar
        self.loss = 0
        self.sess.run(tf.global_variables_initializer())
        
    def get_state_dict(self):
        return copy.deepcopy(self.sess.run([self.w, self.b]))
    
    def calculate_average_gradient(self, clients_grads):
        w_grad = tf.zeros(self.w.shape)
        b_grad = tf.zeros(self.b.shape)
        for c in clients_grads:
            w_grad = tf.add(w_grad, c[0])
            b_grad = tf.add(b_grad, c[1])
        w_grad = tf.div(w_grad, len(clients_grads))
        b_grad = tf.div(b_grad, len(clients_grads))
        return self.sess.run([w_grad, b_grad])
    
    def forward(self, data):
        outputs = tf.add(tf.matmul(data, self.w), self.b)
        return self.sess.run(outputs)
    
    def evaluate(self, data, labels):
        self.loss = tf.reduce_mean(self.loss_func(data, labels))
        return self.sess.run(self.loss)
    
    def loss_func(self, predictions, labels):
        loss = tf.reduce_mean(tf.square(tf.subtract(predictions, labels)), axis=0)
        return self.sess.run(loss)

    def update(self, client_states):
        sum_w = tf.zeros(self.w.shape)
        sum_b = tf.zeros(self.b.shape)
        for c in client_states:
            sum_w = tf.add(sum_w, tf.subtract(c[0], self.w))
            sum_b = tf.add(sum_b, tf.subtract(c[1], self.b))
        sum_w = tf.multiply(self.ag_scalar, tf.div(sum_w, len(client_states)))
        sum_b = tf.multiply(self.ag_scalar, tf.div(sum_b, len(client_states)))
        self.w = tf.add(self.w, sum_w)
        self.b = tf.add(self.b, sum_b)
        return copy.deepcopy(self.sess.run([self.w, self.b]))
            

class ClientLinearModel(object):

    def __init__(self, data, batch_size, labels, input_dim, output_dim, uid, global_iter):
        self.sess = tf.Session()
        with tf.variable_scope("param", reuse=tf.AUTO_REUSE):
            self.w = tf.get_variable("local_weights", shape=[input_dim, output_dim],
                                dtype=tf.float32, trainable=True, initializer=tf.truncated_normal_initializer)
            self.b = tf.get_variable("local_bias", shape=[output_dim], trainable=True, dtype=tf.float32)
            self.sw = tf.get_variable("server_weights", [input_dim, output_dim], trainable=False)
            self.sb = tf.get_variable("server_bias", [output_dim], trainable=False)
            self.data = data
            self.labels = labels
            self.batch_size = batch_size
            self.uid = uid
            self.global_iter = global_iter
            self.loss = None
            self.g_loss = None
        self.sess.run(tf.global_variables_initializer())

    def forward(self, batch_data, mode="local"):
        if mode == "local":
            return tf.add(tf.matmul(batch_data, self.w), self.b)
        elif mode == "global":
            return tf.add(tf.matmul(batch_data, self.sw), self.sb)

    def loss_func(self, batch_predictions, batch_labels):
        loss = tf.reduce_mean(tf.square(tf.subtract(batch_predictions, batch_labels)), axis=0)
        return loss

    def calculate_gradient_by_local_weights(self, batch_data, batch_labels):
        batch_predictions = self.forward(batch_data)
        self.loss = self.loss_func(batch_predictions, batch_labels)
        local_gradients = tf.gradients(self.loss, [self.w, self.b])
        return local_gradients, self.loss

    def calculate_gradient_by_server_weights(self, batch_data, batch_labels, server_weights, server_bias):
        self.sw.assign(server_weights)
        self.sb.assign(server_bias)
        batch_predictions = self.forward(batch_data, "global")
        self.g_loss = self.loss_func(batch_predictions, batch_labels)
        global_gradients = tf.gradients(self.g_loss, [self.sw, self.sb])
        return global_gradients
    
    def calculate_gradient_send_server(self, server_weights, server_bias):
        #t0 = time.time()
        data_batches, label_batches = self.split_data_into_batches()
        t1 = time.time()
        #print("user split data time:", str(t1-t0))
        grad_w = []
        grad_b = []
        for i, (x, y) in enumerate(zip(data_batches, label_batches)):
            local_grad = self.calculate_gradient_by_server_weights(x, y, server_weights, server_bias)
            local_grad[0] = tf.squeeze(local_grad[0], axis=0)
            grad_w.append(local_grad[0])
            local_grad[1] = tf.squeeze(local_grad[1], axis=0)
            grad_b.append(local_grad[1])
        t2 = time.time()
        print("enumerate data time:", str(t2-t1))
        avg_w_grad = tf.reduce_mean(grad_w, axis=0)
        avg_b_grad = tf.reduce_mean(grad_b, axis=0)
        return self.sess.run([avg_w_grad, avg_b_grad])

    def split_data_into_batches(self):
        data_batches = []
        label_batches = []
        assert len(self.data) == len(self.labels)
        full_batch_num = len(self.data) // self.batch_size
        for b in range(full_batch_num):
            data_batches.append(self.data[b:b+self.batch_size])
            label_batches.append(self.labels[b:b+self.batch_size])
        if len(self.data) % self.batch_size != 0:
            res = len(self.data) % self.batch_size
            sam_dat = self.data[0:self.batch_size - res]
            sam_lab = self.labels[0:self.batch_size - res]
            # print(sam_lab)
            # print(self.labels[full_batch_num * self.batch_size:])
            dat = np.concatenate((sam_dat, self.data[full_batch_num*self.batch_size:]))
            lab = np.concatenate((sam_lab, self.labels[full_batch_num * self.batch_size:]))
            data_batches.append(dat)
            label_batches.append(lab)
        #print(len(data_batches[0]))
        data_batches = np.array(data_batches, dtype=np.float32)
        #print(data_batches.shape)
        label_batches = np.array(label_batches, dtype=np.float32)
        return data_batches, label_batches

    def fsvgr_fit(self, local_epochs, step_size, lg_scalar, server_gradients, server_weights, server_bias):
        data_batches, label_batches = self.split_data_into_batches()
        #client_loss = []
        self.w.assign(server_weights)
        self.b.assign(server_bias)
        epoch_loss = []
        for epoch in range(local_epochs):
            batch_loss = []
            for i, (x, y) in enumerate(zip(data_batches, label_batches)):
                local_grad, self.loss = self.calculate_gradient_by_local_weights(x, y)
                glob_grad = self.calculate_gradient_by_server_weights(x, y, server_weights, server_bias)
                self.w = tf.subtract(self.w, tf.multiply(step_size, tf.add(tf.multiply(lg_scalar, tf.subtract(local_grad[0], glob_grad[0])), server_gradients[0])))
                self.b = tf.subtract(self.b, tf.multiply(step_size, tf.add(tf.multiply(lg_scalar, tf.subtract(local_grad[1], glob_grad[1])), server_gradients[1])))
                batch_loss.append(self.loss)
            #if epoch % 5 == 0:
            loss_ep = self.sess.run(tf.reduce_mean(batch_loss))
            epoch_loss.append(loss_ep)
            print("global_iter: {}, client: {} local_epoch {} loss: {}".format(self.global_iter, self.uid, epoch, loss_ep))

        plt.figure()
        plt.plot(range(len(epoch_loss)), epoch_loss)
        plt.ylabel('user_train_loss')
        plt.xlabel('num_local_epoch')
        plt.savefig('./save/global_iter_{}_user_{}_train_loss.png'.format(self.global_iter, self.uid))
        return copy.deepcopy(self.sess.run(self.w)), copy.deepcopy(self.sess.run(self.b))







