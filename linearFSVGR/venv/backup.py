import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import copy
from sample_data import collect_data


class ServerModel(object):
    def __init__(self, input_dim, output_dim, ag_scalar):
        with tf.name_scope("server_model"):
            self.features = tf.placeholder(dtype=tf.float32, name='server_features', shape=[None, input_dim])
            self.labels = tf.placeholder(dtype=tf.float32, name='server_labels', shape=[None, output_dim])
            with tf.variable_scope("server_params", reuse=tf.AUTO_REUSE):
                self.w = tf.get_variable("W", shape=[input_dim, output_dim],
                                         dtype=tf.float32, initializer=tf.truncated_normal_initializer)
                self.b = tf.get_variable("b", shape=output_dim, dtype=tf.float32,
                                         initializer=tf.truncated_normal_initializer)
            with tf.variable_scope("clients_pass_in_params", reuse=tf.AUTO_REUSE):
                self.clients_agg_weights = tf.placeholder(dtype=tf.float32, shape=(input_dim, output_dim), name='Ws')
                self.clients_agg_bias = tf.placeholder(dtype=tf.float32, shape=output_dim, name='bs')
            self.ag_scalar = ag_scalar
            self._loss = None
            self._forward = None
            self._update = None

    @property
    def server_params(self):
        return self.w, self.b

    @property
    def forward(self):
        with tf.name_scope("server_forward"):
            if self._forward is None:
                self._forward = tf.add(tf.matmul(self.features, self.w), self.b)
            return self._forward

    @property
    def loss(self):
        with tf.name_scope("server_loss_calc"):
            if self._loss is None:
                self._loss = tf.reduce_mean(tf.square(tf.subtract(self.forward, self.labels)), axis=0)
            return self._loss

    @property
    def update(self):
        with tf.name_scope("server_param_update"):
            if self._update is None:
                self._update = [self.w.assign(tf.add(self.w, tf.multiply(self.clients_agg_weights, self.ag_scalar))),
                                self.b.assign(tf.add(self.b, tf.multiply(self.clients_agg_bias, self.ag_scalar)))]
            return self._update


class ClientLocalModel(object):
    def __init__(self, input_dim, output_dim, step_size, lg_scalar):
        with tf.name_scope("client_model"):
            self.input_dim = input_dim
            self.output_dim = output_dim
            self.step_size = step_size
            self.lg_scalar = lg_scalar
            self.features = tf.placeholder(dtype=tf.float32, name='client_features', shape=[None, input_dim])
            self.labels = tf.placeholder(dtype=tf.float32, name='client_labels', shape=[None, output_dim])
            with tf.variable_scope('passed_in_params', reuse=tf.AUTO_REUSE):
                self.server_weights = tf.placeholder(dtype=tf.float32, shape=(self.input_dim, self.output_dim),
                                                     name='W')
                self.server_bias = tf.placeholder(dtype=tf.float32, shape=self.output_dim, name='b')
                self.sw_grad = tf.placeholder(dtype=tf.float32, shape=(self.input_dim, self.output_dim), name='W_grad')
                self.sb_grad = tf.placeholder(dtype=tf.float32, shape=self.output_dim, name='b_grad')
            with tf.variable_scope("client_params", reuse=tf.AUTO_REUSE):
                self.w = tf.get_variable('W', shape=(input_dim, output_dim),
                                         initializer=tf.truncated_normal_initializer)
                self.b = tf.get_variable('b', shape=output_dim, initializer=tf.truncated_normal_initializer)
            with tf.variable_scope("global_params_on_local", reuse=tf.AUTO_REUSE):
                self.sw = tf.get_variable("W", shape=(input_dim, output_dim), trainable=False)
                self.sb = tf.get_variable("b", shape=output_dim, trainable=False)
            self._local_forward = None
            self._global_forward = None
            self._local_loss = None
            self._global_loss = None
            self._local_grad = None
            self._global_grad = None
            self._assign_server = None
            self._update_local = None
            self._assign_local = None
            # self._update_local_b = None

    @property
    def local_params(self):
        return self.w, self.b

    @property
    def assign_server(self):
        with tf.name_scope("assign_passed_in_param"):
            if self._assign_server is None:
                self._assign_server = [self.sw.assign(self.server_weights), self.sb.assign(self.server_bias)]
            return self._assign_server

    @property
    def assign_local(self):
        with tf.name_scope("update_local_by_received_params"):
            if self._assign_local is None:
                self._assign_local = [self.w.assign(self.server_weights), self.b.assign(self.server_bias)]
            return self._assign_local

    @property
    def local_param_forward(self):
        with tf.name_scope("local_param_forward"):
            if self._local_forward is None:
                self._local_forward = tf.add(tf.matmul(self.features, self.w), self.b)
            return self._local_forward

    @property
    def global_param_forward(self):
        with tf.name_scope('global_param_forward'):
            if self._global_forward is None:
                self._global_forward = tf.add(tf.matmul(self.features, self.sw), self.sb)
            return self._global_forward

    @property
    def local_loss(self):
        with tf.name_scope('local_param_loss_calc'):
            if self._local_loss is None:
                self._local_loss = tf.reduce_mean(tf.square(tf.subtract(self.local_param_forward, self.labels)), axis=0)
            return self._local_loss

    @property
    def global_loss(self):
        with tf.name_scope('global_param_loss_calc'):
            if self._global_loss is None:
                self._global_loss = tf.reduce_mean(tf.square(tf.subtract(self.global_param_forward, self.labels)),
                                                   axis=0)
            return self._global_loss

    @property
    def calculate_gradient_by_local_params(self):
        with tf.name_scope('local_param_grad_calc'):
            if self._local_grad is None:
                self._local_grad = tf.gradients(self.local_loss, [self.w, self.b])
            return self._local_grad

    @property
    def calculate_gradient_by_server_params(self):
        with tf.name_scope('global_param_grad_calc'):
            if self._global_grad is None:
                self._global_grad = tf.gradients(self.global_loss, [self.sw, self.sb])
            return self._global_grad

    @property
    def optimize_local(self):
        with tf.name_scope('local_optimize'):
            if self._update_local is None:
                self._update_local = [self.w.assign(tf.subtract(self.w, tf.multiply(self.step_size, tf.add(
                    tf.multiply(self.lg_scalar, tf.subtract(self.calculate_gradient_by_local_params[0],
                                                            self.calculate_gradient_by_server_params[0])),
                    self.sw_grad)))),
                                      self.b.assign(tf.subtract(self.b, tf.multiply(self.step_size, tf.add(
                                          tf.multiply(self.lg_scalar,
                                                      tf.subtract(self.calculate_gradient_by_local_params[1],
                                                                  self.calculate_gradient_by_server_params[1])),
                                          self.sb_grad))))]
            return self._update_local


class FederatedClientsModel(object):
    def __init__(self, dict_users, input_dim, output_dim, full_data, full_label, local_epochs, step_size, lg_scalar,
                 batch_size):
        self.dict_users = dict_users
        self.full_data = full_data
        self.full_label = full_label
        self.input_dim = input_dim
        self.output_dim = output_dim
        self.local_epochs = local_epochs
        self.step_size = step_size
        self.lg_scalar = lg_scalar
        self.batch_size = batch_size
        self.users, self.users_batches_data, self.users_batches_label = self.create_users_and_assign_data()

    def create_users_and_assign_data(self):
        users = []
        users_data = {}
        users_label = {}
        users_batches_data = {}
        users_batches_label = {}
        for u in range(len(self.dict_users)):
            users.append(ClientLocalModel(self.input_dim, self.output_dim, self.step_size, self.lg_scalar))
            users_data[u] = collect_data(self.full_data, self.dict_users[u])
            users_label[u] = collect_data(self.full_label, self.dict_users[u])
            users_batches_data[u], users_batches_label[u] = self.split_data_into_batches(users_data[u], users_label[u])
        return users, users_batches_data, users_batches_label

    def split_data_into_batches(self, data, label):
        data_batches = []
        label_batches = []
        assert len(data) == len(label)
        full_batch_num = len(data) // self.batch_size
        for b in range(full_batch_num):
            data_batches.append(data[b * self.batch_size:b * self.batch_size + self.batch_size])
            label_batches.append(label[b * self.batch_size:b * self.batch_size + self.batch_size])
        if len(data) % self.batch_size != 0:
            res = len(data) % self.batch_size
            sam_dat = data[0:self.batch_size - res]
            sam_lab = label[0:self.batch_size - res]
            dat = np.concatenate((sam_dat, data[full_batch_num * self.batch_size:]))
            lab = np.concatenate((sam_lab, label[full_batch_num * self.batch_size:]))
            data_batches.append(dat)
            label_batches.append(lab)
        data_batches = np.array(data_batches, dtype=np.float32)
        label_batches = np.array(label_batches, dtype=np.float32)
        return data_batches, label_batches

    def calc_server_avg_gradients(self, server_weights, server_bias):
        with tf.name_scope("calculate_and_send_server_gradients"):
            users_w_gradients = []
            users_b_gradients = []
            grad_calcs = []
            for user in self.users:
                grad_calcs.append(user.calculate_gradient_by_server_params)
            with tf.Session() as sess:
                sess.run(tf.global_variables_initializer())
                for i in range(len(self.users)):
                    user_w_grad = []
                    user_b_grad = []
                    for idx, (x, y) in enumerate(zip(self.users_batches_data[i], self.users_batches_label[i])):
                        sess.run(self.users[i].assign_server, feed_dict={self.users[i].server_weights: server_weights,
                                                                         self.users[i].server_bias: server_bias})
                        w_grad, b_grad = sess.run(grad_calcs[i], feed_dict={self.users[i].features: x,
                                                                            self.users[i].labels: y,
                                                                            self.users[
                                                                                i].server_weights: server_weights,
                                                                            self.users[i].server_bias: server_bias})
                        user_w_grad.append(w_grad)
                        user_b_grad.append(b_grad)
                    user_data_ratio = len(self.dict_users[i]) / len(self.full_data)
                    users_w_gradients.append(user_data_ratio * np.sum(user_w_grad, axis=0))
                    users_b_gradients.append(user_data_ratio * np.sum(user_b_grad, axis=0))
                if np.sum(users_w_gradients, axis=0).shape == ():
                    avg_server_w_grad = np.expand_dims(np.expand_dims(np.array(np.sum(users_w_gradients, axis=0)), -1),
                                                       -1)
                    avg_server_b_grad = np.expand_dims(np.array(np.sum(users_b_gradients, axis=0)), -1)
                else:
                    avg_server_w_grad = np.array(np.sum(users_w_gradients, axis=0))
                    avg_server_b_grad = np.array(np.sum(users_b_gradients, axis=0))
            return avg_server_w_grad, avg_server_b_grad

    def users_locally_training(self, server_weights, server_bias, global_iter):
        with tf.name_scope('accept_server_average_grads'):
            sw_grad, sb_grad = self.calc_server_avg_gradients(server_weights, server_bias)
        with tf.name_scope('federated_training_users_locally'):
            # Initialize optimizer and loss
            optimizers = []
            loss_calcs = []
            for user in self.users:
                optimizers.append(user.optimize_local)
                loss_calcs.append(user.local_loss)
            # start federated training users locally for 1 global epoch
            with tf.Session() as sess:
                log_dir = './log'
                writer = tf.summary.FileWriter(log_dir + '/local_training', sess.graph)
                sess.run(tf.global_variables_initializer())
                users_loss = []
                agg_w_params = np.zeros_like(server_weights, dtype=np.float32)
                agg_b_params = np.zeros_like(server_bias, dtype=np.float32)
                for i in range(len(self.users)):
                    user_local_epoch_loss = []
                    sess.run(self.users[i].assign_local, feed_dict={self.users[i].server_weights: server_weights,
                                                                    self.users[i].server_bias: server_bias})
                    for e in range(self.local_epochs):
                        batch_loss = []
                        for idx, (x, y) in enumerate(zip(self.users_batches_data[i], self.users_batches_label[i])):
                            _, loss = sess.run([optimizers[i], loss_calcs[i]], feed_dict={self.users[i].features: x,
                                                                                          self.users[i].labels: y,
                                                                                          self.users[
                                                                                              i].server_weights: server_weights,
                                                                                          self.users[
                                                                                              i].server_bias: server_bias,
                                                                                          self.users[
                                                                                              i].sw_grad: sw_grad,
                                                                                          self.users[
                                                                                              i].sb_grad: sb_grad})
                            batch_loss.append(loss)
                        print('global_iter: {}, user: {}, local_epoch: {}, loss: {}'.format(global_iter, i, e,
                                                                                            np.mean(batch_loss)))
                        user_local_epoch_loss.append(np.mean(batch_loss))
                    users_loss.append(np.mean(user_local_epoch_loss))
                    user_w_params, user_b_params = sess.run(self.users[i].local_params)
                    user_data_ratio = len(self.dict_users[i]) / len(self.full_data)
                    agg_w_params = np.add(agg_w_params, user_data_ratio * np.subtract(user_w_params, server_weights))
                    agg_b_params = np.add(agg_b_params, user_data_ratio * np.subtract(user_b_params, server_bias))
                avg_users_loss = np.mean(users_loss)
                print('global_iter: {}, users average loss: {}'.format(global_iter, avg_users_loss))

            return agg_w_params, agg_b_params, avg_users_loss





















