import tensorflow as tf

from modelFSVRG import ModelFSVRG


IMAGE_SIZE = 28


class ClientModelFSVRG(ModelFSVRG):
    def __init__(self, lr, num_classes):
        self.num_classes = num_classes
        super(ClientModelFSVRG, self).__init__(lr)

    def create_model(self):
        """Model function for CNN."""
        features = tf.placeholder(
            tf.float32, shape=[None, IMAGE_SIZE * IMAGE_SIZE], name='features')
        labels = tf.placeholder(tf.int64, shape=[None], name='labels')
        input_layer = tf.reshape(features, [-1, IMAGE_SIZE, IMAGE_SIZE, 1])
        conv1 = tf.layers.conv2d(
          inputs=input_layer,
          filters=32,
          kernel_size=[5, 5],
          padding="same",
          activation=tf.nn.relu)
        pool1 = tf.layers.max_pooling2d(inputs=conv1, pool_size=[2, 2], strides=2)
        conv2 = tf.layers.conv2d(
            inputs=pool1,
            filters=64,
            kernel_size=[5, 5],
            padding="same",
            activation=tf.nn.relu)
        pool2 = tf.layers.max_pooling2d(inputs=conv2, pool_size=[2, 2], strides=2)
        pool2_flat = tf.reshape(pool2, [-1, 7 * 7 * 64])
        dense = tf.layers.dense(inputs=pool2_flat, units=2048, activation=tf.nn.relu)
        logits = tf.layers.dense(inputs=dense, units=self.num_classes)
        predictions = {
          "classes": tf.argmax(input=logits, axis=1),
          "probabilities": tf.nn.softmax(logits, name="softmax_tensor")
        }
        loss = tf.losses.sparse_softmax_cross_entropy(labels=labels, logits=logits)
        self.network_params = tf.get_collection(tf.GraphKeys.MODEL_VARIABLES)

        # Compute gradient a Loss function
        self.grad = tf.gradients(loss, self.network_params)
        # # Pairing grad, vars
        # grads_vars = [(g, v) for g, v in zip(self.grad, self.network_params)]
        # # Use the computed gradients
        # train_op = tf.apply_gradients(grads_vars)

        train_op = self.optimizer.minimize(
            loss=loss,
            global_step=tf.train.get_global_step())

        eval_metric_ops = tf.count_nonzero(tf.equal(labels, predictions["classes"]))
        return features, labels, train_op, eval_metric_ops


    def fetch_grad(self, w, X, y):
        """
        given weight w_i and data {X_i, y_i}, compute gradient.
        """
        self.network_params.load(w, self.sess)
        return self.sess.run(self.grad, feed_dict={self.features:X, self.labels:y})[0]

    # def fetch_weights(self):
    #     return self.sess.run(self.network_params)