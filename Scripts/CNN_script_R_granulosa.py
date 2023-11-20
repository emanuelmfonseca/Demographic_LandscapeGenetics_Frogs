import tensorflow as tf
import keras
from keras.models import Sequential
from keras.layers import Conv2D
from keras.layers import MaxPooling2D
from keras.layers import Flatten
from keras.layers import Dense
from keras.layers import Dropout
from keras.layers import AveragePooling2D

import numpy as np

np.random.seed(123)
tf.random.set_seed(1234)

classifier = Sequential()

classifier.add(Conv2D(1, (3,1), input_shape = (422, 160, 1), activation = 'relu'))
classifier.add(MaxPooling2D(pool_size = (3,1)))

classifier.add(Flatten())
classifier.add(Dense(units = 40, activation = 'relu',kernel_initializer='he_uniform'))
classifier.add(Dense(units = 3, activation = 'softmax'))

classifier.compile(optimizer = "adam", loss = 'categorical_crossentropy', metrics = ['accuracy'])

# Part 2 - Fitting the CNN to the images

from keras.preprocessing.image import ImageDataGenerator

train_datagen = ImageDataGenerator(rescale = 1./255)

test_datagen = ImageDataGenerator(rescale = 1./255)

training_set = train_datagen.flow_from_directory('./Simulated_datasets/training_dataset',
                                                 target_size = (422,160),
                                                 batch_size = 100,
                                                 class_mode = 'categorical',
                                                #class_mode = 'binary',
                                                 color_mode='grayscale',
                                                 shuffle=True)

test_set = test_datagen.flow_from_directory('./Simulated_datasets/test_dataset',
                                            target_size = (422, 160),
                                            batch_size = 100,
                                            class_mode = 'categorical',
                                            #class_mode = 'binary',
                                            color_mode='grayscale',
                                            shuffle=False)

classifier.fit(training_set,
                        steps_per_epoch = int(6000/100),
                         epochs = 1,
                         validation_data = test_set,
                         validation_steps = int(1500/100))


import pandas as pd
import seaborn as sn
import matplotlib.pyplot as plt
import numpy as np

import sklearn.metrics as metrics
from sklearn.metrics import confusion_matrix as cm

batch_size = 100

test_set2=test_set

test_set2.reset()
Y_pred = classifier.predict_generator(test_set2, steps = 1500 // batch_size)
y_pred = np.argmax(Y_pred, axis=-1)

#pd.DataFrame(Y_pred).to_csv('/fs/project/PAA0202/Emanuel/Second_chapter/CNN_calibration.csv')

matrix = cm(test_set2.classes, y_pred)
print(matrix)


import numpy as np
from keras.preprocessing import image
import tensorflow as tf
import os
import glob

print(test_set.class_indices)

folder_path = './Observed/Rhinella_granulosa'
os.chdir(folder_path)
files = glob.glob('*png')

for img in files:
	test_image = tf.keras.preprocessing.image.load_img(img, target_size = (422, 160),color_mode='grayscale')
	test_image = image.img_to_array(test_image)
	test_image = np.expand_dims(test_image, axis = 0)
	test_image = np.vstack([test_image])
	result = classifier.predict(test_image)
	result_ = classifier.predict_proba(test_image)
	if result[0][0] == 1:
		prediction = 'Model1'
	elif result[0][1] == 1:
		prediction = 'Model2'
	elif result[0][2] == 1:
		prediction = 'Model3'

	print(prediction)
	print(result_)


