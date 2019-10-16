import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt
import keras
from keras.layers import Dense
from keras.models import Sequential
from keras.models import load_model
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

# collecting trajectory data
df = pd.read_csv("/Users/admin/Desktop/V-formation/AMPC/newtraj.csv")

# verifying data
print(df.head())

# positions and velocities of all the agents
inputs = ["x1","x2","x3","x4","x5","x6","x7","y1","y2","y3","y4","y5","y6","y7",
          "vx1","vx2","vx3","vx4","vx5","vx6","vx7","vy1","vy2","vy3","vy4","vy5","vy6","vy7","J"]

# acceleration of an agent
outputs = ["ax1","ay1"]

# convert into dataframes
x = df[inputs]
y = df[outputs]

# adjust for NN usage
X = x.values
Y = y.values

# fit data for training
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1)
sc = StandardScaler()
X_train = sc.fit_transform(X_train)

# DNN training
model = Sequential()

# gives the number of input features
n_cols = X.shape[1]

# hidden layers
model.add(Dense(128, activation='sigmoid', input_shape=(n_cols,)))
model.add(Dense(128))
model.add(Dense(128))
model.add(Dense(128))
model.add(Dense(128))

# output layer
model.add(Dense(2))

# optimization
adam = keras.optimizers.Adam(lr=0.0001, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.0, amsgrad=False)
model.compile(optimizer='adam', loss='mean_squared_error', metrics =['accuracy'])
model.summary()

start = time.perf_counter()

# training
history = model.fit(X_train, y_train, validation_data=(X_test,y_test), batch_size=1000 , epochs=100000, verbose = 2)

end = time.perf_counter()

# output training time
print(end-start)

# save the NN model
model.save("nerual_controller.h5")

# training error visualization
plt.plot(history.history['loss'])
plt.title('Model Loss')
plt.ylabel('Error')
plt.xlabel('Epoch')
plt.show()
