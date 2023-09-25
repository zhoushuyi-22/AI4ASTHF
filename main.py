from keras.models import load_model
import scipy.io as sio
import os
from loss_function import *
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

index = 0  # SHF is 0, LHF is 1
if index == 1:
    data = sio.loadmat('./dataset_LHF.mat')
    location = 'LHF_PINN.hdf5'
else:
    data = sio.loadmat('./dataset_SHF.mat')
    location = 'SHF_PINN.hdf5'

test_x = data['test_x']
test_y = data['test_y']
model = load_model(location, custom_objects={'myself_loss': myself_loss})
pred = model.predict(test_x)
mse, mae, r2, bias, std = loss_functions(pred, test_y)


