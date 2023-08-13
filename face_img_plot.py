# -*- coding: utf-8 -*-
"""
re-construct face images from vectors 
"""
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

# import the RFBF obtained from R
data = pd.read_csv('FBF_sphere_test_h10_r20.csv')

# image re-construction begins here
img = np.asarray(data)

print(img.shape)

d,n=img.shape

for i in range(n):
  img_1=img[:,i]
  img_1=img_1.reshape(50,37)
  plt.imshow(img_1,cmap='gray')
  file_name = str(i)+'.png'
  plt.savefig(file_name, bbox_inches='tight')
  plt.axis('off')
plt.show()
