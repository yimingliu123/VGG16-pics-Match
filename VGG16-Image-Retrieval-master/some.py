# -*- coding=utf-8 -*-
import math
from scipy.misc import imread, imresize
# Get number of images to match (default 4)
#dist_type = input("Enter distance algorithm (euc, cos, chev): \n") or "euc"
#print("distance type selected: " + dist_type)
cor_list = []
for x in range(5, 5 + 4):
    cor_list.append("image_" + str(x).zfill(5) + ".jpg")
print('wo shi baba:'+str(cor_list))

img_query1 = imread('image_00001.jpg')
img_query1 = imresize(img_query1, (224, 224))
print(img_query1.shape)

img_query2 = imread('image_00076.jpg')
img_query2 = imresize(img_query2, (224, 224))
print(img_query2.shape)