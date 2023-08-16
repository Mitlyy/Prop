import matplotlib.pyplot as plt
import numpy as np

from PIL import Image
png_filepath = 'photo_2023-05-31_10-01-50.jpg'
png_pil_img = Image.open(png_filepath).convert('L')
# this will print info about the PIL object
print(png_pil_img.format, png_pil_img.size, png_pil_img.mode)
png_np_img = np.asarray(png_pil_img)
png_line = png_np_img.reshape(-1)
fig, ax = plt.subplots(4, 1, figsize=(8,8))
ax[0].imshow(png_np_img, cmap='gray')
ax[1].hist(png_line, bins = 100)
# ccr = np.dot(png_np_img,png_np_img.T)
# ax[2,0].imshow(ccr, cmap='gray')
# ax[3,0].hist(ccr.reshape(-1), bins = 100)

# png_filepath = 'IMG_4800.JPG'
# png_pil_img = Image.open(png_filepath).convert('L')
# # this will print info about the PIL object
# print(png_pil_img.format, png_pil_img.size, png_pil_img.mode)
# png_np_img = np.asarray(png_pil_img)
# png_line = png_np_img.reshape(-1)
# ax[0,1].imshow(png_np_img, cmap='gray')
# ax[1,1].hist(png_line, bins = 100)
# ccr = np.dot(png_np_img,png_np_img.T)
# ax[2,1].imshow(ccr, cmap='gray')
# ax[3,1].hist(ccr.reshape(-1), bins = 100)


def corr(img_1, img_2):
	blank_img = np.likezeros(img_1)
	iter_d = 0
	for x in range(img_1.shape[0]):
		for y in range(img_1.shape[1]):
			blank_img[x,y] += iter_d
			iter_d = 0
			for k in range(img_1.shape[0]):
				for l in range(img_1.shape[1]):
					iter_d += img_1[x+k, y+l]*img_2[k, l]



#
# corr = ...
# fig, ax = plt.subplots(2, 1, figsize=(8,8))
# ax[0].imshow(corr, cmap='gray')
# ax[1].hist(corr, bins = 100)
plt.show()