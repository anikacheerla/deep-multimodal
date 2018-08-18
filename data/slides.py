
import numpy as np
import matplotlib.pyplot as plt

import os, sys, random, yaml, shutil, time
import IPython

import openslide
from scipy.misc import imresize

i = 0
def open_slide(slide_image_path, allow_copy=True):	
	if "cloud-data" in slide_image_path:
		if not allow_copy: return None
		shutil.copy("data/" + slide_image_path, "cache.svs")
		slide_image_path = "cache" + str(i) + ".svs"
		print ("Had to copy")
		i += 1

	return openslide.open_slide(slide_image_path)
	

def sample_from_slides(slides, window_size=500, view_size=1000, num=10, max_num=40, tiling=2):

	samples = []
	i = 0

	while len(samples) < num:
		slide = random.choice(slides)
		XMAX, YMAX = slide.dimensions[0], slide.dimensions[1]

		tile_size = window_size*tiling
		xv, yv = random.randint(0, XMAX - tile_size - 5), random.randint(0, YMAX - tile_size - 5)
		window = slide.read_region((xv, yv), 0, (tile_size, tile_size))
		window = np.array(window)

		i += 1
		if i == max_num:
			return None

		if window.mean() > 175 or window.mean() < 50:
			continue

		if window[:, :, 2].mean() > window[:, :, 0].mean() + 30.0:
			continue

		window = imresize(np.array(window), (view_size*tiling, view_size*tiling, 4))
		
		for sx in range(0, view_size*tiling, view_size):
			for sy in range(0, view_size*tiling, view_size):
				tile = window[sx:(sx+view_size), sy:(sy+view_size), 0:3]
				if tile.mean() > 165 or tile.mean() < 80:
					continue

				#print (tile[:, :, 2].mean(), tile[:, :, 0].mean())
				#plt.imshow(tile); plt.show()
				samples.append(tile)

	samples = samples[0:num]
	return np.array(samples)




if __name__ == "__main__":
	IPython.embed()
	for case in sorted(data.keys()):
		print (case, i)
		sample_from_patient(case)

