import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import skimage.io

from models.Graph_Model_v2 import GraphCut
from models.Graph_Model import Kruskal


def file_write(data, file_name):
    if len(data.shape) == 3:
        for channel in range(data.shape[2]):
            f = open(file=file_name + " channel" + str(channel) + ".data", mode="w+")
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    f.write(str(data[i, j, channel]) + " ")
                f.write("\n")
            f.close()
    elif len(data.shape) == 2:
        f = open(file=file_name + ".data", mode="w+")
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                f.write(str(data[i, j]) + " ")
            f.write("\n")
        f.close()
    else:
        print("Not enough dimensions\n")

def rgb2gray(rgb):
    return np.dot(rgb[..., :3], [0.299, 0.587, 0.114])


def Main():
    image = skimage.io.imread(fname="./resource/depthMap_Color.png")#image.jpg
    image = image / 16.0
    image = (image[:, :].astype(np.int))
    # noise = np.random.normal(loc=0.0, scale=64.0, size=(image.shape[0], image.shape[1]))
    # q = 0.9
    # noise = np.random.random(size = image.size).reshape(image.shape) > q
    noise = np.random.randn(image.shape[0], image.shape[1])
    sigma = 1.0
    noise = sigma * noise
    noise = noise.astype(int)
    #noise = image + noise
	#noise = image

    #obj = Kruskal(noise, 20, 1)
    processed = image#obj.get_image()
    fig, axes = plt.subplots(figsize=(10, 6))
    axes.imshow(processed, cmap=cm.gray, aspect="equal", interpolation="none", vmin=-16.0, vmax=16.0)
    plt.show()

    noise_int = np.random.randint(low=0, high=32, size=(noise.shape[0], noise.shape[1]))
    graph = GraphCut(10, image, processed)
    graph.solve()

    print("Solved")

if __name__ == "__main__":
    Main()
