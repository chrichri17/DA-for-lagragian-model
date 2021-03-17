import os
import cv2
import matplotlib.pyplot as plt

__file__ = "contour_detection.py"
dir_path = os.path.dirname(os.path.realpath(__file__))


image = cv2.imread(dir_path+"/fire.jpg")
# convert to RGB
image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
# convert to grayscale
gray = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)

# create a binary thresholded image
_, binary = cv2.threshold(gray, 225, 255, cv2.THRESH_BINARY_INV)
# show it
plt.imshow(binary, cmap="gray")
plt.show()

# show the image with the drawn contours
plt.imshow(image)
plt.show()