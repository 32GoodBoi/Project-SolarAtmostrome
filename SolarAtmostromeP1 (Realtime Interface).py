#------------------------------Libraries---------------------------------
import cv2
from skimage import io
from matplotlib import pyplot as plt
import numpy as np

#------------------------------Define part-------------------------------

A193 = io.imread('https://sdo.gsfc.nasa.gov/assets/img/latest/latest_1024_0193.jpg')
OGSun = io.imread('https://sdo.gsfc.nasa.gov/assets/img/latest/latest_1024_0304.jpg')
Sursun = io.imread('https://sdo.gsfc.nasa.gov/assets/img/latest/f_HMImag_171_1024.jpg')
HMIB = io.imread('https://sdo.gsfc.nasa.gov/assets/img/latest/latest_1024_HMIBC.jpg')
graph1 = io.imread('https://sdo.gsfc.nasa.gov/assets/img/latest/latest_ar_map.png')
graph2 = io.imread('https://sdo.gsfc.nasa.gov/assets/img/latest/latest_composite_map.png')
graph3 = io.imread('https://stereo-ssc.nascom.nasa.gov/beacon/euvi_195_heliographic.gif')
(h,w) = A193.shape[:2]
cenX, cenY = (w//2), (h//2)
tl = A193[0:cenY, 0:cenX]
tr = A193[0:cenY, cenX:w]
bl = A193[cenY:h, 0:cenX]
br = A193[cenY:h, cenX:w]


#--------------------------conditions------------------------------------
conclusion = 'Solar wind expected soon (in 3 to 6 hrs)'


#------------------------------------------------------------------------
plt.subplot(2,2,2)
nOGSun = cv2.cvtColor(OGSun, cv2.COLOR_BGR2GRAY)
OGhist = cv2.calcHist([nOGSun],[0],None,[256],[0,256])

fig = plt.figure()
ax = fig.add_subplot()
fig.set_figheight(20)
fig.set_figwidth(20)
fig. suptitle("Solar Atmostrome; Dataset(live)", fontsize=15)
 
ax1 = plt.subplot2grid(shape=(4, 4), loc=(0, 0), colspan=1) #ogsun
ax2 = plt.subplot2grid(shape=(4, 4), loc=(1, 0), colspan=1) #a193
ax3 = plt.subplot2grid(shape=(4, 4), loc=(1, 1), rowspan=1) #t1
ax4 = plt.subplot2grid((4, 4), (2, 0)) #hmib
ax5 = plt.subplot2grid((4, 4), (2, 1), colspan=1)#sursun
ax6 = plt.subplot2grid((3, 3), (0, 2), colspan=2)#graph1
ax7 = plt.subplot2grid((3, 3), (1, 2), colspan=2)#graph2
ax8 = plt.subplot2grid((3, 3), (2, 2), colspan=2)#graph3
ax9 = plt.subplot2grid((4, 4), (3, 0), colspan=2)#conclusion

# initializing x,y axis value
x = np.arange(0, 10, 0.1)
y = np.cos(x)
 
# plotting subplots
ax1.imshow(OGSun)
ax1.set_title('What you normally see')
ax2.imshow(A193)
ax2.set_title('This is how it is beyond the vision')
ax3.imshow(tl)
ax3.set_title('Incoming flare status;')
ax4.imshow(HMIB)
ax4.set_title("Magnetic Index")
ax5.imshow(Sursun)
ax5.set_title('solar radiation at 171 Å')
ax6.imshow(graph1)
ax6.set_title('Ar Mapping of data')
ax7.imshow(graph2)
ax7.set_title('Composite Mapping')
ax8.imshow(graph3)
ax8.set_title('Heliographic map')


ax9.text(1/10, 1/2, conclusion, fontsize=15,  fontweight="bold")
ax9.set_title('Prediction')

# automatically adjust padding horizontally
# as well as vertically.
plt.tight_layout()
 
# display plot
plt.show()







