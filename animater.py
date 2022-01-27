import os
import sys
from the_fluid import The_Fluid
import os; cur_path = os.getcwd()
if sys.version_info.major == 3:
  import _pickle as cPickle
else:
  import cPickle
# end if/else
my_fluid = cPickle.load(open("my_fluid.p", "rb"))

FRAMERATE = 60

input_name  = my_fluid.PLOT_DIR + "/hat_periodic_boundaries_cycle_%09d.png"
output_name = my_fluid.cmap.name + ".mp4"

os.system("ffmpeg -y -framerate " + str(FRAMERATE) + " -i " + input_name 
          + " -pix_fmt yuv420p " + output_name)
## end animater.py
