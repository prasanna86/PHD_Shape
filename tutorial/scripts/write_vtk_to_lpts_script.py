import os
import math
from array import *
from vtk import *
import numpy

def vtk_to_array(vtk_array):
  #this is slow. numpy.zeros would be faster.
  r = array('f', [0]*vtk_array.GetSize())
  vtk_array.ExportToVoidPointer(r)
  return r

def main():

  ####################################################################
  ###        Ensure these paths are correct for your system        ###
  ####################################################################
  reg_path = "../reg_at_obs_time_pts/"
  output_base_path = "../shape_lpts_for_stats/"

  if not (os.path.isdir(output_base_path)):
    os.system('mkdir ' + output_base_path)

  subjects = ['50015', '50352', '50567', '50855', '50983', '51034', '51211', '51706', '51855', '51888', '51909', '52598', '52710', '52850']
  times = [['61.66461328', '63.18685832', '64.15605749'], ['61.87268994', '63.93976728'], ['58.09445585', '60.31485284'], ['59.56468172', '60.73100616', '61.5578371', '63.55099247', '64.44900753', '65.65639973'], ['55.2991102', '56.25188227', '57.25119781', '58.30527036', '59.41409993'], ['56.11225188', '57.27036277', '58.05886379'], ['56.10403833', '58.0971937', '59.05270363', '60.07118412'], ['55.84668036', '58.091718', '58.88569473', '59.79466119'], ['58.89664613', '59.89322382', '62.9431896'], ['58.37097878', '60.70910335'], ['54.97878166', '58.02600958', '59.00342231', '59.75085558'], ['57.94661191', '58.92128679', '59.88774812', '60.97193703'], ['59.55920602', '61.60985626'], ['58.90485969', '59.9945243', '60.93360712', '61.89459274', '62.7761807']]

  segs = ['left_caudate_time', 'right_caudate_time', 'left_putamen_time', 'right_putamen_time']

  for i in range(0, len(subjects)): 
    
    print "/======================================================="
    print "| Working on subject " + subjects[i]
    print "\======================================================="

    num_cur_timepts = len(times[i])
    min_time = float(times[i][0])
    time_range = float(times[i][num_cur_timepts-1]) - float(times[i][0])
    timepts_to_copy = []
  
    # We know we want the first and last time points, 
    # but also intermediate, which we must calculate
    for j in range(0, num_cur_timepts):

      time_pt = int(round(((float(times[i][j]) - min_time) / time_range)*29.0))
      timepts_to_copy.append(time_pt)

    cur_regression_dir = reg_path + subjects[i] + "/"
    cur_output_dir = output_base_path + subjects[i] + "/"

    # Create the output directory if needed
    if not (os.path.isdir(cur_output_dir)):
      os.makedirs(cur_output_dir)

    # Copy the data we need over from the regression directory
    for j in range(0, len(timepts_to_copy)):

      # There are several shapes
      for k in range(0, len(segs)):
        
        data_in = "%s%s_%0.2d.vtk" %(cur_regression_dir, segs[k], j)
        data_out = "%s%s_%0.2d.lpts" %(cur_output_dir, segs[k], j)
                                
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(data_in)
        reader.Update()
        polydata = reader.GetOutput()

        #Grab a scalar from the vtk file
        points_array = vtk_to_array(polydata.GetPoints().GetData())
        x = numpy.reshape(points_array, (642, 3))
        numpy.savetxt(data_out, x) 
    
if __name__ == "__main__":
  main()

