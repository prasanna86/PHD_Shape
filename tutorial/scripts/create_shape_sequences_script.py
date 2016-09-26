import os
import math
import csv
import sys
import vtk

def main():

  ####################################################################
  ###        Ensure these paths are correct for your system        ###
  ####################################################################
  input_path = '../mixed_effects_results/'
  output_path = '../output_shape_seq_vtk/'
  vtk_template_path = '../reg_at_obs_time_pts/'

  if not (os.path.isdir(output_path)):
    os.system('mkdir ' + output_path)

  struct = ['left_caudate', 'right_caudate', 'left_putamen', 'right_putamen']
  num_shape_seq = 41
  subjects = ['50015', '50352', '50567', '50855', '50983', '51034', '51211', '51706', '51855', '51888', '51909', '52598', '52710', '52850']
  times = [['61.66461328', '63.18685832', '64.15605749'], ['61.87268994', '63.93976728'], ['58.09445585', '60.31485284'], ['59.56468172', '60.73100616', '61.5578371', '63.55099247', '64.44900753', '65.65639973'], ['55.2991102', '56.25188227', '57.25119781', '58.30527036', '59.41409993'], ['56.11225188', '57.27036277', '58.05886379'], ['56.10403833', '58.0971937', '59.05270363', '60.07118412'], ['55.84668036', '58.091718', '58.88569473', '59.79466119'], ['58.89664613', '59.89322382', '62.9431896'], ['58.37097878', '60.70910335'], ['54.97878166', '58.02600958', '59.00342231', '59.75085558'], ['57.94661191', '58.92128679', '59.88774812', '60.97193703'], ['59.55920602', '61.60985626'], ['58.90485969', '59.9945243', '60.93360712', '61.89459274', '62.7761807']]
        
  ## Time range is min to max of time data: you can change this
  min_time = float(min(min(times)))
  max_time = float(max(max(times)))

  # lets create shape sequences by structure
  for k in range(0, len(struct)):

    print "/======================================================="
    print "| Working on the " + struct[k]
    print "\======================================================="

    cur_input_dir = input_path + struct[k] + '/'
    cur_output_dir = output_path + struct[k] + '/'

    # Create the output directory if needed
    if not (os.path.isdir(cur_output_dir)):
            os.makedirs(cur_output_dir)

    template_file = '%s%s_time_00.vtk' %(vtk_template_path, struct[k])
    head = []
    bottom = []
    with open(template_file) as input_data:
      # Skips text before the beginning of the interesting block:
      for line in input_data:
        head.append(line)  # block_of_lines.append(line)
        if line.strip() == 'POINTS 642 float':
          break

      for line in input_data:
        if line.strip() == 'POLYGONS 1280 5120':
          bottom.append(line)
          break

      for line in input_data:
        bottom.append(line)  # block_of_lines.append(line)
        if line.strip() == 'VECTORS velocity double':
          break
  
    # read in the fixed-effects file
    fe_file = '%sfixed_effects_%s.txt' %(cur_input_dir, struct[k])
    with open(fe_file, 'r') as f:
      content = f.readlines()

    for j in range(0, num_shape_seq):

      ctrl_seq_file = '%sctrl_%s.vtk' %(cur_output_dir, str(j))
      hd_seq_file = '%shd_%s.vtk' %(cur_output_dir, str(j))
      time_pt = min_time + float(j * (max_time - min_time)) / float(num_shape_seq)
      
      # Now, print the data back to a file how you'd like
      ctrlout = open(ctrl_seq_file,'w')
      hdout = open(hd_seq_file,'w')

      for line in head:
        ctrlout.write(line)
        hdout.write(line)
      
      count = 0
      
      for line in content:
              
        if(count == 0):
          row1 = line.split()
        elif(count == 1):
          row2 = line.split()
        elif(count == 2):
          row3 = line.split()

        count = count + 1

        if(count == 3):
          ctrl1 = float(row1[0]) + time_pt * float(row1[1])
          ctrl2 = float(row2[0]) + time_pt * float(row2[1])
          ctrl3 = float(row3[0]) + time_pt * float(row3[1])
          hd1 = float(row1[0]) + float(row1[2]) + time_pt * (float(row1[1]) + float(row1[3]))
          hd2 = float(row2[0]) + float(row2[2]) + time_pt * (float(row2[1]) + float(row2[3]))
          hd3 = float(row3[0]) + float(row3[2]) + time_pt * (float(row3[1]) + float(row3[3]))
          # For specific columns
          ctrlout.write('{0} {1} {2}\n'.format(ctrl1, ctrl2, ctrl3))
          hdout.write('{0} {1} {2}\n'.format(hd1, hd2, hd3))
          count = 0

      ctrlout.write('\n')
      hdout.write('\n')

      for line in bottom:
        ctrlout.write(line)
        hdout.write(line)

      for line in content:
              
        if(count == 0):
          row1 = line.split()
        elif(count == 1):
          row2 = line.split()
        elif(count == 2):
          row3 = line.split()

        count = count + 1

        if(count == 3):
          ctrl1 = float(row1[1])
          ctrl2 = float(row2[1])
          ctrl3 = float(row3[1])
          hd1 = float(row1[1]) + float(row1[3])
          hd2 = float(row2[1]) + float(row2[3])
          hd3 = float(row3[1]) + float(row3[3])
          # For specific columns
          ctrlout.write('{0} {1} {2}\n'.format(ctrl1, ctrl2, ctrl3))
          hdout.write('{0} {1} {2}\n'.format(hd1, hd2, hd3))
          count = 0

      ctrlout.close()
      hdout.close()

      # Calculate Point Normals for each of these vtk files
      polyReader = vtk.vtkPolyDataReader()
      polyReader.SetFileName( ctrl_seq_file )
      polyReader.Update()
      ctrlData = polyReader.GetOutput()

      polyReader2 = vtk.vtkPolyDataReader()
      polyReader2.SetFileName( hd_seq_file )
      polyReader2.Update()
      hdData = polyReader2.GetOutput()
      
      polyNormals_ctrl = vtk.vtkPolyDataNormals()
      polyNormals_ctrl.SetInputData( ctrlData )
      polyNormals_ctrl.ComputeCellNormalsOff()
      polyNormals_ctrl.ComputePointNormalsOn()
      polyNormals_ctrl.ConsistencyOn()
      # Splitting Off to make the number of normals consistent with 'velocity' vectors
      polyNormals_ctrl.SplittingOff()
      polyNormals_ctrl.Update()
      ctrlData_wPtNormals = polyNormals_ctrl.GetOutput()

      polyNormals_hd = vtk.vtkPolyDataNormals()
      polyNormals_hd.SetInputData( hdData )
      polyNormals_hd.ComputeCellNormalsOff()
      polyNormals_hd.ComputePointNormalsOn()
      polyNormals_hd.ConsistencyOn()
      # Splitting Off to make the number of normals consistent with 'velocity' vectors
      polyNormals_hd.SplittingOff()
      polyNormals_hd.Update()
      hdData_wPtNormals = polyNormals_hd.GetOutput()

      # Calculate expansion by computing dot product between velocity and normals
      ExpansionFactor_ctrl = vtk.vtkFloatArray()
      ExpansionFactor_ctrl.SetName( "Expansion" )

      ExpansionFactor_hd = vtk.vtkFloatArray()
      ExpansionFactor_hd.SetName( "Expansion" )

      # Read Point Normals 
      Normals_ctrl = ctrlData_wPtNormals.GetPointData().GetArray( "Normals" )
      Normals_hd = hdData_wPtNormals.GetPointData().GetArray( "Normals" )

      # Read Point Normals 
      velocity_ctrl = ctrlData_wPtNormals.GetPointData().GetArray( "velocity" )
      velocity_hd = hdData_wPtNormals.GetPointData().GetArray( "velocity" )

      for i in xrange(ctrlData_wPtNormals.GetNumberOfPoints() ):
        Nc = Normals_ctrl.GetTuple( i )
        Nh = Normals_hd.GetTuple( i )
        vc = velocity_ctrl.GetTuple( i )
        vh = velocity_hd.GetTuple( i )

        e_ctrl = ( vc[ 0 ] * Nc[ 0 ] ) + ( vc[ 1 ] * Nc[ 1 ] ) + ( vc[ 2 ] * Nc[ 2 ] ) 
        ExpansionFactor_ctrl.InsertNextValue( e_ctrl )

        e_hd = ( vh[ 0 ] * Nh[ 0 ] ) + ( vh[ 1 ] * Nh[ 1 ] ) + ( vh[ 2 ] * Nh[ 2 ] ) 
        ExpansionFactor_hd.InsertNextValue( e_hd )

      ctrlData_wPtNormals.GetPointData().AddArray( ExpansionFactor_ctrl )
      #ctrlData_wPtNormals.Update()
      hdData_wPtNormals.GetPointData().AddArray( ExpansionFactor_hd )
      #hdData_wPtNormals.Update()

      polyWriter_ctrl = vtk.vtkPolyDataWriter()
      polyWriter_ctrl.SetFileName( ctrl_seq_file )
      polyWriter_ctrl.SetInputData( ctrlData_wPtNormals )
      polyWriter_ctrl.Write()
      polyWriter_ctrl.Update()

      polyWriter_hd = vtk.vtkPolyDataWriter()
      polyWriter_hd.SetFileName( hd_seq_file )
      polyWriter_hd.SetInputData( hdData_wPtNormals )
      polyWriter_hd.Write()
      polyWriter_hd.Update()

    time_seq_file = '%stime_seq.txt' %(cur_output_dir)
    fout = open(time_seq_file, 'w')
    for j in range(0, num_shape_seq):
      time_pt = min_time + float(j * (max_time - min_time)) / float(num_shape_seq)
      fout.write('{0}\n'.format(time_pt))
    fout.close()

if __name__ == "__main__":
  main()
