import os
import math

def main():

  ####################################################################
  ###        Ensure these paths are correct for your system        ###
  ####################################################################
  input_path = '../param_files/'
  output_path = '../mixed_effects_results/'
  hypothesis_test_path = '../../bin/shape-lme/hypothesis-testing-tutorial'
  template_file_dir = '../shape_lpts_for_stats/50855/'

  subjects = ['50015', '50352', '50567', '50855', '50983', '51034', '51211', '51706', '51855', '51888', '51909', '52598', '52710', '52850']
  num_perm = 100

  segs = ['left_caudate_time_', 'right_caudate_time_', 'left_putamen_time_', 'right_putamen_time_']
  struct = ['left_caudate', 'right_caudate', 'left_putamen', 'right_putamen']

  # lets do stats structure by structure
  for k in range(0, len(struct)):

    print "/======================================================="
    print "| Working on the " + struct[k]
    print "\======================================================="
    
    param_file = '%s%s_param_file.txt' %(input_path, segs[k])
    cur_io_dir = output_path + struct[k] + '/'
    # Create the output directory if needed
    if not (os.path.isdir(cur_io_dir)):
      os.makedirs(cur_io_dir)

    # Create the lpts_file
    fixed = '%sfixed_effects_%s.txt' %(cur_io_dir, struct[k])
    random = '%srandom_effects_%s.txt' %(cur_io_dir, struct[k])
    template = '%s%s_time_00.lpts' %(template_file_dir, struct[k])
    out_file = '%shypothesis_ctrl_high.txt' %(cur_io_dir)

    # Create test path
    test_path = hypothesis_test_path + ' ' + param_file + ' ' + template + ' ' + fixed + ' ' + random + ' ' + str(num_perm) + ' ' + out_file + ' &'
    os.system(test_path)
                
if __name__ == "__main__":
  main()

