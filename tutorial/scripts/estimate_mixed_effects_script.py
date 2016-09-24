import os
import math

def main():

  ####################################################################
  ###        Ensure these paths are correct for your system        ###
  ####################################################################
  input_path = '../param_files/'
  output_path = '../mixed_effects_results/'
  me_test_path = '../../bin/shape-lme/estimate-covariate-mixed-effects'

  if not (os.path.isdir(output_path)):
    os.system('mkdir ' + output_path)

  subjects = ['50015', '50352', '50567', '50855', '50983', '51034', '51211', '51706', '51855', '51888', '51909', '52598', '52710', '52850']
    
  segs = ['left_caudate_time_', 'right_caudate_time_', 'left_putamen_time_', 'right_putamen_time_']
  struct = ['left_caudate', 'right_caudate', 'left_putamen', 'right_putamen']

  # lets do stats structure by structure
  for k in range(0, len(struct)):

    print "/======================================================="
    print "| Working on the " + struct[k]
    print "\======================================================="
    
    param_file = '%s%s_param_file.txt' %(input_path, segs[k])
    template_file = '%s%s/%s_00.lpts' %(input_path, subjects[3], segs[k])

    cur_output_dir = output_path + struct[k] + '/'
    # Create the output directory if needed
    if not (os.path.isdir(cur_output_dir)):
      os.makedirs(cur_output_dir)

    # Create the lpts_file
    out_fixed_effects = '%sfixed_effects_%s.txt' %(cur_output_dir, struct[k])
    out_random_effects = '%srandom_effects_%s.txt' %(cur_output_dir, struct[k])

    os.system(me_test_path + ' ' + param_file + ' ' + template_file + ' ' + out_fixed_effects + ' ' + out_random_effects)

if __name__ == "__main__":
  main()

