import os
import math

def main():

  #######################################################################
  ### Ensure these paths are correct for your system / Edit if needed ###
  #######################################################################
  paramfile_path = '../param_files/'
  template_path = '../shape_lpts_for_stats/50855/'
  output_path = '../mixed_effects_results/'
  me_test_path = '../../bin/shape-lme/estimate-covariate-mixed-effects'

  if not (os.path.isdir(output_path)):
    os.system('mkdir ' + output_path)

  segs = ['left_caudate_time_', 'right_caudate_time_', 'left_putamen_time_', 'right_putamen_time_']
  struct = ['left_caudate', 'right_caudate', 'left_putamen', 'right_putamen']

  # lets do stats structure by structure
  for k in range(0, len(struct)):

    print "/======================================================="
    print "| Working on the " + struct[k]
    print "\======================================================="
    
    param_file = '%s%sparam_file.txt' %(paramfile_path, segs[k])
    template_file = '%s%s00.lpts' %(template_path, segs[k])

    # Create the output directory if needed
    cur_output_dir = output_path + struct[k] + '/'
    if not (os.path.isdir(cur_output_dir)):
      os.makedirs(cur_output_dir)

    # Create output files
    out_fixed_effects = '%sfixed_effects_%s.txt' %(cur_output_dir, struct[k])
    out_random_effects = '%srandom_effects_%s.txt' %(cur_output_dir, struct[k])

    os.system(me_test_path + ' ' + param_file + ' ' + template_file + ' ' + out_fixed_effects + ' ' + out_random_effects)

if __name__ == "__main__":
  main()

