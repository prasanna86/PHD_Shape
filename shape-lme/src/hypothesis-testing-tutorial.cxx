#include "primitives.hxx"
#include <fstream>
#include <iostream>

typedef vector<MatrixType> MatListType;

// permute function
void permute(Col<unsigned int> & x)
{
  unsigned int i, j, sz;
  unsigned int swapTemp;
  double unifRand;

  sz = x.size();
  for(i = 0; i < sz - 1; i++)
  {
    unifRand = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    j = static_cast<unsigned int>(static_cast<double>(sz - i) * unifRand) + i;

    swapTemp = x[i];
    x[i] = x[j];
    x[j] = swapTemp;
  }
}

MatrixType procrustes(MatrixType & p, MatrixType & q)
{
  MatrixType m = p * strans(q);

  MatrixType U, V; VectorType s;
  svd_econ(U, s, V, m, 'b');

  MatrixType rotation = U * strans(V);
  MatrixType qRot = rotation * q;
  return qRot;
}

void runTests(MatrixType & X, MatrixType & design, TimePointsType & tp, 
	      MatrixType & fixed, MatrixType & random)
{
  //cout << "Estimating params using ME" << endl;
  VectorType expl; expl.zeros(design.n_rows);
  expl = design.col(1);

  unsigned int num_fixed = design.n_cols;

  // Number of samples
  unsigned int nsh = static_cast<double>(X.n_cols);
  unsigned int nr = X.n_rows; //number of points*3
  unsigned int nsub = tp.size();

  MatrixType identity_2;
  identity_2.eye(2,2);
  
  MatrixType * Ws = NULL, * Vs = NULL, * identity_n = NULL, * Xp = NULL, * Zp = NULL;
  Ws = new MatrixType[nsub];
  Vs = new MatrixType[nsub];
  identity_n = new MatrixType[nsub];
  Xp = new MatrixType[nsub];  // Xp: Design matrix. for fixed part
  Zp = new MatrixType[nsub];  // Zp: Design matrix. for random part
  
  VectorType * y = NULL, * residual = NULL;
  y = new VectorType[nsub];
  residual = new VectorType[nsub];

  int obs_thus_far = 0;

  for(int k = 0; k < nsub; k++)
  {
    Vs[k].set_size(tp[k], tp[k]);
    Ws[k].set_size(tp[k], tp[k]);
    
    identity_n[k].eye(tp[k], tp[k]);
    Xp[k].set_size(tp[k], num_fixed);
    Zp[k].set_size(tp[k], 2);

    y[k].zeros(tp[k]);
    residual[k].zeros(tp[k]);
  
    for (int l = 0; l < tp[k]; l++)
    {
      Zp[k](l, 0) = 1;
      Zp[k](l, 1) = expl(obs_thus_far + l);

      Xp[k](l, 0) = 1;
      for(int m = 1; m < num_fixed; m++)
      {
	Xp[k](l, m) = design(obs_thus_far + l, m);
      }
    }
 
    obs_thus_far += tp[k];
  }
  
  MatrixType sum_mat1; sum_mat1.zeros(num_fixed, num_fixed);
  VectorType sum_mat2; sum_mat2.zeros(num_fixed);
  double ecorr = 0.0;
  double tracevar = 0.0;
  MatrixType bscorr; bscorr.zeros(2, 2);
  MatrixType bsvar; bsvar.zeros(2, 2);
  VectorType tempvect_fixed; tempvect_fixed.zeros(num_fixed);
  VectorType tempvect_rand; tempvect_rand.zeros(num_fixed);

  double sigma2s = 0, sigma2s_new = 1; // variance of error
  MatrixType Ds, Ds_new; // covariance matrix of random parameters (2x2)
  Ds.zeros(2,2); Ds_new.eye(2,2);
  double epsilon1 = 1.0e-09, epsilon2 = 1.0e-05;
  double sigmaDiff = 2 * epsilon1, Ds_norm_diff = 2 * epsilon2;
  unsigned int em = 0;

  for(int i = 0; i < nr; i++) //for all points (x,y,z coordinates)
  {
    em = 0;
    obs_thus_far = 0;

    for(int k = 0; k < nsub; k++)
    {
      y[k].fill(0.0);
  
      for (int l = 0; l < tp[k]; l++)
      	y[k](l) = X(i, obs_thus_far + l);
 
      obs_thus_far += tp[k];
    }

    while(em < 500 && (sigmaDiff > epsilon1 || Ds_norm_diff > epsilon2)) // EM iterations
    {
      sigma2s = sigma2s_new;
      Ds = Ds_new;

      sum_mat1.fill(0.0); sum_mat2.fill(0.0);
      ecorr = 0.0;
      tracevar = 0.0;
      bscorr.fill(0.0);
      bsvar.fill(0.0);
     
      for (int k = 0; k < nsub; k++)
      {
	Vs[k] = (identity_n[k] * sigma2s) + Zp[k] * Ds * strans(Zp[k]);
        Ws[k] = inv(Vs[k]);
       
        sum_mat1 += strans(Xp[k]) * Ws[k] * Xp[k];
        sum_mat2 += strans(Xp[k]) * Ws[k] * y[k];
      }

      tempvect_fixed = solve(sum_mat1, sum_mat2);
      fixed.col(i) = tempvect_fixed;

      for (int k = 0; k < nsub; k++)
      {
        tempvect_rand = Ds * strans(Zp[k]) * Ws[k] * (y[k] - (Xp[k] * fixed.col(i)));
        random.col(i * nsub + k) = tempvect_rand;
	residual[k] = y[k] - (Xp[k] * fixed.col(i)) - (Zp[k] * random.col(i * nsub + k));	
        ecorr += dot(residual[k], residual[k]);
        tracevar += (tp[k] - sigma2s * trace(Ws[k]));
        bscorr += random.col(i * nsub + k) * strans(random.col(i * nsub + k));
        bsvar += (identity_2 - (strans(Zp[k]) * Ws[k] * Zp[k] * Ds));
      }
      
      sigma2s_new = (ecorr + sigma2s * tracevar) / static_cast<double>(nsh);
      Ds_new = (bscorr + Ds * bsvar) / static_cast<double>(nsub);

      sigmaDiff = abs(sigma2s_new - sigma2s);
      Ds_norm_diff = norm(Ds_new - Ds, "fro");

      em += 1;
    } //endfor EM iterations

    sigmaDiff = 2 * epsilon1;
    Ds_norm_diff = 2 * epsilon2;

    //cout << "particle, em iter, Ds_new norm, sigma2 = " << i+1 << " " << em << " " << norm(Ds_new, "fro") << " " << sigma2s_new << endl;
  } //endfor all points on shape (x,y & z)
  
  delete [] Vs;
  delete [] Ws;
  delete [] identity_n;
  delete [] Xp;
  delete [] Zp;
  delete [] y;
  delete [] residual;
}

int main(int argc, char ** argv)
{
  if(argc != 7)
  {
    cout << "Usage src/hypothesis-testing-tutorial param_file template_file" 
	 << "est_fixed est_random num_perm out_file"
	 << endl;
    return EXIT_FAILURE;
  }

  ofstream out_file; out_file.open(argv[6]);

  ifstream input_file, shape_file; string param_line, shape_line;
  unsigned int np = 0, nsh = 0;
  // read param file and count number of files
  input_file.open(argv[1]);
  while(getline(input_file, param_line))
    ++nsh;
  input_file.close();

  // read template file and count number of particles
  shape_file.open(argv[2]);
  while(getline(shape_file, shape_line))
    ++np;
  shape_file.close();

  out_file << "np, nsh = " << np << " " << nsh << endl;

  // variables to read param file
  string filename;
  unsigned int nsub = 0, num_fixed = 2, t = 0;
  ScalarType age, px, py, pz;
  ScalarType gid, sex, subj_num = 50015, prev_subj_num;

  // The Lists
  ScalarListType subj_nums; subj_nums.resize(nsh);
  ScalarListType groupIDs; groupIDs.resize(nsh);
  ScalarListType sex_of_subjs; sex_of_subjs.resize(nsh);
  ScalarListType ages; ages.resize(nsh);
  TimePointsType time_pts;
  MatrixType X; X.zeros(3*np, nsh);

  // open param file again
  input_file.open(argv[1]);
  istringstream ss, ss2;

  // initialize subj number
  prev_subj_num = subj_num; int i, j = 0;
  // reading param file
  out_file << "Reading param file: " << endl;
  while(getline(input_file, param_line))
  {
    ss.str(param_line);
    ss >> subj_num >> gid >> sex >> age >> filename;
    out_file << subj_num << " " << gid << " " << sex 
	     << " " << age << " " << filename << endl;
    
    subj_nums[j] = subj_num;
    groupIDs[j] = gid;
    sex_of_subjs[j] = sex;
    ages[j] = age;

    // opening points file
    shape_file.open(filename.c_str());
    out_file << "Reading shape file: " << filename << endl;
    i = 0;
    while(getline(shape_file, shape_line))
    {
      ss2.str(shape_line);
      ss2 >> px >> py >> pz;
      X(3 * i, j) = px;
      X(3 * i + 1, j) = py;
      X(3 * i + 2, j) = pz;
      ++i;
      ss2.clear();
    }
    shape_file.close();
    ++j; // incrementing shape column

    ++t; // incrementing time point
    if(subj_num != prev_subj_num)
    {
      time_pts.push_back(t); ++nsub;
      // re-intializing
      t = 0;
    }

    ss.clear();
    prev_subj_num = subj_num;
  }

  // if file ends, update timepoints with last subject
  time_pts.push_back(t); ++nsub; 
  
  out_file << "number of subjects = " << nsub << endl;

  // deciding number of fixed effects parameters based on categories
  ScalarListType sexes = unique(sex_of_subjs);
  ScalarListType groups = unique(groupIDs);

  // factor of 2 = intercept, slope. you can add quadratic, cubic as
  // well if that makes sense in which case the factor would be 3 or 4
  if(groups.size() != 1)
    num_fixed += 2 * (groups.size() - 1);
  if(sexes.size() != 1)
    num_fixed += 2 * (sexes.size() - 1);

  out_file << "number of fixed effects = " << num_fixed << endl;

  // creating design matrix
  MatrixType design; design.zeros(nsh, num_fixed);
  // first column: ones - coefficient: intercept
  design.col(0).ones();
  // second column: ages - coefficient: slope
  design.col(1) = ages;

  if(groups.size() == 2) // gid = 0, 1
  {
    design.col(2) = groupIDs;
    design.col(3) = groupIDs % ages;
  }
  else if(groups.size() > 2) // gid = 0, 1, 2, etc
  {
    for(i = 1; i < groups.size(); i++)
    {
      design.col(2*i) = groupIDs;
      VectorType v = design.col(2*i);
      v.elem(find(v != i)).zeros();
      v.elem(find(v == i)).ones();
      design.col(2*i) = v;
      design.col(2*i+1) = v % ages;
    }
  }

  if(sexes.size() == 2) // sex = female (0), male (1)
  {
    design.col(num_fixed-2) = sex_of_subjs;
    design.col(num_fixed-1) = sex_of_subjs % ages;
  }

  out_file << "design = " << design << endl;
 
  // align shapes to template
  MatrixType template_shape; template_shape.zeros(3, np);
  shape_file.open(argv[2]);
  i = 0;
  while(shape_file >> px >> py >> pz)
  {
    template_shape(0, i) = px;
    template_shape(1, i) = py;
    template_shape(2, i) = pz;
    ++i;
  }
  shape_file.close();

  // out_file << "template = " << endl << strans(template_shape) << endl;

  // removing translations from template
  VectorType mean_xyz; mean_xyz.zeros(3); 
  for(i = 0; i < np; i++)
    mean_xyz += template_shape.col(i);

  mean_xyz /= static_cast<double>(np);

  for(i = 0; i < np; i++)
    template_shape.col(i) = template_shape.col(i) - mean_xyz;

  MatrixType Y; Y.zeros(3*np, nsh);
  // ...and align everything to the template
  for(j = 0; j < nsh; j++)
  {
    MatrixType current_shape = X.col(j);
    current_shape.reshape(3, np);

    // removing translations from current shape
    mean_xyz.zeros(3); 
    for(i = 0; i < np; i++)
      mean_xyz += current_shape.col(i);

    mean_xyz /= static_cast<double>(np);

    for(i = 0; i < np; i++)
      current_shape.col(i) = current_shape.col(i) - mean_xyz;

    // procrustes align to template
    current_shape = procrustes(template_shape, current_shape);
    
    // reshape back
    current_shape.reshape(3 * np, 1);
    Y.col(j) = current_shape;
  }

  // parameters of interest
  MatrixType fixed; // slopes + intercepts for all points
  MatrixType random; // slopes + intercepts for all groups, for all points
  fixed.set_size(num_fixed, 3 * np);
  random.set_size(2, 3 * np * nsub);
  
  // estimate model with aligned shapes
  runTests(Y, design, time_pts, fixed, random);

  // output fixed params
  ofstream out;
  out.open(argv[3]);
 
  for(i = 0; i < np * 3; i++)
  {
    for(j = 0; j < num_fixed; j++)
      out << fixed(j, i) << " "; 
    out << endl;
  }

  out.close();
  // output random params
  out.open(argv[4]);

  for(i = 0; i < nsub * np * 3; i++)
  {
    for(j = 0; j < 2; j++)
      out << random(j, i) << " ";
    out << endl;
  }

  out.close();

  // ------------ Hypothesis testing: control vs hd --------------
  int numPermutes = atoi(argv[5]);

  double origTestStat_int = norm(fixed.row(2), "fro");
  double origTestStat_slope = norm(fixed.row(3), "fro");
  out_file << "orig test stat = " << origTestStat_int 
	   << " " << origTestStat_slope << endl;

  MatrixType subj_group_info; subj_group_info.zeros(nsub, 2);
  ScalarType obs_thus_far = 0;
  for(i = 0; i < nsub; i++)
  {
    subj_group_info(i, 0) = groupIDs[obs_thus_far];
    subj_group_info(i, 1) = obs_thus_far;
    obs_thus_far += time_pts[i];
  }

  uvec ctrl = find(subj_group_info.col(0) == 0);
  uvec hd = find(subj_group_info.col(0) == 1);

  // out_file << "subj_group_info = " << endl << subj_group_info << endl;
  // out_file << "ctrl = " << strans(ctrl) << endl;
  // out_file << "hd = " << strans(hd) << endl;

  Col<unsigned int> x; x.zeros(ctrl.size() + hd.size());
  for(i = 0; i < x.size(); i++)
  {
    if(i < ctrl.size())
      x[i] = ctrl[i];
    else
      x[i] = hd[i - ctrl.size()];
  }

  unsigned int numBelow_int = 0, numBelow_slope = 0;
  MatrixType new_fixed, new_random;
  new_fixed.set_size(num_fixed, 3 * np);
  new_random.set_size(2, 3 * np * nsub);

  srand(time(NULL));

  for(i = 0; i < numPermutes; i++)
  {
    permute(x);
    MatrixType new_design; new_design.zeros(nsh, num_fixed);
    new_design = design;

    for(j = 0; j < ctrl.size(); j++)
    {
      for(int k = 0; k < time_pts[x[j]]; k++)
      {
	new_design(subj_group_info(x[j], 1) + k, 2) = 0;
	new_design(subj_group_info(x[j], 1) + k, 3) = 0;
      }
    }

    for(j = ctrl.size(); j < x.size(); j++)
    {
      for(int k = 0; k < time_pts[x[j]]; k++)
      {
	new_design(subj_group_info(x[j], 1) + k, 2) = 1;
	new_design(subj_group_info(x[j], 1) + k, 3) = new_design(subj_group_info(x[j], 1) + k, 1);
      }
    }

    runTests(Y, new_design, time_pts, new_fixed, new_random);
    double testStat_int = norm(new_fixed.row(2), "fro");
    double testStat_slope = norm(new_fixed.row(3), "fro");

    // intercept
    if(testStat_int < origTestStat_int)
      numBelow_int++;
    out_file << i << ": Test stat, orig (int) = " << testStat_int << "\t" << origTestStat_int << "\t"
	    << 1 - ((double) numBelow_int / (double)(i + 1)) << endl;

    // slope
    if(testStat_slope < origTestStat_slope)
      numBelow_slope++;
    out_file << i << ": Test stat, orig (slope) = " << testStat_slope << "\t" << origTestStat_slope << "\t"
	    << 1 - ((double) numBelow_slope / (double)(i + 1)) << endl;
  }

  out_file << "p value (int) = " << 1 - ((double)(numBelow_int) / numPermutes) << endl;
  out_file << "p value (slope) = " << 1 - ((double)(numBelow_slope) / numPermutes) << endl;
  out_file.close();

  return EXIT_SUCCESS;
}
