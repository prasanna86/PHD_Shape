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
    Vs[k].set_size(tp(k), tp(k));
    Ws[k].set_size(tp(k), tp(k));
    
    identity_n[k].eye(tp(k), tp(k));
    Xp[k].set_size(tp(k), num_fixed);
    Zp[k].set_size(tp(k), 2);

    y[k].zeros(tp(k));
    residual[k].zeros(tp(k));
  
    for (int l = 0; l < tp(k); l++)
    {
      Zp[k](l, 0) = 1;
      Zp[k](l, 1) = expl(obs_thus_far + l);

      Xp[k](l, 0) = 1;
      for(int m = 1; m < num_fixed; m++)
      {
	Xp[k](l, m) = design(obs_thus_far + l, m);
      }
    }
 
    obs_thus_far += tp(k);
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
  
      for (int l = 0; l < tp(k); l++)
      	y[k](l) = X(i, obs_thus_far + l);
 
      obs_thus_far += tp(k);
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
        tracevar += (tp(k) - sigma2s * trace(Ws[k]));
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
  if(argc != 9)
  {
    cout << "Usage hypothesis-testing-tutorial lpts_paramfile design_file timepoints_file template" 
	 << " fixed random num_perm out_file"
	 << endl;
    return EXIT_FAILURE;
  }

  // reading lpts param file
  ifstream lpts_paramfile, design_file, timepoints_file, shapefile; 
  ofstream outfile; outfile.open(argv[8]);
  string filename;
  double px, py, pz, entry;
  unsigned int t;

  // loading lpts file
  unsigned int nsh = 0;
  lpts_paramfile.open(argv[1]);
  std::string line;
  while (std::getline(lpts_paramfile, line))
    ++nsh;
  lpts_paramfile.close();

  // std::cout << "nsh: " << nsh << "\n";

  // loading design file
  unsigned int num_fixed = 0;
  design_file.open(argv[2]);
  std::getline(design_file, line);
  std::istringstream iss(line);
  do
  {
    std::string sub;
    iss >> sub;
    if (sub.length())
      ++num_fixed;
  }
  while(iss);
  design_file.close();

  // outfile << "nfixed = " << num_fixed << endl;

  // counting number of subjects from timepts file
  unsigned int nsub = 0;
  timepoints_file.open(argv[3]);
  while (std::getline(timepoints_file, line))
    ++nsub;

  timepoints_file.close();

  // outfile << "nsub = " << nsub << endl;

  // counting number of particles
  unsigned int np = 0;
  lpts_paramfile.open(argv[1]);
  lpts_paramfile >> filename;
  shapefile.open(filename.c_str());
  while (std::getline(shapefile, line))
    ++np; 

  shapefile.close();
  lpts_paramfile.close();

  // outfile << "np = " << np << endl;

  // creating a subjects id for all shapes
  VectorType subj_ids; subj_ids.set_size(nsh);
  MatrixType subj_group; subj_group.zeros(nsub, 2);

  MatrixType X, Y; X.zeros(3*np, nsh); Y.zeros(3*np, nsh);
  MatrixType design; design.zeros(nsh, num_fixed);
  TimePointsType tp; tp.set_size(nsub);
  
  lpts_paramfile.open(argv[1]);

  unsigned int i = 0, j = 0, k = 0;
  while(lpts_paramfile >> filename)
  {
    shapefile.open(filename.c_str());
    outfile << "Reading file: " << filename << endl;
    i = 0;
    while(shapefile >> px >> py >> pz)
    {
      X(3 * i, j) = px;
      X(3 * i + 1, j) = py;
      X(3 * i + 2, j) = pz;
      ++i;
    }
    shapefile.close();
    ++j;
  }

  lpts_paramfile.close();

  // loading timepoints
  i = 0; int obs_thus_far = 0;
  timepoints_file.open(argv[3]);
  while(timepoints_file >> t)
  {    
    tp[i] = t;

    for(k = 0; k < tp[j]; k++)
      subj_ids[obs_thus_far + k] = i;

    subj_group(i, 1) = obs_thus_far;
    obs_thus_far += tp[i];
    ++i;
  }

  timepoints_file.close();
  
  // loading design matrix
  i = 0; j = 0; k = 0; obs_thus_far = tp[0];
  design_file.open(argv[2]);

  int int_col = 2, slope_col = 3;

  while(design_file >> entry)
  {
    design(j, k) = entry;
    ++k;
    if(k == num_fixed)
    {
      ++j;
      k = 0;
    }

    if(j == obs_thus_far)
    {
      subj_group(i, 0) = design(j-1, int_col);
      if(i < nsub - 1)
      {
	++i;
	obs_thus_far += tp[i];
      }
    }
  }

  design_file.close();

  // align shapes to template
  MatrixType template_shape; template_shape.zeros(3, np);
  // load in a nice template shape
  shapefile.open(argv[4]);
  j = 0;
  while(shapefile >> px >> py >> pz)
  {
    template_shape(0, j) = px;
    template_shape(1, j) = py;
    template_shape(2, j) = pz;
    ++j;
  }
  shapefile.close();

  // removing translations from template
  VectorType mean_xyz; mean_xyz.zeros(3); 
  for(k = 0; k < np; k++)
    mean_xyz += template_shape.col(k);

  mean_xyz /= static_cast<double>(np);

  for(k = 0; k < np; k++)
    template_shape.col(k) = template_shape.col(k) - mean_xyz;

  j = 0; // reinitializing

  // ...and align everything to the template
  for(j = 0; j < nsh; j++)
  {
    MatrixType current_shape = X.col(j);
    current_shape.reshape(3, np);

    // removing translations from current shape
    mean_xyz.zeros(3); 
    for(k = 0; k < np; k++)
      mean_xyz += current_shape.col(k);

    mean_xyz /= static_cast<double>(np);

    for(k = 0; k < np; k++)
      current_shape.col(k) = current_shape.col(k) - mean_xyz;

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

  // reading in fixed file
  shapefile.open(argv[5]);
  j = 0; k = 0;
  while(shapefile >> entry)
  {
    fixed(j, k) = entry;
    ++j;
    if(j == num_fixed)
    {
      ++k;
      j = 0;
    }
  }

  shapefile.close();

  // reading in random file
  shapefile.open(argv[6]);
  j = 0; k = 0;
  while(shapefile >> entry)
  {
    random(j, k) = entry;
    ++j;
    if(j == 2)
    {
      ++k;
      j = 0;
    }
  }

  shapefile.close();

  // ------------ Hypothesis testing: control vs hd --------------
  int numPermutes = atoi(argv[7]);

  double origTestStat_int = norm(fixed.row(int_col), "fro");
  double origTestStat_slope = norm(fixed.row(slope_col), "fro");
  outfile << "orig test stat = " << origTestStat_int << " " << origTestStat_slope << endl;

  uvec ctrl = find(subj_group.col(0) == 0);
  uvec hd = find(subj_group.col(0) == 1);

  Col<unsigned int> x; x.zeros(ctrl.size() + hd.size());
  for(int i = 0; i < x.size(); i++)
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

  for(int i = 0; i < numPermutes; i++)
  {
    permute(x);
    MatrixType new_design; new_design.zeros(nsh, num_fixed);
    new_design = design;

    for(j = 0; j < ctrl.size(); j++)
    {
      for(int k = 0; k < tp[x[j]]; k++)
      {
	new_design(subj_group(x[j], 1) + k, int_col) = 0;
	new_design(subj_group(x[j], 1) + k, slope_col) = 0;
      }
    }

    for(j = ctrl.size(); j < x.size(); j++)
    {
      for(int k = 0; k < tp[x[j]]; k++)
      {
	new_design(subj_group(x[j], 1) + k, int_col) = 1;
	new_design(subj_group(x[j], 1) + k, slope_col) = new_design(subj_group(x[j], 1) + k, 1);
      }
    }

    // outfile << "new design = " << endl << new_design << endl; 
  
    runTests(Y, new_design, tp, new_fixed, new_random);
    double testStat_int = norm(new_fixed.row(int_col), "fro");
    double testStat_slope = norm(new_fixed.row(slope_col), "fro");

    // intercept
    if(testStat_int < origTestStat_int)
      numBelow_int++;
    outfile << i << ": Test stat, orig (int) = " << testStat_int << "\t" << origTestStat_int << "\t"
	    << 1 - ((double) numBelow_int / (double)(i + 1)) << endl;

    // slope
    if(testStat_slope < origTestStat_slope)
      numBelow_slope++;
    outfile << i << ": Test stat, orig (slope) = " << testStat_slope << "\t" << origTestStat_slope << "\t"
	    << 1 - ((double) numBelow_slope / (double)(i + 1)) << endl;
  }

  outfile << "p value (int) = " << 1 - ((double)(numBelow_int) / numPermutes) << endl;
  outfile << "p value (slope) = " << 1 - ((double)(numBelow_slope) / numPermutes) << endl;
  outfile.close();

  return EXIT_SUCCESS;
}
