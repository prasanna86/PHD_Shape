#include "primitives.hxx"
#include <fstream>
#include <iostream>

MatrixType procrustes(MatrixType & p, MatrixType & q)
{
  MatrixType m = p * strans(q);

  MatrixType U, V; VectorType s;
  svd_econ(U, s, V, m, 'b');

  MatrixType rotation = U * strans(V);
  MatrixType qRot = rotation * q;
  return qRot;
}

void runTests(MatrixType & X, MatrixType & design, TimePointsType & tp, MatrixType & fixed, MatrixType & random)
{
  cout << "Estimating params using ME" << endl;
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
  double epsilon = 1.0e-05;
  double sigmaDiff = 2 * epsilon, Ds_norm_diff = 2 * epsilon;
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

    while(em < 500 && (sigmaDiff > epsilon || Ds_norm_diff > epsilon)) // EM iterations
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

    sigmaDiff = 2 * epsilon;
    Ds_norm_diff = 2 * epsilon;

    cout << "particle, em iter, Ds_new norm, sigma2 = " << i+1 << " " << em << " " << norm(Ds_new, "fro") << " " << sigma2s_new << endl;
  } //endfor all points on shape (x,y & z)
  
  delete [] Vs;
  delete [] Ws;
  delete [] identity_n;
  delete [] Xp;
  delete [] Zp;
  delete [] y;
  delete [] residual;
}

void print_stuff(MatrixType & design, MatrixType & fixed, MatrixType & random)
{
  int np = fixed.n_cols / 3, j;
  int nsub = random.n_cols / (3 * np);

  // output fixed params to file
  double min_expl = min(design.col(1));
  double max_expl = max(design.col(1));
  int step = 10;
  
  VectorType control_female_shapes; control_female_shapes.set_size(np * 3);
  VectorType control_male_shapes; control_male_shapes.set_size(np * 3);
  VectorType low_female_shapes; low_female_shapes.set_size(np * 3);
  VectorType low_male_shapes; low_male_shapes.set_size(np * 3);
  VectorType med_female_shapes; med_female_shapes.set_size(np * 3);
  VectorType med_male_shapes; med_male_shapes.set_size(np * 3);
  VectorType high_female_shapes; high_female_shapes.set_size(np * 3);
  VectorType high_male_shapes; high_male_shapes.set_size(np * 3);

  VectorType control_female_slope; control_female_slope.set_size(np * 3);
  VectorType control_male_slope; control_male_slope.set_size(np * 3);
  VectorType low_female_slope; low_female_slope.set_size(np * 3);
  VectorType low_male_slope; low_male_slope.set_size(np * 3);
  VectorType med_female_slope; med_female_slope.set_size(np * 3);
  VectorType med_male_slope; med_male_slope.set_size(np * 3);
  VectorType high_female_slope; high_female_slope.set_size(np * 3);
  VectorType high_male_slope; high_male_slope.set_size(np * 3);

  string fname;
  ofstream out;
  int counter = 0;
  
  control_female_slope = strans(fixed.row(1));
  low_female_slope = strans(fixed.row(1) + fixed.row(5));
  med_female_slope = strans(fixed.row(1) + fixed.row(7)); 
  high_female_slope = strans(fixed.row(1) + fixed.row(9)); 

  control_male_slope = strans(fixed.row(1) + fixed.row(3));
  low_male_slope = strans(fixed.row(1) + fixed.row(3) + fixed.row(5)); 
  med_male_slope = strans(fixed.row(1) + fixed.row(3) + fixed.row(7)); 
  high_male_slope = strans(fixed.row(1) + fixed.row(3) + fixed.row(9)); 
  
  for(int i = 0; i < step; i++)
  {
    double value = min_expl + static_cast<double>(i) * ((max_expl - min_expl) / static_cast<double>(step - 1));
    cout << value << endl;

    control_female_shapes = strans(fixed.row(0) + value * fixed.row(1));
    low_female_shapes = strans((fixed.row(0) + fixed.row(4)) + value * (fixed.row(1) + fixed.row(5)));
    med_female_shapes = strans((fixed.row(0) + fixed.row(6)) + value * (fixed.row(1) + fixed.row(7))); 
    high_female_shapes = strans((fixed.row(0) + fixed.row(8)) + value * (fixed.row(1) + fixed.row(9))); 

    control_male_shapes = strans((fixed.row(0) + fixed.row(2)) + value * (fixed.row(1) + fixed.row(3)));
    low_male_shapes = strans((fixed.row(0) + fixed.row(2) + fixed.row(4)) + value * (fixed.row(1) + fixed.row(3) + fixed.row(5))); 
    med_male_shapes = strans((fixed.row(0) + fixed.row(2) + fixed.row(6)) + value * (fixed.row(1) + fixed.row(3) + fixed.row(7))); 
    high_male_shapes = strans((fixed.row(0) + fixed.row(2) + fixed.row(8)) + value * (fixed.row(1) + fixed.row(3) + fixed.row(9))); 
    
    // printing control_female_shapes
    ostringstream ss;
    ss << i;
    fname = "control_female" + ss.str() + ".lpts";
    out.open(fname.c_str());
    counter = 0;
    for(j = 0; j < np * 3; j++)
    {
      //cout << control_female_shapes(j) << ' ';
      out << control_female_shapes(j) << ' ';
      counter++;
      if(counter == 3)
      {
        //cout << endl;
	out << endl;
	counter = 0;
      }
    }
    out.close();

    // printing low_female
    fname = "low_female" + ss.str() + ".lpts";
    out.open(fname.c_str());
    counter = 0;
    for(j = 0; j < np * 3; j++)
    {
      out << low_female_shapes[j] << ' ';
      counter++;
      if(counter == 3)
      {
	out << endl;
	counter = 0;
      }
    }

    out.close();

    // printing med_female
    fname = "med_female" + ss.str() + ".lpts";
    out.open(fname.c_str());
    counter = 0;
    for(j = 0; j < np * 3; j++)
    {
      out << med_female_shapes[j] << ' ';
      counter++;
      if(counter == 3)
      {
	out << endl;
	counter = 0;
      }
    }

    out.close();

    // printing high_female
    fname = "high_female" + ss.str() + ".lpts";
    out.open(fname.c_str());
    counter = 0;
    for(j = 0; j < np * 3; j++)
    {
      out << high_female_shapes[j] << ' ';
      counter++;
      if(counter == 3)
      {
	out << endl;
	counter = 0;
      }
    }

    out.close();

    // printing control_male
    fname = "control_male" + ss.str() + ".lpts";
    out.open(fname.c_str());
    counter = 0;
    for(j = 0; j < np * 3; j++)
    {
      out << control_male_shapes[j] << ' ';
      counter++;
      if(counter == 3)
      {
	out << endl;
	counter = 0;
      }
    }

    out.close();

    // printing low_male
    fname = "low_male" + ss.str() + ".lpts";
    out.open(fname.c_str());
    counter = 0;
    for(j = 0; j < np * 3; j++)
    {
      out << low_male_shapes[j] << ' ';
      counter++;
      if(counter == 3)
      {
	out << endl;
	counter = 0;
      }
    }

    out.close();

    // printing med_male
    fname = "med_male" + ss.str() + ".lpts";
    out.open(fname.c_str());
    counter = 0;
    for(j = 0; j < np * 3; j++)
    {
      out << med_male_shapes[j] << ' ';
      counter++;
      if(counter == 3)
      {
	out << endl;
	counter = 0;
      }
    }

    out.close();

    // printing high_male
    fname = "high_male" + ss.str() + ".lpts";
    out.open(fname.c_str());
    counter = 0;
    for(j = 0; j < np * 3; j++)
    {
      out << high_male_shapes[j] << ' ';
      counter++;
      if(counter == 3)
      {
	out << endl;
	counter = 0;
      }
    }
    
    out.close();
  }

  // printing slopes:

  // control_female_slope
  fname = "control_female_slope.txt";
  out.open(fname.c_str());
  counter = 0;
  for(j = 0; j < np * 3; j++)
  {
    out << control_female_slope[j] << ' ';
    counter++;
    if(counter == 3)
    {
      out << endl;
      counter = 0;
    }
  }

  out.close();

  // low_female_slope
  fname = "low_female_slope.txt";
  out.open(fname.c_str());
  counter = 0;
  for(j = 0; j < np * 3; j++)
  {
    out << low_female_slope[j] << ' ';
    counter++;
    if(counter == 3)
    {
      out << endl;
      counter = 0;
    }
  }

  out.close();

  // printing med_female_slope
  fname = "med_female_slope.txt";
  out.open(fname.c_str());
  counter = 0;
  for(j = 0; j < np * 3; j++)
  {
    out << med_female_slope[j] << ' ';
    counter++;
    if(counter == 3)
    {
      out << endl;
      counter = 0;
    }
  }

  out.close();

  // printing high_female_slope
  fname = "high_female_slope.txt";
  out.open(fname.c_str());
  counter = 0;
  for(j = 0; j < np * 3; j++)
  {
    out << high_female_slope[j] << ' ';
    counter++;
    if(counter == 3)
    {
      out << endl;
      counter = 0;
    }
  }

  out.close();

  // printing control_male_slope
  fname = "control_male_slope.txt";
  out.open(fname.c_str());
  counter = 0;
  for(j = 0; j < np * 3; j++)
  {
    out << control_male_slope[j] << ' ';
    counter++;
    if(counter == 3)
    {
      out << endl;
      counter = 0;
    }
  }

  out.close();

  // printing low_male_slope
  fname = "low_male_slope.txt";
  out.open(fname.c_str());
  counter = 0;
  for(j = 0; j < np * 3; j++)
  {
    out << low_male_slope[j] << ' ';
    counter++;
    if(counter == 3)
    {
      out << endl;
      counter = 0;
    }
  }

  out.close();

  // printing med_male_slope
  fname = "med_male_slope.txt";
  out.open(fname.c_str());
  counter = 0;
  for(j = 0; j < np * 3; j++)
  {
    out << med_male_slope[j] << ' ';
    counter++;
    if(counter == 3)
    {
      out << endl;
      counter = 0;
    }
  }

  out.close();

  // printing high_male_slope
  fname = "high_male_slope.txt";
  out.open(fname.c_str());
  counter = 0;
  for(j = 0; j < np * 3; j++)
  {
    out << high_male_slope[j] << ' ';
    counter++;
    if(counter == 3)
    {
      out << endl;
      counter = 0;
    }
  }
    
  out.close();
}

int main(int argc, char ** argv)
{
  if(argc != 5)
  {
    cout << "Usage mixed-effects-estimation lpts_paramfile design_file timepoints_file print_stuff (yes: 1, no: 0)"
	 << endl;
    return EXIT_FAILURE;
  }

  // reading lpts param file
  ifstream lpts_paramfile, design_file, timepoints_file, shapefile;
  string filename;
  unsigned int np = 642, nsh = 49, nsub = 14, num_fixed = 4, num_cat = 8;
  unsigned int j = 0;
  double px, py, pz, age, g1, g2, g3, sex;
  unsigned int t;
  MatrixType X, Y; X.zeros(3*np, nsh); Y.zeros(3*np, nsh);
  MatrixType design; design.zeros(nsh, num_fixed);
  TimePointsType tp; tp.set_size(nsub);

  // loading pts file
  lpts_paramfile.open(argv[1]);
  while(lpts_paramfile >> filename)
  {
    shapefile.open(filename.c_str());
    cout << "Reading file: " << filename << endl;
    int i = 0;
    while(shapefile >> px >> py >> pz)
    {
      X(3 * i, j) = px;
      X(3 * i + 1, j) = py;
      X(3 * i + 2, j) = pz;
      i += 1;
    }
    shapefile.close();
    j += 1;
  }

  // loading design matrix
  j = 0;
  design_file.open(argv[2]);
  while(design_file >> sex >> g1 >> g2 >> g3 >> age)
  {
    design(j, 0) = 1;
    design(j, 1) = age;
    design(j, 2) = sex;
    design(j, 3) = age * sex;
    design(j, 4) = g1;
    design(j, 5) = age * g1;
    design(j, 6) = g2;
    design(j, 7) = age * g2;
    design(j, 8) = g3;
    design(j, 9) = age * g3;

    j += 1;
  }

  // loading timepoints
  j = 0;
  timepoints_file.open(argv[3]);
  while(timepoints_file >> t)
  {    
    tp[j] = t;
    j += 1;
  }

  j = 0;

  // parameters of interest
  MatrixType fixed; // slopes + intercepts for all points
  MatrixType random; // slopes + intercepts for all groups, for all points
  fixed.set_size(num_fixed, 3 * np);
  random.set_size(2, 3 * np * nsub);

  // align shapes to template
  MatrixType template_shape; template_shape.zeros(3, np);
  // load in a nice template shape (control female)...
  shapefile.open("template_control_female.lpts");
  j = 0;
  while(shapefile >> px >> py >> pz)
  {
    template_shape(0, j) = px;
    template_shape(1, j) = py;
    template_shape(2, j) = pz;
    j += 1;
  }
  
  // removing translations from template
  VectorType mean_xyz; mean_xyz.zeros(3); 
  for(int k = 0; k < np; k++)
    mean_xyz += template_shape.col(k);

  mean_xyz /= static_cast<double>(np);

  for(int k = 0; k < np; k++)
    template_shape.col(k) = template_shape.col(k) - mean_xyz;

  j = 0; // reinitializing
  shapefile.close();

  // ...and align everything to the template
  for(int k = 0; k < nsh; k++)
  {
    MatrixType current_shape = X.col(k);
    current_shape.reshape(3, np);

    // removing translations from current shape
    mean_xyz.zeros(3); 
    for(int k = 0; k < np; k++)
      mean_xyz += current_shape.col(k);

    mean_xyz /= static_cast<double>(np);

    for(int k = 0; k < np; k++)
      current_shape.col(k) = current_shape.col(k) - mean_xyz;

    // procrustes align to template
    current_shape = procrustes(template_shape, current_shape);
    
    // reshape back
    current_shape.reshape(3 * np, 1);
    Y.col(k) = current_shape;
  }
    
  // estimate model with aligned shapes
  runTests(Y, design, tp, fixed, random);

  if(atoi(argv[4]) == 1)
    print_stuff(design, fixed, random);

  return EXIT_SUCCESS;
}
