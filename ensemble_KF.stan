functions{
  real KalmanFilter_seq_em(vector rhos, vector dee, vector R, matrix Q, matrix C,// parameters
  matrix P, vector xhat, int Time, matrix y,matrix obs){
  real log_like = 0;
  matrix [num_elements(xhat),num_elements(xhat)] P_ = P;
  matrix [num_elements(xhat),num_elements(xhat)] A = rhos * rhos';
  vector [num_elements(xhat)] xhat_=xhat;
  matrix [num_elements(xhat),cols(y)] G;
  matrix [num_elements(xhat),cols(y)] G_;
  real er;
  real est;
  real Qs_inv;
  //matrix [cols(y),cols(y)] F;
  for (i in 1:Time)
  {
    P_ = P_ .* A + Q; // l 455
    xhat_ = rhos .* xhat_; // l 454
    for (j in 1:cols(y))
    {
      if (obs[i,j]==1.0)
      {
        G_[,j] = P_ * C[j,]';
        Qs_inv = 1 / (C[j,] * G_[,j] + R[j]);
        est = dot_product((xhat_ + dee)' ,  C[j,]) ; // careful we are required to make a transformation in dee space!! bigM * dee
        er = y[i,j] - est; 
        G[,j] = Qs_inv * G_[,j];
        xhat_ += G[,j] * (y[i,j] - est);
        P_ -=  G[,j] * G_[,j]';
        log_like -= 0.5 * (- log(Qs_inv) + er^2 * Qs_inv);
      }
    }
  }
  return log_like;
  }
  
  // future
  matrix KalmanFilter_back(vector rhos, vector dee, vector R, matrix Q, matrix C,// parameters
  matrix P, vector xhat, int Time, matrix y,matrix obs){
  matrix [Time,num_elements(xhat)] xhat_b;
  matrix [num_elements(xhat),num_elements(xhat)] P_ = P;
  matrix [num_elements(xhat),num_elements(xhat)] P_s [Time];
  matrix [num_elements(xhat),num_elements(xhat)] A = rhos * rhos';
  vector [num_elements(xhat)] xhat_=xhat;
  matrix [Time,num_elements(xhat)] xhat_s;
  matrix [num_elements(xhat),cols(y)] G[Time];
  matrix [num_elements(xhat),cols(y)] G_;
  matrix [Time,cols(y)] er;
  real est;
  matrix[Time,cols(y)] Qs_inv;
  /// smoothing bit
  matrix [num_elements(xhat),num_elements(xhat)] Identity = diag_matrix(rep_vector(1.0,num_elements(xhat)));
  matrix [num_elements(xhat),num_elements(xhat)] L;
  vector [num_elements(xhat)] r = rep_vector(0,num_elements(xhat));
  int k_i;
  int l_j;
  for (i in 1:Time)
  {
    P_ = P_ .* A + Q; // l 455
    xhat_ = rhos .* xhat_; // l 454
    P_s[i,,] = P_ ;
    xhat_s[i,] = xhat_' ;
    for (j in 1:cols(y)) // change this to some input that lists the numbers of the elements of y that we are interested in, maybe another value that says how many we are interested in.
    {
      if (obs[i,j]==1.0)
      {
        G_[,j] = P_ * C[j,]';
        Qs_inv[i,j] = 1 / (C[j,] * G_[,j] + R[j]);
        est = dot_product((xhat_ + dee)' ,  C[j,]) ;
        er[i,j] = y[i,j] - est; 
        G[i,,j] = Qs_inv[i,j] * G_[,j];
        xhat_ += G[i,,j] * er[i,j];
        P_ -=  G[i,,j] * G_[,j]';
      }
    }
  }
  for (i in 1:Time)
  {
    k_i = Time + 1 - i;
    for (j in 1:cols(y))
    {
      l_j = cols(y) + 1 - j;
      if (obs[k_i,l_j] == 1.0)
      {
        L = Identity - G[k_i,,l_j] * C[l_j,];
        r = C[l_j,]'  * Qs_inv[k_i,l_j] * er[k_i,l_j] + L' * r;
      }
    }
    xhat_b[k_i,] = xhat_s[k_i,] + (P_s[k_i,,] * r)';
    r = rhos .* r;
  }
  return xhat_b;
  }
}

data {
  int<lower=0> N; // years
  int<lower=0> M; // models
  matrix[N,M] M_runs;
  int N_os;
  int N_oe; // past years
  vector[N_oe - N_os + 1] obs;
  real obs_sd; // maybe we can learn this
  real l_disc_p; // priors
  real init_t_p; //
  real s_ldisc_p; // 
  real tr_shape;
  real tr_rate;
}

transformed data{
  matrix [M+1,M+2] M_mat = rep_matrix(0.0,M+1,M+2);
  matrix [N,M+1] which_obs = rep_matrix(1.0,N,M+1);
  matrix [N,M+1] y = rep_matrix(0.0,N,M+1);
  which_obs[,1] = rep_vector(0.0,N);
  for (i in N_os:N_oe){
    y[i,1] = obs[i - N_os + 1];
    which_obs[i,1] = 1.0;
  }
  y[,2:(M+1)] = M_runs;
  // this has to be automated
  M_mat[,1] = rep_vector(1.0,M+1);
  M_mat[,2] = rep_vector(1.0,M+1);
  M_mat[1,2] = 0.0;
  M_mat[2:(M+1),3:(M+2)] = diag_matrix(rep_vector(1.0,M));
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real shared_l_disc;
  //vector[N] shared_s_disc_raw;
  real<lower=0> shared_s_disc_sig_sq;
  real<lower=-1,upper=1> rho_s;
  real<lower=0> l_disc_var;
  vector[M] l_disc_raw;
  vector<lower=-1,upper=1>[M] rho_i;
  real<lower=0> sd_tr_sq;
  vector<lower=0>[M] sd_ind_sq;
}
transformed parameters{
  // 
  vector[M+2] A = rep_vector(0.0,M+2); // rhos
  vector[M+2] dee = rep_vector(0.0,M+2); //
  vector[M+1] R = rep_vector(0.0,M+1); // obs_var
  matrix[M + 2, M +2] SIGMA = rep_matrix(0.0,M+2,M+2);
  matrix[M+2, M+2] SIGMA_init = rep_matrix(0.0,M+2,M+2);
  vector[M] l_disc = l_disc_raw * sqrt(l_disc_var);
  //
  A = append_row([0.0,rho_s]',rho_i);
  dee = append_row([0.0,shared_l_disc]',l_disc);
  R = append_row(obs_sd^2,rep_vector(0.0,M));
  SIGMA = diag_matrix(append_row([sd_tr_sq,shared_s_disc_sig_sq]',sd_ind_sq));
  // C is going to be calculated in data -- matrix
  SIGMA_init = diag_matrix(append_row([init_t_p^2,shared_s_disc_sig_sq / (1.0 - rho_s^2.0)]',sd_ind_sq ./ (1.0 - rho_i .* rho_i)));
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  shared_s_disc_sig_sq ~ cauchy(0,10); // the lower the better
  sd_ind_sq ~ cauchy(0,10); // the lower the better
  sd_tr_sq ~ inv_gamma(tr_shape,tr_rate); // I have some knowledge of this
  shared_l_disc ~ normal(0.0,s_ldisc_p);
  l_disc_var ~ cauchy(0,10);
  l_disc_raw ~ normal(0.0,1.0); // test sensitvity of this
  //
  target += KalmanFilter_seq_em(A,dee, R, SIGMA, M_mat,SIGMA_init, rep_vector(0.0,M+2), N, y,which_obs);
}

generated quantities{
  //clever thing by Durbin and Koopman
  matrix [N,M+2] new_x;
  matrix [N,M+2] latent;
  matrix [N,M+2] latent_mean;//
  matrix [N,M+1] new_y; // new observations
  // initial parameters
  new_x[1,] = (multi_normal_rng(rep_vector(0.0,M+2), SIGMA_init))' ;
  for (t in 2:N){
    new_x[t,] = (multi_normal_rng(A .* new_x[t-1,]', SIGMA))';
  }
  for (t in 1:N){
    for (i in 1:(M+1)){
      if (R[i] > 0.0){
        new_y[t,i] = normal_rng(dot_product((new_x[t,]' + dee),M_mat[i,]), sqrt(R[i])); // need a dot product here? -- don't think so although it might help to make it a real
      }else{
        new_y[t,i] = dot_product((new_x[t,]' + dee),M_mat[i,]); // need to add a dee on here. Some problems because of that!! might need a dot product
      }
    }
  }
  latent = KalmanFilter_back(A,dee, R, SIGMA, M_mat,SIGMA_init, rep_vector(0.0,M+2), N, y,which_obs) - KalmanFilter_back(A,dee, R, SIGMA, M_mat,SIGMA_init, rep_vector(0.0,M+2), N, new_y,which_obs) + new_x;
  latent_mean = KalmanFilter_back(A,dee, R, SIGMA, M_mat,SIGMA_init, rep_vector(0.0,M+2), N, y,which_obs);
}
