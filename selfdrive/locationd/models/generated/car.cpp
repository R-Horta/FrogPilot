#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4709938673122520820) {
   out_4709938673122520820[0] = delta_x[0] + nom_x[0];
   out_4709938673122520820[1] = delta_x[1] + nom_x[1];
   out_4709938673122520820[2] = delta_x[2] + nom_x[2];
   out_4709938673122520820[3] = delta_x[3] + nom_x[3];
   out_4709938673122520820[4] = delta_x[4] + nom_x[4];
   out_4709938673122520820[5] = delta_x[5] + nom_x[5];
   out_4709938673122520820[6] = delta_x[6] + nom_x[6];
   out_4709938673122520820[7] = delta_x[7] + nom_x[7];
   out_4709938673122520820[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_5966965373853949397) {
   out_5966965373853949397[0] = -nom_x[0] + true_x[0];
   out_5966965373853949397[1] = -nom_x[1] + true_x[1];
   out_5966965373853949397[2] = -nom_x[2] + true_x[2];
   out_5966965373853949397[3] = -nom_x[3] + true_x[3];
   out_5966965373853949397[4] = -nom_x[4] + true_x[4];
   out_5966965373853949397[5] = -nom_x[5] + true_x[5];
   out_5966965373853949397[6] = -nom_x[6] + true_x[6];
   out_5966965373853949397[7] = -nom_x[7] + true_x[7];
   out_5966965373853949397[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_6615511429899936093) {
   out_6615511429899936093[0] = 1.0;
   out_6615511429899936093[1] = 0;
   out_6615511429899936093[2] = 0;
   out_6615511429899936093[3] = 0;
   out_6615511429899936093[4] = 0;
   out_6615511429899936093[5] = 0;
   out_6615511429899936093[6] = 0;
   out_6615511429899936093[7] = 0;
   out_6615511429899936093[8] = 0;
   out_6615511429899936093[9] = 0;
   out_6615511429899936093[10] = 1.0;
   out_6615511429899936093[11] = 0;
   out_6615511429899936093[12] = 0;
   out_6615511429899936093[13] = 0;
   out_6615511429899936093[14] = 0;
   out_6615511429899936093[15] = 0;
   out_6615511429899936093[16] = 0;
   out_6615511429899936093[17] = 0;
   out_6615511429899936093[18] = 0;
   out_6615511429899936093[19] = 0;
   out_6615511429899936093[20] = 1.0;
   out_6615511429899936093[21] = 0;
   out_6615511429899936093[22] = 0;
   out_6615511429899936093[23] = 0;
   out_6615511429899936093[24] = 0;
   out_6615511429899936093[25] = 0;
   out_6615511429899936093[26] = 0;
   out_6615511429899936093[27] = 0;
   out_6615511429899936093[28] = 0;
   out_6615511429899936093[29] = 0;
   out_6615511429899936093[30] = 1.0;
   out_6615511429899936093[31] = 0;
   out_6615511429899936093[32] = 0;
   out_6615511429899936093[33] = 0;
   out_6615511429899936093[34] = 0;
   out_6615511429899936093[35] = 0;
   out_6615511429899936093[36] = 0;
   out_6615511429899936093[37] = 0;
   out_6615511429899936093[38] = 0;
   out_6615511429899936093[39] = 0;
   out_6615511429899936093[40] = 1.0;
   out_6615511429899936093[41] = 0;
   out_6615511429899936093[42] = 0;
   out_6615511429899936093[43] = 0;
   out_6615511429899936093[44] = 0;
   out_6615511429899936093[45] = 0;
   out_6615511429899936093[46] = 0;
   out_6615511429899936093[47] = 0;
   out_6615511429899936093[48] = 0;
   out_6615511429899936093[49] = 0;
   out_6615511429899936093[50] = 1.0;
   out_6615511429899936093[51] = 0;
   out_6615511429899936093[52] = 0;
   out_6615511429899936093[53] = 0;
   out_6615511429899936093[54] = 0;
   out_6615511429899936093[55] = 0;
   out_6615511429899936093[56] = 0;
   out_6615511429899936093[57] = 0;
   out_6615511429899936093[58] = 0;
   out_6615511429899936093[59] = 0;
   out_6615511429899936093[60] = 1.0;
   out_6615511429899936093[61] = 0;
   out_6615511429899936093[62] = 0;
   out_6615511429899936093[63] = 0;
   out_6615511429899936093[64] = 0;
   out_6615511429899936093[65] = 0;
   out_6615511429899936093[66] = 0;
   out_6615511429899936093[67] = 0;
   out_6615511429899936093[68] = 0;
   out_6615511429899936093[69] = 0;
   out_6615511429899936093[70] = 1.0;
   out_6615511429899936093[71] = 0;
   out_6615511429899936093[72] = 0;
   out_6615511429899936093[73] = 0;
   out_6615511429899936093[74] = 0;
   out_6615511429899936093[75] = 0;
   out_6615511429899936093[76] = 0;
   out_6615511429899936093[77] = 0;
   out_6615511429899936093[78] = 0;
   out_6615511429899936093[79] = 0;
   out_6615511429899936093[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_416114564767500801) {
   out_416114564767500801[0] = state[0];
   out_416114564767500801[1] = state[1];
   out_416114564767500801[2] = state[2];
   out_416114564767500801[3] = state[3];
   out_416114564767500801[4] = state[4];
   out_416114564767500801[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_416114564767500801[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_416114564767500801[7] = state[7];
   out_416114564767500801[8] = state[8];
}
void F_fun(double *state, double dt, double *out_6271263489731407948) {
   out_6271263489731407948[0] = 1;
   out_6271263489731407948[1] = 0;
   out_6271263489731407948[2] = 0;
   out_6271263489731407948[3] = 0;
   out_6271263489731407948[4] = 0;
   out_6271263489731407948[5] = 0;
   out_6271263489731407948[6] = 0;
   out_6271263489731407948[7] = 0;
   out_6271263489731407948[8] = 0;
   out_6271263489731407948[9] = 0;
   out_6271263489731407948[10] = 1;
   out_6271263489731407948[11] = 0;
   out_6271263489731407948[12] = 0;
   out_6271263489731407948[13] = 0;
   out_6271263489731407948[14] = 0;
   out_6271263489731407948[15] = 0;
   out_6271263489731407948[16] = 0;
   out_6271263489731407948[17] = 0;
   out_6271263489731407948[18] = 0;
   out_6271263489731407948[19] = 0;
   out_6271263489731407948[20] = 1;
   out_6271263489731407948[21] = 0;
   out_6271263489731407948[22] = 0;
   out_6271263489731407948[23] = 0;
   out_6271263489731407948[24] = 0;
   out_6271263489731407948[25] = 0;
   out_6271263489731407948[26] = 0;
   out_6271263489731407948[27] = 0;
   out_6271263489731407948[28] = 0;
   out_6271263489731407948[29] = 0;
   out_6271263489731407948[30] = 1;
   out_6271263489731407948[31] = 0;
   out_6271263489731407948[32] = 0;
   out_6271263489731407948[33] = 0;
   out_6271263489731407948[34] = 0;
   out_6271263489731407948[35] = 0;
   out_6271263489731407948[36] = 0;
   out_6271263489731407948[37] = 0;
   out_6271263489731407948[38] = 0;
   out_6271263489731407948[39] = 0;
   out_6271263489731407948[40] = 1;
   out_6271263489731407948[41] = 0;
   out_6271263489731407948[42] = 0;
   out_6271263489731407948[43] = 0;
   out_6271263489731407948[44] = 0;
   out_6271263489731407948[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_6271263489731407948[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_6271263489731407948[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6271263489731407948[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6271263489731407948[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_6271263489731407948[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_6271263489731407948[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_6271263489731407948[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_6271263489731407948[53] = -9.8000000000000007*dt;
   out_6271263489731407948[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_6271263489731407948[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_6271263489731407948[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6271263489731407948[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6271263489731407948[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_6271263489731407948[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_6271263489731407948[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_6271263489731407948[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6271263489731407948[62] = 0;
   out_6271263489731407948[63] = 0;
   out_6271263489731407948[64] = 0;
   out_6271263489731407948[65] = 0;
   out_6271263489731407948[66] = 0;
   out_6271263489731407948[67] = 0;
   out_6271263489731407948[68] = 0;
   out_6271263489731407948[69] = 0;
   out_6271263489731407948[70] = 1;
   out_6271263489731407948[71] = 0;
   out_6271263489731407948[72] = 0;
   out_6271263489731407948[73] = 0;
   out_6271263489731407948[74] = 0;
   out_6271263489731407948[75] = 0;
   out_6271263489731407948[76] = 0;
   out_6271263489731407948[77] = 0;
   out_6271263489731407948[78] = 0;
   out_6271263489731407948[79] = 0;
   out_6271263489731407948[80] = 1;
}
void h_25(double *state, double *unused, double *out_5463292229834969125) {
   out_5463292229834969125[0] = state[6];
}
void H_25(double *state, double *unused, double *out_6598408891029782616) {
   out_6598408891029782616[0] = 0;
   out_6598408891029782616[1] = 0;
   out_6598408891029782616[2] = 0;
   out_6598408891029782616[3] = 0;
   out_6598408891029782616[4] = 0;
   out_6598408891029782616[5] = 0;
   out_6598408891029782616[6] = 1;
   out_6598408891029782616[7] = 0;
   out_6598408891029782616[8] = 0;
}
void h_24(double *state, double *unused, double *out_1453363659542421936) {
   out_1453363659542421936[0] = state[4];
   out_1453363659542421936[1] = state[5];
}
void H_24(double *state, double *unused, double *out_6995054195383642111) {
   out_6995054195383642111[0] = 0;
   out_6995054195383642111[1] = 0;
   out_6995054195383642111[2] = 0;
   out_6995054195383642111[3] = 0;
   out_6995054195383642111[4] = 1;
   out_6995054195383642111[5] = 0;
   out_6995054195383642111[6] = 0;
   out_6995054195383642111[7] = 0;
   out_6995054195383642111[8] = 0;
   out_6995054195383642111[9] = 0;
   out_6995054195383642111[10] = 0;
   out_6995054195383642111[11] = 0;
   out_6995054195383642111[12] = 0;
   out_6995054195383642111[13] = 0;
   out_6995054195383642111[14] = 1;
   out_6995054195383642111[15] = 0;
   out_6995054195383642111[16] = 0;
   out_6995054195383642111[17] = 0;
}
void h_30(double *state, double *unused, double *out_5188098167550463236) {
   out_5188098167550463236[0] = state[4];
}
void H_30(double *state, double *unused, double *out_318281450461834139) {
   out_318281450461834139[0] = 0;
   out_318281450461834139[1] = 0;
   out_318281450461834139[2] = 0;
   out_318281450461834139[3] = 0;
   out_318281450461834139[4] = 1;
   out_318281450461834139[5] = 0;
   out_318281450461834139[6] = 0;
   out_318281450461834139[7] = 0;
   out_318281450461834139[8] = 0;
}
void h_26(double *state, double *unused, double *out_8256972638498242735) {
   out_8256972638498242735[0] = state[7];
}
void H_26(double *state, double *unused, double *out_8106831863805712776) {
   out_8106831863805712776[0] = 0;
   out_8106831863805712776[1] = 0;
   out_8106831863805712776[2] = 0;
   out_8106831863805712776[3] = 0;
   out_8106831863805712776[4] = 0;
   out_8106831863805712776[5] = 0;
   out_8106831863805712776[6] = 0;
   out_8106831863805712776[7] = 1;
   out_8106831863805712776[8] = 0;
}
void h_27(double *state, double *unused, double *out_6406802082736631560) {
   out_6406802082736631560[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1856481861338590772) {
   out_1856481861338590772[0] = 0;
   out_1856481861338590772[1] = 0;
   out_1856481861338590772[2] = 0;
   out_1856481861338590772[3] = 1;
   out_1856481861338590772[4] = 0;
   out_1856481861338590772[5] = 0;
   out_1856481861338590772[6] = 0;
   out_1856481861338590772[7] = 0;
   out_1856481861338590772[8] = 0;
}
void h_29(double *state, double *unused, double *out_186424899889617591) {
   out_186424899889617591[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3569844588208141805) {
   out_3569844588208141805[0] = 0;
   out_3569844588208141805[1] = 1;
   out_3569844588208141805[2] = 0;
   out_3569844588208141805[3] = 0;
   out_3569844588208141805[4] = 0;
   out_3569844588208141805[5] = 0;
   out_3569844588208141805[6] = 0;
   out_3569844588208141805[7] = 0;
   out_3569844588208141805[8] = 0;
}
void h_28(double *state, double *unused, double *out_5640900473559070910) {
   out_5640900473559070910[0] = state[0];
}
void H_28(double *state, double *unused, double *out_8652243605277672379) {
   out_8652243605277672379[0] = 1;
   out_8652243605277672379[1] = 0;
   out_8652243605277672379[2] = 0;
   out_8652243605277672379[3] = 0;
   out_8652243605277672379[4] = 0;
   out_8652243605277672379[5] = 0;
   out_8652243605277672379[6] = 0;
   out_8652243605277672379[7] = 0;
   out_8652243605277672379[8] = 0;
}
void h_31(double *state, double *unused, double *out_2911814669133737213) {
   out_2911814669133737213[0] = state[8];
}
void H_31(double *state, double *unused, double *out_6567762929152822188) {
   out_6567762929152822188[0] = 0;
   out_6567762929152822188[1] = 0;
   out_6567762929152822188[2] = 0;
   out_6567762929152822188[3] = 0;
   out_6567762929152822188[4] = 0;
   out_6567762929152822188[5] = 0;
   out_6567762929152822188[6] = 0;
   out_6567762929152822188[7] = 0;
   out_6567762929152822188[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_4709938673122520820) {
  err_fun(nom_x, delta_x, out_4709938673122520820);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5966965373853949397) {
  inv_err_fun(nom_x, true_x, out_5966965373853949397);
}
void car_H_mod_fun(double *state, double *out_6615511429899936093) {
  H_mod_fun(state, out_6615511429899936093);
}
void car_f_fun(double *state, double dt, double *out_416114564767500801) {
  f_fun(state,  dt, out_416114564767500801);
}
void car_F_fun(double *state, double dt, double *out_6271263489731407948) {
  F_fun(state,  dt, out_6271263489731407948);
}
void car_h_25(double *state, double *unused, double *out_5463292229834969125) {
  h_25(state, unused, out_5463292229834969125);
}
void car_H_25(double *state, double *unused, double *out_6598408891029782616) {
  H_25(state, unused, out_6598408891029782616);
}
void car_h_24(double *state, double *unused, double *out_1453363659542421936) {
  h_24(state, unused, out_1453363659542421936);
}
void car_H_24(double *state, double *unused, double *out_6995054195383642111) {
  H_24(state, unused, out_6995054195383642111);
}
void car_h_30(double *state, double *unused, double *out_5188098167550463236) {
  h_30(state, unused, out_5188098167550463236);
}
void car_H_30(double *state, double *unused, double *out_318281450461834139) {
  H_30(state, unused, out_318281450461834139);
}
void car_h_26(double *state, double *unused, double *out_8256972638498242735) {
  h_26(state, unused, out_8256972638498242735);
}
void car_H_26(double *state, double *unused, double *out_8106831863805712776) {
  H_26(state, unused, out_8106831863805712776);
}
void car_h_27(double *state, double *unused, double *out_6406802082736631560) {
  h_27(state, unused, out_6406802082736631560);
}
void car_H_27(double *state, double *unused, double *out_1856481861338590772) {
  H_27(state, unused, out_1856481861338590772);
}
void car_h_29(double *state, double *unused, double *out_186424899889617591) {
  h_29(state, unused, out_186424899889617591);
}
void car_H_29(double *state, double *unused, double *out_3569844588208141805) {
  H_29(state, unused, out_3569844588208141805);
}
void car_h_28(double *state, double *unused, double *out_5640900473559070910) {
  h_28(state, unused, out_5640900473559070910);
}
void car_H_28(double *state, double *unused, double *out_8652243605277672379) {
  H_28(state, unused, out_8652243605277672379);
}
void car_h_31(double *state, double *unused, double *out_2911814669133737213) {
  h_31(state, unused, out_2911814669133737213);
}
void car_H_31(double *state, double *unused, double *out_6567762929152822188) {
  H_31(state, unused, out_6567762929152822188);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
