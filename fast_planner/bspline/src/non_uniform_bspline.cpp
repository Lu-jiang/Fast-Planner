/**
* This file is part of Fast-Planner.
*
* Copyright 2019 Boyu Zhou, Aerial Robotics Group, Hong Kong University of Science and Technology, <uav.ust.hk>
* Developed by Boyu Zhou <bzhouai at connect dot ust dot hk>, <uv dot boyuzhou at gmail dot com>
* for more information see <https://github.com/HKUST-Aerial-Robotics/Fast-Planner>.
* If you use this code, please cite the respective publications as
* listed on the above website.
*
* Fast-Planner is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Fast-Planner is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with Fast-Planner. If not, see <http://www.gnu.org/licenses/>.
*/



#include "bspline/non_uniform_bspline.h"
#include <ros/ros.h>

namespace fast_planner {

NonUniformBspline::NonUniformBspline(const Eigen::MatrixXd& points, const int& order,
                                     const double& interval) {
  setUniformBspline(points, order, interval);
}

NonUniformBspline::~NonUniformBspline() {}


/*
输入参数为，点，阶数，时间间隔
功能：初始化成员变量：控制点矩阵、B样条次数、时间间隔、多样条阶数、节点向量
*/
void NonUniformBspline::setUniformBspline(const Eigen::MatrixXd& points, const int& order,
                                          const double& interval) {
  //control_points_是动态大小的矩阵，元素为双浮点数，每行代表一个控制点，维度由列号确定即N*3矩阵
  control_points_ = points;
  //p为B样条次数
  p_              = order;
  //节点之间的是检查，即delta_t
  interval_       = interval;
 //B样条控制点个数-1得到阶数
  n_ = points.rows() - 1;
  //节点数为阶数+次数
  m_ = n_ + p_ + 1;
//初始化存储节点的向量
  u_ = Eigen::VectorXd::Zero(m_ + 1);
  //B样条节点初始化，初始化为均匀B样条，即间隔相等。
  /*
  前p_+1个节点被称为开头节点或者起始节点，他们决定了B样条的起始点，这样定义出来的起始节点为小于0的区间
  中间的p_ + 1 到 m_ - p_ 之间节点为中间节点，由累加而成。确保连续
  后m_ - p_ 个节点称为尾节点（trailing nodes）或结束节点，它们决定了 B 样条曲线的结束点，初始化时，这些节点的值与前一个节点相等，确保在这一段区间内的节点是连续的。
  这个区间的长度是 (m_ - p_) 
  自己写出来很好理解
  */
  for (int i = 0; i <= m_; ++i) {

    if (i <= p_) {
      u_(i) = double(-p_ + i) * interval_;
    } else if (i > p_ && i <= m_ - p_) {
      u_(i) = u_(i - 1) + interval_;
    } else if (i > m_ - p_) {
      u_(i) = u_(i - 1) + interval_;
    }
  }
}

void NonUniformBspline::setKnot(const Eigen::VectorXd& knot) { this->u_ = knot; }

Eigen::VectorXd NonUniformBspline::getKnot() { return this->u_; }

void NonUniformBspline::getTimeSpan(double& um, double& um_p) {
  um   = u_(p_);
  um_p = u_(m_ - p_);
}

Eigen::MatrixXd NonUniformBspline::getControlPoint() { return control_points_; }

pair<Eigen::VectorXd, Eigen::VectorXd> NonUniformBspline::getHeadTailPts() {
  Eigen::VectorXd head = evaluateDeBoor(u_(p_));
  Eigen::VectorXd tail = evaluateDeBoor(u_(m_ - p_));
  return make_pair(head, tail);
}


/*De Boor递归算法，计算传入参数u处的B样条曲线值，返回B-样条曲线在参数 u 处的B样条曲线函数*/
Eigen::VectorXd NonUniformBspline::evaluateDeBoor(const double& u) {
//限制参数 u 的范围： 将参数 u 限制在 B-样条曲线有效范围内，即 [u_(p_), u_(m_ - p_)]。
  double ub = min(max(u_(p_), u), u_(m_ - p_));
 
  // determine which [ui,ui+1] lay in
  //找到输入的u是第几个节点处
  int k = p_;
  while (true) {
    if (u_(k + 1) >= ub) break;
    ++k;
  }
 
//通过k找到与当前节点相关的控制点k-p_到k
  /* deBoor's alg */
  vector<Eigen::VectorXd> d;
  for (int i = 0; i <= p_; ++i) {
    d.push_back(control_points_.row(k - p_ + i));
    // cout << d[i].transpose() << endl;
  }

//De Boor递归计算控制点，计算过程中使用了混合参数 alpha，该参数用于混合现有控制点以生成新的控制点。
  for (int r = 1; r <= p_; ++r) {
    for (int i = p_; i >= r; --i) {
      double alpha = (ub - u_[i + k - p_]) / (u_[i + 1 + k - r] - u_[i + k - p_]);
      // cout << "alpha: " << alpha << endl;
      d[i] = (1 - alpha) * d[i - 1] + alpha * d[i];
    }
  }

  return d[p_];
}

Eigen::VectorXd NonUniformBspline::evaluateDeBoorT(const double& t) {
  return evaluateDeBoor(t + u_(p_));
}

Eigen::MatrixXd NonUniformBspline::getDerivativeControlPoints() {
  // The derivative of a b-spline is also a b-spline, its order become p_-1
  // control point Qi = p_*(Pi+1-Pi)/(ui+p_+1-ui+1)
  Eigen::MatrixXd ctp = Eigen::MatrixXd::Zero(control_points_.rows() - 1, control_points_.cols());
  for (int i = 0; i < ctp.rows(); ++i) {
    ctp.row(i) =
        p_ * (control_points_.row(i + 1) - control_points_.row(i)) / (u_(i + p_ + 1) - u_(i + 1));
  }
  return ctp;
}

NonUniformBspline NonUniformBspline::getDerivative() {
  Eigen::MatrixXd   ctp = getDerivativeControlPoints();
  NonUniformBspline derivative(ctp, p_ - 1, interval_);

  /* cut the first and last knot */
  Eigen::VectorXd knot(u_.rows() - 2);
  knot = u_.segment(1, u_.rows() - 2);
  derivative.setKnot(knot);

  return derivative;
}

double NonUniformBspline::getInterval() { return interval_; }

void NonUniformBspline::setPhysicalLimits(const double& vel, const double& acc) {
  limit_vel_   = vel;
  limit_acc_   = acc;
  limit_ratio_ = 1.1;
}

/*轨迹安全性检查，并更新超阈值比例*/
bool NonUniformBspline::checkFeasibility(bool show) {
  bool fea = true;//可行性标记位
  // SETY << "[Bspline]: total points size: " << control_points_.rows() << endl;
 
  Eigen::MatrixXd P         = control_points_;//取出控制点
  int             dimension = control_points_.cols();//求控制点列数即维度，x，y，z
 
  /* check vel feasibility and insert points */
  double max_vel = -1.0;//初始化最大速度为负
  for (int i = 0; i < P.rows() - 1; ++i) {//遍历每一个控制点
    Eigen::VectorXd vel = p_ * (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1));//根据控制点计算速度
 
    if (fabs(vel(0)) > limit_vel_ + 1e-4 || fabs(vel(1)) > limit_vel_ + 1e-4 ||//判断是否有任一维度上的速度超过阈值
        fabs(vel(2)) > limit_vel_ + 1e-4) {
 
      if (show) cout << "[Check]: Infeasible vel " << i << " :" << vel.transpose() << endl;
      fea = false;//标记为不可行
 
      for (int j = 0; j < dimension; ++j) {
        max_vel = max(max_vel, fabs(vel(j)));//采用计算出来的速度更新最大速度
      }
    }
  }
 
  /* acc feasibility */
  double max_acc = -1.0;
  for (int i = 0; i < P.rows() - 2; ++i) {//遍历n-2个点
 
    Eigen::VectorXd acc = p_ * (p_ - 1) *
        ((P.row(i + 2) - P.row(i + 1)) / (u_(i + p_ + 2) - u_(i + 2)) -
         (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1))) /
        (u_(i + p_ + 1) - u_(i + 2));//计算加速度
 
    if (fabs(acc(0)) > limit_acc_ + 1e-4 || fabs(acc(1)) > limit_acc_ + 1e-4 ||
        fabs(acc(2)) > limit_acc_ + 1e-4) {//判断是否有任一维度上的加速度超过阈值
 
      if (show) cout << "[Check]: Infeasible acc " << i << " :" << acc.transpose() << endl;
      fea = false;
 
      for (int j = 0; j < dimension; ++j) {//记录更新最大加速度
        max_acc = max(max_acc, fabs(acc(j)));
      }
    }
  }
 
  double ratio = max(max_vel / limit_vel_, sqrt(fabs(max_acc) / limit_acc_));//计算最大速度或加速度超过约束的
 
  return fea;
}

double NonUniformBspline::checkRatio() {
  Eigen::MatrixXd P         = control_points_;
  int             dimension = control_points_.cols();

  // find max vel
  double max_vel = -1.0;
  for (int i = 0; i < P.rows() - 1; ++i) {
    Eigen::VectorXd vel = p_ * (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1));
    for (int j = 0; j < dimension; ++j) {
      max_vel = max(max_vel, fabs(vel(j)));
    }
  }
  // find max acc
  double max_acc = -1.0;
  for (int i = 0; i < P.rows() - 2; ++i) {
    Eigen::VectorXd acc = p_ * (p_ - 1) *
        ((P.row(i + 2) - P.row(i + 1)) / (u_(i + p_ + 2) - u_(i + 2)) -
         (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1))) /
        (u_(i + p_ + 1) - u_(i + 2));
    for (int j = 0; j < dimension; ++j) {
      max_acc = max(max_acc, fabs(acc(j)));
    }
  }
  double ratio = max(max_vel / limit_vel_, sqrt(fabs(max_acc) / limit_acc_));
  ROS_ERROR_COND(ratio > 2.0, "max vel: %lf, max acc: %lf.", max_vel, max_acc);

  return ratio;
}

bool NonUniformBspline::reallocateTime(bool show) {
  // SETY << "[Bspline]: total points size: " << control_points_.rows() << endl;
  // cout << "origin knots:\n" << u_.transpose() << endl;
  bool fea = true;

  Eigen::MatrixXd P         = control_points_;
  int             dimension = control_points_.cols();

  double max_vel, max_acc;

  /* check vel feasibility and insert points */
  for (int i = 0; i < P.rows() - 1; ++i) {
    Eigen::VectorXd vel = p_ * (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1));

    if (fabs(vel(0)) > limit_vel_ + 1e-4 || fabs(vel(1)) > limit_vel_ + 1e-4 ||
        fabs(vel(2)) > limit_vel_ + 1e-4) {

      fea = false;
      if (show) cout << "[Realloc]: Infeasible vel " << i << " :" << vel.transpose() << endl;

      max_vel = -1.0;
      for (int j = 0; j < dimension; ++j) {
        max_vel = max(max_vel, fabs(vel(j)));
      }

      double ratio = max_vel / limit_vel_ + 1e-4;
      if (ratio > limit_ratio_) ratio = limit_ratio_;

      double time_ori = u_(i + p_ + 1) - u_(i + 1);
      double time_new = ratio * time_ori;
      double delta_t  = time_new - time_ori;
      double t_inc    = delta_t / double(p_);

      for (int j = i + 2; j <= i + p_ + 1; ++j) {
        u_(j) += double(j - i - 1) * t_inc;
        if (j <= 5 && j >= 1) {
          // cout << "vel j: " << j << endl;
        }
      }

      for (int j = i + p_ + 2; j < u_.rows(); ++j) {
        u_(j) += delta_t;
      }
    }
  }

  /* acc feasibility */
  for (int i = 0; i < P.rows() - 2; ++i) {

    Eigen::VectorXd acc = p_ * (p_ - 1) *
        ((P.row(i + 2) - P.row(i + 1)) / (u_(i + p_ + 2) - u_(i + 2)) -
         (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1))) /
        (u_(i + p_ + 1) - u_(i + 2));

    if (fabs(acc(0)) > limit_acc_ + 1e-4 || fabs(acc(1)) > limit_acc_ + 1e-4 ||
        fabs(acc(2)) > limit_acc_ + 1e-4) {

      fea = false;
      if (show) cout << "[Realloc]: Infeasible acc " << i << " :" << acc.transpose() << endl;

      max_acc = -1.0;
      for (int j = 0; j < dimension; ++j) {
        max_acc = max(max_acc, fabs(acc(j)));
      }

      double ratio = sqrt(max_acc / limit_acc_) + 1e-4;
      if (ratio > limit_ratio_) ratio = limit_ratio_;
      // cout << "ratio: " << ratio << endl;

      double time_ori = u_(i + p_ + 1) - u_(i + 2);
      double time_new = ratio * time_ori;
      double delta_t  = time_new - time_ori;
      double t_inc    = delta_t / double(p_ - 1);

      if (i == 1 || i == 2) {
        // cout << "acc i: " << i << endl;
        for (int j = 2; j <= 5; ++j) {
          u_(j) += double(j - 1) * t_inc;
        }

        for (int j = 6; j < u_.rows(); ++j) {
          u_(j) += 4.0 * t_inc;
        }

      } else {

        for (int j = i + 3; j <= i + p_ + 1; ++j) {
          u_(j) += double(j - i - 2) * t_inc;
          if (j <= 5 && j >= 1) {
            // cout << "acc j: " << j << endl;
          }
        }

        for (int j = i + p_ + 2; j < u_.rows(); ++j) {
          u_(j) += delta_t;
        }
      }
    }
  }

  return fea;
}

void NonUniformBspline::lengthenTime(const double& ratio) {
  int num1 = 5;
  int num2 = getKnot().rows() - 1 - 5;

  double delta_t = (ratio - 1.0) * (u_(num2) - u_(num1));
  double t_inc   = delta_t / double(num2 - num1);
  for (int i = num1 + 1; i <= num2; ++i) u_(i) += double(i - num1) * t_inc;
  for (int i = num2 + 1; i < u_.rows(); ++i) u_(i) += delta_t;
}

void NonUniformBspline::recomputeInit() {}

/*
函数参数：ts：轨迹执行时间、point_set：原轨迹点集合、start_end_derivative：起点与终点的高阶约束
输出：更新ctrl_pts控制点矩阵的值
函数功能：将给定点集和起始/终止导数转换为 B-样条曲线的控制点矩阵，通过对原始轨迹的拟合得到B样条轨迹的控制点
*/
void NonUniformBspline::parameterizeToBspline(const double& ts, const vector<Eigen::Vector3d>& point_set,
                                              const vector<Eigen::Vector3d>& start_end_derivative,
                                              Eigen::MatrixXd&               ctrl_pts) {
  //判断输入的时间是否合理
  if (ts <= 0) {
    cout << "[B-spline]:time step error." << endl;
    return;
  }
  //判断输入的轨迹点是否合理
  if (point_set.size() < 2) {
    cout << "[B-spline]:point set have only " << point_set.size() << " points." << endl;
    return;
  }
//起点与终点的导数限制（即一阶和二阶导数，有此信息才可以确保得到一个具有足够信息的线性方程组，从而可以求解出 B-样条曲线的控制点。）
  if (start_end_derivative.size() != 4) {
    cout << "[B-spline]:derivatives error." << endl;
  }
//记录原曲线路经点个数
  int K = point_set.size();
 
  // write A
  //该矩阵用来构建线性方程组， B-样条曲线参数化过程中求解控制点
  //一阶 B-样条基函数的系数向量、一阶 B-样条基函数的导数系数向量、二阶 B-样条基函数的系数向量
  //取该数值的原因是B样条曲线矩阵表示论文中提出的，以满足 B-样条曲线的一些良好性质。
  Eigen::Vector3d prow(3), vrow(3), arow(3);
  prow << 1, 4, 1;    // 位置基函数的系数（对应“曲线位置”与控制点的关系）
  vrow << -1, 0, 1;   // 速度基函数的系数（对应“曲线速度”与控制点的关系）
  arow << 1, -2, 1;   // 加速度基函数的系数（对应“曲线加速度”与控制点的关系）
 
//K+4是因为添加了四行导数约束信息
  // 行数（K+4）：拟合目标的数量 = 原始路径点数量（K） + 导数约束数量（4，起点速度、终点速度、起点加速度、终点加速度）；
  // 列数（K+2）：控制点的数量 = 原始路径点数量（K） + 2（3 次 B 样条需要比路径点多 2 个控制点才能覆盖所有段，数学上的必要条件）。
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(K + 4, K + 2);
 
  for (int i = 0; i < K; ++i) A.block(i, i, 1, 3) = (1 / 6.0) * prow.transpose();
 
  A.block(K, 0, 1, 3)         = (1 / 2.0 / ts) * vrow.transpose();
  A.block(K + 1, K - 1, 1, 3) = (1 / 2.0 / ts) * vrow.transpose();
 
  A.block(K + 2, 0, 1, 3)     = (1 / ts / ts) * arow.transpose();
  A.block(K + 3, K - 1, 1, 3) = (1 / ts / ts) * arow.transpose();
  // cout << "A:\n" << A << endl;
 
  // A.block(0, 0, K, K + 2) = (1 / 6.0) * A.block(0, 0, K, K + 2);
  // A.block(K, 0, 2, K + 2) = (1 / 2.0 / ts) * A.block(K, 0, 2, K + 2);
  // A.row(K + 4) = (1 / ts / ts) * A.row(K + 4);
  // A.row(K + 5) = (1 / ts / ts) * A.row(K + 5);
 
  // write b
  // 初始化x、y、z三个方向的b向量
  Eigen::VectorXd bx(K + 4), by(K + 4), bz(K + 4);

  // 前K个元素：原始路径点的x/y/z坐标（拟合位置目标）
  for (int i = 0; i < K; ++i) {
    bx(i) = point_set[i](0);  // 第i个路径点的x坐标
    by(i) = point_set[i](1);  // 第i个路径点的y坐标
    bz(i) = point_set[i](2);  // 第i个路径点的z坐标
  }

  // 后4个元素：起点/终点的导数约束（速度、加速度的x/y/z分量）
  for (int i = 0; i < 4; ++i) {
    bx(K + i) = start_end_derivative[i](0);  // 第i个导数约束的x分量
    by(K + i) = start_end_derivative[i](1);  // 第i个导数约束的y分量
    bz(K + i) = start_end_derivative[i](2);  // 第i个导数约束的z分量
  }
 
  // solve Ax = b
  //求解线性方程
  // 用ColPivHouseholderQr分解法求解线性方程组（数值稳定，适合中等维度矩阵）
  Eigen::VectorXd px = A.colPivHouseholderQr().solve(bx);  // x方向控制点
  Eigen::VectorXd py = A.colPivHouseholderQr().solve(by);  // y方向控制点
  Eigen::VectorXd pz = A.colPivHouseholderQr().solve(bz);  // z方向控制点
 
  // convert to control pts
  //控制点赋值
  ctrl_pts.resize(K + 2, 3);  // 控制点数量是K+2，每个控制点有3个坐标（x/y/z）
  ctrl_pts.col(0) = px;       // 第0列：所有控制点的x坐标
  ctrl_pts.col(1) = py;       // 第1列：所有控制点的y坐标
  ctrl_pts.col(2) = pz;       // 第2列：所有控制点的z坐标
 
  // cout << "[B-spline]: parameterization ok." << endl;
}

double NonUniformBspline::getTimeSum() {
  double tm, tmp;
  getTimeSpan(tm, tmp);
  return tmp - tm;
}

double NonUniformBspline::getLength(const double& res) {
  double          length = 0.0;
  double          dur    = getTimeSum();
  Eigen::VectorXd p_l    = evaluateDeBoorT(0.0), p_n;
  for (double t = res; t <= dur + 1e-4; t += res) {
    p_n = evaluateDeBoorT(t);
    length += (p_n - p_l).norm();
    p_l = p_n;
  }
  return length;
}

double NonUniformBspline::getJerk() {
  NonUniformBspline jerk_traj = getDerivative().getDerivative().getDerivative();

  Eigen::VectorXd times     = jerk_traj.getKnot();
  Eigen::MatrixXd ctrl_pts  = jerk_traj.getControlPoint();
  int             dimension = ctrl_pts.cols();

  double jerk = 0.0;
  for (int i = 0; i < ctrl_pts.rows(); ++i) {
    for (int j = 0; j < dimension; ++j) {
      jerk += (times(i + 1) - times(i)) * ctrl_pts(i, j) * ctrl_pts(i, j);
    }
  }

  return jerk;
}

void NonUniformBspline::getMeanAndMaxVel(double& mean_v, double& max_v) {
  NonUniformBspline vel = getDerivative();
  double            tm, tmp;
  vel.getTimeSpan(tm, tmp);

  double max_vel = -1.0, mean_vel = 0.0;
  int    num = 0;
  for (double t = tm; t <= tmp; t += 0.01) {
    Eigen::VectorXd vxd = vel.evaluateDeBoor(t);
    double          vn  = vxd.norm();

    mean_vel += vn;
    ++num;
    if (vn > max_vel) {
      max_vel = vn;
    }
  }

  mean_vel = mean_vel / double(num);
  mean_v   = mean_vel;
  max_v    = max_vel;
}

void NonUniformBspline::getMeanAndMaxAcc(double& mean_a, double& max_a) {
  NonUniformBspline acc = getDerivative().getDerivative();
  double            tm, tmp;
  acc.getTimeSpan(tm, tmp);

  double max_acc = -1.0, mean_acc = 0.0;
  int    num = 0;
  for (double t = tm; t <= tmp; t += 0.01) {
    Eigen::VectorXd axd = acc.evaluateDeBoor(t);
    double          an  = axd.norm();

    mean_acc += an;
    ++num;
    if (an > max_acc) {
      max_acc = an;
    }
  }

  mean_acc = mean_acc / double(num);
  mean_a   = mean_acc;
  max_a    = max_acc;
}
}  // namespace fast_planner
