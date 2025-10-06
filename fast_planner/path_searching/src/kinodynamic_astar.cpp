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

#include <path_searching/kinodynamic_astar.h>
#include <sstream>
#include <plan_env/sdf_map.h>

using namespace std;
using namespace Eigen;

namespace fast_planner
{
KinodynamicAstar::~KinodynamicAstar()
{
  for (int i = 0; i < allocate_num_; i++)
  {
    delete path_node_pool_[i];
  }
}

int KinodynamicAstar::search(Eigen::Vector3d start_pt, Eigen::Vector3d start_v, Eigen::Vector3d start_a,
                             Eigen::Vector3d end_pt, Eigen::Vector3d end_v, 
                             bool init, bool dynamic, double time_start)
{
  //函数传入参数为起始点的位置、速度、加速度、终点位置、速度、初始化标志位、动态（可行）标志位？、起始时间
  start_vel_ = start_v;//取出传入的起始点的速度
  start_acc_ = start_a;//传入起始点的加速度
  //path_node_pool_ 可能在初始化时被分配一定数量的节点，这些节点在运行时被重复使用，而不是每次需要节点时都动态地分配内存。这有助于提高性能，避免频繁的内存分配和释放。
  PathNodePtr cur_node = path_node_pool_[0];//取出第一个路径点赋给当前节点
  cur_node->parent = NULL;//父节点
  cur_node->state.head(3) = start_pt;//state矩阵前3列记录位置
  cur_node->state.tail(3) = start_v;//state矩阵后三列记录速度
  cur_node->index = posToIndex(start_pt);//获取全剧坐标系下的位置索引
  cur_node->g_score = 0.0;//记录当前点成本代价值
 
  Eigen::VectorXd end_state(6);//初始化目标点状态的
  Eigen::Vector3i end_index;//记录终点索引值
  double time_to_goal;//路径规划时间
 
  end_state.head(3) = end_pt;
  end_state.tail(3) = end_v;
  end_index = posToIndex(end_pt);
  cur_node->f_score = lambda_heu_ * estimateHeuristic(cur_node->state, end_state, time_to_goal);//记录当前节点的代价函数
  cur_node->node_state = IN_OPEN_SET;//标记为待探索列表
  open_set_.push(cur_node);//将当前节点添加到openlist中
  use_node_num_ += 1;//已探索点个数记录

  if (dynamic)
  {
    time_origin_ = time_start;
    cur_node->time = time_start;
    cur_node->time_idx = timeToIndex(time_start);
    expanded_nodes_.insert(cur_node->index, cur_node->time_idx, cur_node);
    // cout << "time start: " << time_start << endl;
  }
  else
    expanded_nodes_.insert(cur_node->index, cur_node);

  PathNodePtr neighbor = NULL;
  PathNodePtr terminate_node = NULL;
  bool init_search = init;
  const int tolerance = ceil(1 / resolution_);

  while (!open_set_.empty())
  {
    cur_node = open_set_.top();


    // 1. 终止条件
    // Terminate?
    /*标志位在这个搜索算法中的主要作用是为了提前终止搜索。它表示当前节点是否已经离起点足够远，
    超过了预定的 horizon_ 阈值。这个设计可能是为了在搜索空间较大时，避免过度的搜索，以提高算法的效率。
    ，当 reach_horizon 为真时，意味着当前节点已经离起点足够远，可能是因为在搜索空间中没有更好的路径，
    或者当前路径不太可能变得更好。在这种情况下，算法可以提前结束搜索，以减少计算时间。
    */

    bool reach_horizon = (cur_node->state.head(3) - start_pt).norm() >= horizon_;
    bool near_end = abs(cur_node->index(0) - end_index(0)) <= tolerance &&
                    abs(cur_node->index(1) - end_index(1)) <= tolerance &&
                    abs(cur_node->index(2) - end_index(2)) <= tolerance;

    if (reach_horizon || near_end)
    {
      terminate_node = cur_node;
      retrievePath(terminate_node);
      if (near_end)
      {
        // Check whether shot traj exist
        estimateHeuristic(cur_node->state, end_state, time_to_goal);    // //计算优化的最小时间
        computeShotTraj(cur_node->state, end_state, time_to_goal);      // 带入优化的时间判断轨迹是否合理
        if (init_search)
          ROS_ERROR("Shot in first search loop!");
      }
    }
    if (reach_horizon)  // //路径动态不可行（规划失败）
    {
      if (is_shot_succ_) //搜索到的路径安全可行的
      {
        std::cout << "reach end" << std::endl;  
        // //此处终点是最远点（其实也不算是规划成功，因为超出距离后到达的是最远处路径可行，并非是原来要找的终点）
        return REACH_END;
      }
      else
      {
        std::cout << "reach horizon" << std::endl;  // //路径动态不可行（规划失败）
        return REACH_HORIZON;
      }
    }

    if (near_end)
    {
      if (is_shot_succ_)
      {
        std::cout << "reach end" << std::endl;
        return REACH_END;
      }
      else if (cur_node->parent != NULL)
      {
        std::cout << "near end" << std::endl;   //规划路径动态不可行，但是上一节点是可行的，此时为规划到接近终点的状态
        return NEAR_END;
      }
      else
      {
        std::cout << "no path" << std::endl;
        return NO_PATH;
      }
    }
    
     /*
  上述代码皆是单个循环内单个openlist节点内部的终止的判定条件，类似于递归代码的逻辑
  下方代码可以看作是递归单个节点不满足逻辑时候需要处理的步骤（反正我是这么理解的）
  */ 
    
    // 2. 节点扩张
    open_set_.pop();    // 弹出该节点
    cur_node->node_state = IN_CLOSE_SET;    //将已经探索的节点加入到close_list中，标记为已探索
    iter_num_ += 1;     // 搜索完的点+1

    double res = 1 / 2.0, time_res = 1 / 1.0, time_res_init = 1 / 20.0; //初始化时间常数
    Eigen::Matrix<double, 6, 1> cur_state = cur_node->state;    // xyz vx vy vz
    Eigen::Matrix<double, 6, 1> pro_state;  // 下一个节点状态
    vector<PathNodePtr> tmp_expand_nodes;   // 临时存扩展节点列表
    Eigen::Vector3d um;                     // //离散的控制量，这里指的是三维上加速度
    double pro_t;                           // 拓展节点的时间戳
    vector<Eigen::Vector3d> inputs;//输入
    vector<double> durations;//单位输入控制量的持续时间
    /*
    这里大家要首先明确一点，init_search标志位用来标志是否用不连续的初始状态重试搜索。
    本人的理解是如果有初始的节点输入状态存储，但是由于路径规划失败导致无法到达终点
    所以需要拓展新的节点重新规划路径。
    */

    if (init_search)    // //具有上次搜索失败的离散状态量
    {
      inputs.push_back(start_acc_); // //存储起始状态的加速度
      for (double tau = time_res_init * init_max_tau_; tau <= init_max_tau_ + 1e-3;
           tau += time_res_init * init_max_tau_)
        durations.push_back(tau);
      init_search = false;
    }
    else    //没有初始离散状态，未初始化成功，通过最大最小加速度初始化离散加速度后再离散时间
    {
      for (double ax = -max_acc_; ax <= max_acc_ + 1e-3; ax += max_acc_ * res)
        for (double ay = -max_acc_; ay <= max_acc_ + 1e-3; ay += max_acc_ * res)
          for (double az = -max_acc_; az <= max_acc_ + 1e-3; az += max_acc_ * res)
          {
            um << ax, ay, az;
            inputs.push_back(um);
          }
      for (double tau = time_res * max_tau_; tau <= max_tau_; tau += time_res * max_tau_)
        durations.push_back(tau);
    }

    // cout << "cur state:" << cur_state.head(3).transpose() << endl;
    //组合不同的时间与加速度组合，得到不一样的下一状态集合
    for (int i = 0; i < inputs.size(); ++i)
      for (int j = 0; j < durations.size(); ++j)
      {
        um = inputs[i];
        double tau = durations[j];
        stateTransit(cur_state, pro_state, um, tau);    //输入时间与加速度离散值，更新pro_state
        pro_t = cur_node->time + tau;   //拓展状态时间记录

        Eigen::Vector3d pro_pos = pro_state.head(3);    //拓展位置向量

        // Check if in close set
        Eigen::Vector3i pro_id = posToIndex(pro_pos);   // //拓展位置世界坐标系下的索引
        int pro_t_id = timeToIndex(pro_t);              // //拓展位置时间id

        ////根据是否带时间戳记录拓展节点
        PathNodePtr pro_node = dynamic ? expanded_nodes_.find(pro_id, pro_t_id) : expanded_nodes_.find(pro_id);
        if (pro_node != NULL && pro_node->node_state == IN_CLOSE_SET)   // //当存在于closet中，跳过此次拓展
        {
          if (init_search)
            std::cout << "close" << std::endl;
          continue;
        }

        // 此时，//剩下NULL或不存在与closelist中的点
        // Check maximal velocity
        Eigen::Vector3d pro_v = pro_state.tail(3);  //检查速度限制超过则跳过该次检查
        if (fabs(pro_v(0)) > max_vel_ || fabs(pro_v(1)) > max_vel_ || fabs(pro_v(2)) > max_vel_)
        {
          if (init_search)
            std::cout << "vel" << std::endl;
          continue;
        }

        // 此时， //剩下不超过速度限制的点，且不存在于close中，但可能为空
        // Check not in the same voxel
        Eigen::Vector3i diff = pro_id - cur_node->index;    //检查拓展的栅格是否为当前栅格相同
        int diff_time = pro_t_id - cur_node->time_idx;
        if (diff.norm() == 0 && ((!dynamic) || diff_time == 0)) // 带时间戳时时间相同，或者时间相同跳出循环
        {
          if (init_search)
            std::cout << "same" << std::endl;
          continue;
        }

        // 剩下不超过速度限制、且不是重复的点，且不存在于close中，但可能为空或者存在同一栅格中（但花费的时间不同的点）
        // Check safety
        Eigen::Vector3d pos;
        Eigen::Matrix<double, 6, 1> xt;
        bool is_occ = false;
        for (int k = 1; k <= check_num_; ++k)   //分辨率提高检查点，比如a到b我检查5个点，判定每个点是否占用
        {
          double dt = tau * double(k) / double(check_num_);
          stateTransit(cur_state, xt, um, dt);  //前向积分得到扩展节点的位置和速度
          pos = xt.head(3);
          if (edt_environment_->sdf_map_->getInflateOccupancy(pos) == 1 )   // //检查是否占用
          {
            is_occ = true;  // //占用标志位记为tru
            break;
          }
        }
        if (is_occ)
        {
          if (init_search)
            std::cout << "safe" << std::endl;   //当检测到被占用时，跳过
          continue;
        }

        //更新拓展状态的代价
        double time_to_goal, tmp_g_score, tmp_f_score;
        tmp_g_score = (um.squaredNorm() + w_time_) * tau + cur_node->g_score;
        tmp_f_score = tmp_g_score + lambda_heu_ * estimateHeuristic(pro_state, end_state, time_to_goal);



        // 3. 剪枝与添加
              /*
                prune为剪枝标志位
                如果 prune 为 true，表示需要对当前生成的节点进行剪枝。这说明当前节点与之前扩展的节点在相同的体素内（位置索引相同且时间索引相同，如果是动态规划）。
                如果当前节点的代价（tmp_f_score）比之前扩展的相同节点的代价更小，则更新之前扩展的节点的信息，以保留更优的路径信息。
                如果 prune 为 false，表示不需要对当前生成的节点进行剪枝。这说明当前节点与之前扩展的节点在不同的体素内。
                如果当前节点是一个新的节点（pro_node == NULL），将其添加到搜索队列中，并更新相关信息。
                如果当前节点已经在开放集合中，检查新路径的代价是否更小，如果是，则更新节点信息。
                */
        // Compare nodes expanded from the same parent
        bool prune = false;
        //遍历所有临时拓展节点列表，用来检查拓展的点中有没有同个节点内的点，有的话标记剪枝操作
        for (int j = 0; j < tmp_expand_nodes.size(); ++j)
        {
          PathNodePtr expand_node = tmp_expand_nodes[j];
          //当前遍历节点与本次循环拓展节点在相同的体素内
          if ((pro_id - expand_node->index).norm() == 0 && ((!dynamic) || pro_t_id == expand_node->time_idx))
          {
            prune = true;   //标记需要剪枝（即比较代价函数）
            if (tmp_f_score < expand_node->f_score) //比较代价函数，选取代价小的节点并保存
            {
              expand_node->f_score = tmp_f_score;
              expand_node->g_score = tmp_g_score;
              expand_node->state = pro_state;
              expand_node->input = um;
              expand_node->duration = tau;
              if (dynamic)
                expand_node->time = cur_node->time + tau;
            }
            break;
          }
        }
        
        // This node end up in a voxel different from others
        if (!prune) //当前节点为全新节点
        {
          if (pro_node == NULL) //如果 pro_node 为空，说明这是一个新的节点，需要将其添加到搜索队列中。更新节点的信息，包括索引、状态、代价等。
          {
            pro_node = path_node_pool_[use_node_num_];
            pro_node->index = pro_id;
            pro_node->state = pro_state;
            pro_node->f_score = tmp_f_score;
            pro_node->g_score = tmp_g_score;
            pro_node->input = um;
            pro_node->duration = tau;
            pro_node->parent = cur_node;
            pro_node->node_state = IN_OPEN_SET;
            if (dynamic)
            {
              pro_node->time = cur_node->time + tau;
              pro_node->time_idx = timeToIndex(pro_node->time);
            }
            open_set_.push(pro_node);

            if (dynamic)
              expanded_nodes_.insert(pro_id, pro_node->time, pro_node);
            else
              expanded_nodes_.insert(pro_id, pro_node);

            //加入临时探索列表中
            tmp_expand_nodes.push_back(pro_node);

            use_node_num_ += 1;
            if (use_node_num_ == allocate_num_) ////超过原定义好的ptr存储序列
            {
              cout << "run out of memory." << endl;
              return NO_PATH;
            }
          }
          else if (pro_node->node_state == IN_OPEN_SET) //如果 pro_node 已经在开放集合中，检查新路径的代价是否更小，如果是，则更新节点信息。
          {
            if (tmp_g_score < pro_node->g_score)
            {
              // pro_node->index = pro_id;
              pro_node->state = pro_state;
              pro_node->f_score = tmp_f_score;
              pro_node->g_score = tmp_g_score;
              pro_node->input = um;
              pro_node->duration = tau;
              pro_node->parent = cur_node;
              if (dynamic)
                pro_node->time = cur_node->time + tau;
            }
          }
          else
          {
            cout << "error type in searching: " << pro_node->node_state << endl;
          }
        }
      }
    // init_search = false;
  }

  cout << "open set empty, no path!" << endl;
  cout << "use node num: " << use_node_num_ << endl;
  cout << "iter num: " << iter_num_ << endl;
  return NO_PATH;
}

void KinodynamicAstar::setParam(ros::NodeHandle& nh)
{
  nh.param("search/max_tau", max_tau_, -1.0);
  nh.param("search/init_max_tau", init_max_tau_, -1.0);
  nh.param("search/max_vel", max_vel_, -1.0);
  nh.param("search/max_acc", max_acc_, -1.0);
  nh.param("search/w_time", w_time_, -1.0);
  nh.param("search/horizon", horizon_, -1.0);
  nh.param("search/resolution_astar", resolution_, -1.0);
  nh.param("search/time_resolution", time_resolution_, -1.0);
  nh.param("search/lambda_heu", lambda_heu_, -1.0);
  nh.param("search/allocate_num", allocate_num_, -1);
  nh.param("search/check_num", check_num_, -1);
  nh.param("search/optimistic", optimistic_, true);
  tie_breaker_ = 1.0 + 1.0 / 10000;

  double vel_margin;
  nh.param("search/vel_margin", vel_margin, 0.0);
  max_vel_ += vel_margin;
}

void KinodynamicAstar::retrievePath(PathNodePtr end_node)
{
  PathNodePtr cur_node = end_node;
  path_nodes_.push_back(cur_node);

  while (cur_node->parent != NULL)
  {
    cur_node = cur_node->parent;
    path_nodes_.push_back(cur_node);
  }

  reverse(path_nodes_.begin(), path_nodes_.end());
}
double KinodynamicAstar::estimateHeuristic(Eigen::VectorXd x1, Eigen::VectorXd x2, double& optimal_time)
{
  // estimateHeuristic用于计算两点之间的启发函数代价，传入参数为起点与终点的状态
  // 返回值为计算所得的代价值并且同时更新到达终点的最优时间optimal_time
  // 计算两个机器人状态（x1 为起点状态，x2 为终点状态）之间的启发式代价（即 “预估的最小代价”）。
  // Kinodynamic A* 和普通 A * 的核心区别是：
  //    它会考虑机器人的动力学约束（比如速度、加速度限制），而不是只算 “几何距离”。
  //    这段代码正是通过 “假设一条满足动力学的最优轨迹（这里是多项式轨迹）”，来估算两点间的最小代价,
  //    同时输出这条轨迹的 “最优时间”。

  // 1. 提取状态中的 “位置” 和 “速度”（基础数据准备）
  const Vector3d dp = x2.head(3) - x1.head(3);//位置变化量
  // head(size)：提取前 size 个元素（等价于segment(0, size)）
  // tail(size)：提取后 size 个元素（等价于segment(总长度 - size, size)）

  const Vector3d v0 = x1.segment(3, 3);//起点速度矩阵
  const Vector3d v1 = x2.segment(3, 3);//终点速度矩阵
  // segment()的作用是 “从第start个元素开始，提取连续size个元素，组成一个子向量”。


  // 2. 定义四次多项式的系数（核心：假设 “最优轨迹是多项式轨迹”）
  // Kinodynamic A* 中，启发式代价的计算通常会 “假设一条满足动力学的理想轨迹”（比如多项式轨迹，因为它能平滑衔接速度 / 加速度）。
  // 这里选择的是四次多项式轨迹（时间 t 为自变量，位置是 t 的四次函数），代码先计算该多项式对应的 “代价函数系数”：
  // c1~c5是四次多项式轨迹的代价函数系数（推导来自轨迹优化的数学模型）
  double c1 = -36 * dp.dot(dp);          // dp.dot(dp)：位置变化量的平方（即空间距离的平方）
  double c2 = 24 * (v0 + v1).dot(dp);    // (v0+v1)与dp的点积（关联速度和位置变化的关系）
  double c3 = -4 * (v0.dot(v0) + v0.dot(v1) + v1.dot(v1));  // 起点/终点速度的平方与点积（动力学约束项）
  double c4 = 0;                         // 四次项系数为0（实际是简化后的结果，不影响核心逻辑）
  double c5 = w_time_;                   // w_time_：“时间权重”（超参数，控制“时间代价”在总代价中的重要性）
  

  // 3. 求解 “最优时间 t” 的候选值（找到代价最小的轨迹时间）
  // 轨迹的 “总代价” 是时间 t 的函数（时间越长，时间成本越高；但时间太短，可能需要更大的加速度，能量成本越高）。
  // 我们需要找到 “让总代价最小的 t”，这一步就是求解 t 的候选值： 
  // 3.1. 调用quartic函数，求解四次方程 c5*t^4 + c4*t^3 + c3*t^2 + c2*t + c1 = 0 的根
  // （物理意义：代价函数对t求导后为0的点，即“代价最小/最大的候选时间”）
  std::vector<double> ts = quartic(c5, c4, c3, c2, c1);

  // 3.2. 计算“最短时间限制t_bar”（物理意义：不考虑加速度，只按最大速度走的最短时间）
  double v_max = max_vel_ * 0.5;  // max_vel_是机器人最大速度，乘以0.5是留安全余量
  // lpNorm<Infinity>()：无穷范数（即x/y/z三个方向中位置变化最大的那个值，保证不超过最大速度）
  double t_bar = (x1.head(3) - x2.head(3)).lpNorm<Infinity>() / v_max;
  ts.push_back(t_bar);  // 把“最短时间”也加入候选集（防止数学解超出物理极限）


  // 4. 从候选时间中筛选 “最优解”（计算最小代价）
  // 初始化：cost设为极大值（确保第一次比较能更新），t_d是最终的最优时间
  double cost = 100000000;
  double t_d = t_bar;

  // 遍历所有候选时间t
  for (auto t : ts)
  {
    if (t < t_bar)  // 时间小于“物理最短时间”，违反速度约束，跳过
      continue;
  
    // 计算当前t对应的总代价（推导自四次多项式轨迹的代价函数，已化简）
    // 代价由两部分组成：1. 能量相关项（-c1/(3t³) -c2/(2t²) -c3/t）；2. 时间成本（w_time_ * t）
    double c = -c1 / (3 * t * t * t) - c2 / (2 * t * t) - c3 / t + w_time_ * t;
  
    if (c < cost)  // 找到更小的代价，更新最优值
    {
      cost = c;
      t_d = t;
    }
  }

  // 输出最优时间（通过引用参数optimal_time返回）
  optimal_time = t_d;

  // 返回启发式代价：1.0*(1+tie_breaker_)*cost
  // tie_breaker_：“打破平局”的小系数（防止两个节点代价相同，影响A*的排序效率）
  return 1.0 * (1 + tie_breaker_) * cost;
  // 总结：通过 “假设一条连接 x1 和 x2 的四次多项式轨迹”，
  // 求解 “满足机器人最大速度约束、且总代价（能量 + 时间）最小” 的轨迹参数，
  // 将这个最小代价作为 Kinodynamic A * 的启发式值，同时输出这条最优轨迹的时间。
}

bool KinodynamicAstar::computeShotTraj(Eigen::VectorXd state1, Eigen::VectorXd state2, double time_to_goal)
{
  /* ---------- get coefficient ---------- */
  const Vector3d p0 = state1.head(3);
  const Vector3d dp = state2.head(3) - p0;
  const Vector3d v0 = state1.segment(3, 3);
  const Vector3d v1 = state2.segment(3, 3);
  const Vector3d dv = v1 - v0;
  double t_d = time_to_goal;
  MatrixXd coef(3, 4);
  end_vel_ = v1;

  Vector3d a = 1.0 / 6.0 * (-12.0 / (t_d * t_d * t_d) * (dp - v0 * t_d) + 6 / (t_d * t_d) * dv);
  Vector3d b = 0.5 * (6.0 / (t_d * t_d) * (dp - v0 * t_d) - 2 / t_d * dv);
  Vector3d c = v0;
  Vector3d d = p0;

  // 1/6 * alpha * t^3 + 1/2 * beta * t^2 + v0
  // a*t^3 + b*t^2 + v0*t + p0
  coef.col(3) = a, coef.col(2) = b, coef.col(1) = c, coef.col(0) = d;

  Vector3d coord, vel, acc;
  VectorXd poly1d, t, polyv, polya;
  Vector3i index;

  Eigen::MatrixXd Tm(4, 4);
  Tm << 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 3, 0, 0, 0, 0;

  /* ---------- forward checking of trajectory ---------- */
  double t_delta = t_d / 10;
  for (double time = t_delta; time <= t_d; time += t_delta)
  {
    t = VectorXd::Zero(4);
    for (int j = 0; j < 4; j++)
      t(j) = pow(time, j);

    for (int dim = 0; dim < 3; dim++)
    {
      poly1d = coef.row(dim);
      coord(dim) = poly1d.dot(t);
      vel(dim) = (Tm * poly1d).dot(t);
      acc(dim) = (Tm * Tm * poly1d).dot(t);

      if (fabs(vel(dim)) > max_vel_ || fabs(acc(dim)) > max_acc_)
      {
        // cout << "vel:" << vel(dim) << ", acc:" << acc(dim) << endl;
        // return false;
      }
    }

    if (coord(0) < origin_(0) || coord(0) >= map_size_3d_(0) || coord(1) < origin_(1) || coord(1) >= map_size_3d_(1) ||
        coord(2) < origin_(2) || coord(2) >= map_size_3d_(2))
    {
      return false;
    }

    // if (edt_environment_->evaluateCoarseEDT(coord, -1.0) <= margin_) {
    //   return false;
    // }
    if (edt_environment_->sdf_map_->getInflateOccupancy(coord) == 1)
    {
      return false;
    }
  }
  coef_shot_ = coef;
  t_shot_ = t_d;
  is_shot_succ_ = true;
  return true;
}

vector<double> KinodynamicAstar::cubic(double a, double b, double c, double d)
{
  vector<double> dts;

  double a2 = b / a;
  double a1 = c / a;
  double a0 = d / a;

  double Q = (3 * a1 - a2 * a2) / 9;
  double R = (9 * a1 * a2 - 27 * a0 - 2 * a2 * a2 * a2) / 54;
  double D = Q * Q * Q + R * R;
  if (D > 0)
  {
    double S = std::cbrt(R + sqrt(D));
    double T = std::cbrt(R - sqrt(D));
    dts.push_back(-a2 / 3 + (S + T));
    return dts;
  }
  else if (D == 0)
  {
    double S = std::cbrt(R);
    dts.push_back(-a2 / 3 + S + S);
    dts.push_back(-a2 / 3 - S);
    return dts;
  }
  else
  {
    double theta = acos(R / sqrt(-Q * Q * Q));
    dts.push_back(2 * sqrt(-Q) * cos(theta / 3) - a2 / 3);
    dts.push_back(2 * sqrt(-Q) * cos((theta + 2 * M_PI) / 3) - a2 / 3);
    dts.push_back(2 * sqrt(-Q) * cos((theta + 4 * M_PI) / 3) - a2 / 3);
    return dts;
  }
}

vector<double> KinodynamicAstar::quartic(double a, double b, double c, double d, double e)
{
  //四次方程的费拉里解法，https://zh.wikipedia.org/wiki/%E5%9B%9B%E6%AC%A1%E6%96%B9%E7%A8%8B，降次为三次方程电泳cubic
  // 四次方程无法直接求解，必须通过 “降次” 转化为三次方程（费拉里解法），三次方程再通过卡尔达诺公式或三角函数法求解。
  // 四次方程的求解核心是费拉里解法：通过引入一个新变量y，将四次方程分解为两个二次方程，
  //    而y的取值需要先解一个三次方程（这就是为什么quartic函数会调用cubic函数）。

  // 费拉里解法的核心思想是：将四次方程写成  的形式，展开后对比系数，发现p, q, r, s的求解依赖一个三次方程（即 “辅助三次方程”）。

  vector<double> dts;

  double a3 = b / a;
  double a2 = c / a;
  double a1 = d / a;
  double a0 = e / a;

  vector<double> ys = cubic(1, -a2, a1 * a3 - 4 * a0, 4 * a2 * a0 - a1 * a1 - a3 * a3 * a0);
  double y1 = ys.front();
  double r = a3 * a3 / 4 - a2 + y1;
  if (r < 0)
    return dts;

  double R = sqrt(r);
  double D, E;
  if (R != 0)
  {
    D = sqrt(0.75 * a3 * a3 - R * R - 2 * a2 + 0.25 * (4 * a3 * a2 - 8 * a1 - a3 * a3 * a3) / R);
    E = sqrt(0.75 * a3 * a3 - R * R - 2 * a2 - 0.25 * (4 * a3 * a2 - 8 * a1 - a3 * a3 * a3) / R);
  }
  else
  {
    D = sqrt(0.75 * a3 * a3 - 2 * a2 + 2 * sqrt(y1 * y1 - 4 * a0));
    E = sqrt(0.75 * a3 * a3 - 2 * a2 - 2 * sqrt(y1 * y1 - 4 * a0));
  }

  if (!std::isnan(D))
  {
    dts.push_back(-a3 / 4 + R / 2 + D / 2);
    dts.push_back(-a3 / 4 + R / 2 - D / 2);
  }
  if (!std::isnan(E))
  {
    dts.push_back(-a3 / 4 - R / 2 + E / 2);
    dts.push_back(-a3 / 4 - R / 2 - E / 2);
  }

  return dts;
}

void KinodynamicAstar::init()
{
  /* ---------- map params ---------- */
  this->inv_resolution_ = 1.0 / resolution_;
  inv_time_resolution_ = 1.0 / time_resolution_;
  edt_environment_->sdf_map_->getRegion(origin_, map_size_3d_);

  cout << "origin_: " << origin_.transpose() << endl;
  cout << "map size: " << map_size_3d_.transpose() << endl;

  /* ---------- pre-allocated node ---------- */
  path_node_pool_.resize(allocate_num_);
  for (int i = 0; i < allocate_num_; i++)
  {
    path_node_pool_[i] = new PathNode;
  }

  phi_ = Eigen::MatrixXd::Identity(6, 6);
  use_node_num_ = 0;
  iter_num_ = 0;
}

void KinodynamicAstar::setEnvironment(const EDTEnvironment::Ptr& env)
{
  this->edt_environment_ = env;
}

void KinodynamicAstar::reset()
{
  expanded_nodes_.clear();
  path_nodes_.clear();

  std::priority_queue<PathNodePtr, std::vector<PathNodePtr>, NodeComparator> empty_queue;
  open_set_.swap(empty_queue);

  for (int i = 0; i < use_node_num_; i++)
  {
    PathNodePtr node = path_node_pool_[i];
    node->parent = NULL;
    node->node_state = NOT_EXPAND;
  }

  use_node_num_ = 0;
  iter_num_ = 0;
  is_shot_succ_ = false;
  has_path_ = false;
}

std::vector<Eigen::Vector3d> KinodynamicAstar::getKinoTraj(double delta_t)
{
  vector<Vector3d> state_list;

  /* ---------- get traj of searching ---------- */
  PathNodePtr node = path_nodes_.back();
  Matrix<double, 6, 1> x0, xt;

  while (node->parent != NULL)
  {
    Vector3d ut = node->input;
    double duration = node->duration;
    x0 = node->parent->state;

    for (double t = duration; t >= -1e-5; t -= delta_t)
    {
      stateTransit(x0, xt, ut, t);
      state_list.push_back(xt.head(3));
    }
    node = node->parent;
  }
  reverse(state_list.begin(), state_list.end());
  /* ---------- get traj of one shot ---------- */
  if (is_shot_succ_)
  {
    Vector3d coord;
    VectorXd poly1d, time(4);

    for (double t = delta_t; t <= t_shot_; t += delta_t)
    {
      for (int j = 0; j < 4; j++)
        time(j) = pow(t, j);

      for (int dim = 0; dim < 3; dim++)
      {
        poly1d = coef_shot_.row(dim);
        coord(dim) = poly1d.dot(time);
      }
      state_list.push_back(coord);
    }
  }

  return state_list;
}

void KinodynamicAstar::getSamples(double& ts, vector<Eigen::Vector3d>& point_set,
                                  vector<Eigen::Vector3d>& start_end_derivatives)
{
  /* ---------- path duration ---------- */
  double T_sum = 0.0;
  if (is_shot_succ_)
    T_sum += t_shot_;
  PathNodePtr node = path_nodes_.back();
  while (node->parent != NULL)
  {
    T_sum += node->duration;
    node = node->parent;
  }
  // cout << "duration:" << T_sum << endl;

  // Calculate boundary vel and acc
  Eigen::Vector3d end_vel, end_acc;
  double t;
  if (is_shot_succ_)
  {
    t = t_shot_;
    end_vel = end_vel_;
    for (int dim = 0; dim < 3; ++dim)
    {
      Vector4d coe = coef_shot_.row(dim);
      end_acc(dim) = 2 * coe(2) + 6 * coe(3) * t_shot_;
    }
  }
  else
  {
    t = path_nodes_.back()->duration;
    end_vel = node->state.tail(3);
    end_acc = path_nodes_.back()->input;
  }

  // Get point samples
  int seg_num = floor(T_sum / ts);
  seg_num = max(8, seg_num);
  ts = T_sum / double(seg_num);
  bool sample_shot_traj = is_shot_succ_;
  node = path_nodes_.back();

  for (double ti = T_sum; ti > -1e-5; ti -= ts)
  {
    if (sample_shot_traj)
    {
      // samples on shot traj
      Vector3d coord;
      Vector4d poly1d, time;

      for (int j = 0; j < 4; j++)
        time(j) = pow(t, j);

      for (int dim = 0; dim < 3; dim++)
      {
        poly1d = coef_shot_.row(dim);
        coord(dim) = poly1d.dot(time);
      }

      point_set.push_back(coord);
      t -= ts;

      /* end of segment */
      if (t < -1e-5)
      {
        sample_shot_traj = false;
        if (node->parent != NULL)
          t += node->duration;
      }
    }
    else
    {
      // samples on searched traj
      Eigen::Matrix<double, 6, 1> x0 = node->parent->state;
      Eigen::Matrix<double, 6, 1> xt;
      Vector3d ut = node->input;

      stateTransit(x0, xt, ut, t);

      point_set.push_back(xt.head(3));
      t -= ts;

      // cout << "t: " << t << ", t acc: " << T_accumulate << endl;
      if (t < -1e-5 && node->parent->parent != NULL)
      {
        node = node->parent;
        t += node->duration;
      }
    }
  }
  reverse(point_set.begin(), point_set.end());

  // calculate start acc
  Eigen::Vector3d start_acc;
  if (path_nodes_.back()->parent == NULL)
  {
    // no searched traj, calculate by shot traj
    start_acc = 2 * coef_shot_.col(2);
  }
  else
  {
    // input of searched traj
    start_acc = node->input;
  }

  start_end_derivatives.push_back(start_vel_);
  start_end_derivatives.push_back(end_vel);
  start_end_derivatives.push_back(start_acc);
  start_end_derivatives.push_back(end_acc);
}

std::vector<PathNodePtr> KinodynamicAstar::getVisitedNodes()
{
  vector<PathNodePtr> visited;
  visited.assign(path_node_pool_.begin(), path_node_pool_.begin() + use_node_num_ - 1);
  return visited;
}

Eigen::Vector3i KinodynamicAstar::posToIndex(Eigen::Vector3d pt)
{
  Vector3i idx = ((pt - origin_) * inv_resolution_).array().floor().cast<int>();

  // idx << floor((pt(0) - origin_(0)) * inv_resolution_), floor((pt(1) -
  // origin_(1)) * inv_resolution_),
  //     floor((pt(2) - origin_(2)) * inv_resolution_);

  return idx;
}

int KinodynamicAstar::timeToIndex(double time)
{
  int idx = floor((time - time_origin_) * inv_time_resolution_);
  return idx;
}

void KinodynamicAstar::stateTransit(Eigen::Matrix<double, 6, 1>& state0, Eigen::Matrix<double, 6, 1>& state1,
                                    Eigen::Vector3d um, double tau)
{
  for (int i = 0; i < 3; ++i)
    phi_(i, i + 3) = tau;

  Eigen::Matrix<double, 6, 1> integral;
  integral.head(3) = 0.5 * pow(tau, 2) * um;
  integral.tail(3) = tau * um;

  state1 = phi_ * state0 + integral;
}

}  // namespace fast_planner
